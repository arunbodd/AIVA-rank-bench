#!/usr/bin/env python3
"""
benchmark_gene_ranking.py — Gene Prioritization Benchmark (Full-Information Mode)

Single ablation: hpo+vcf+genes+disease (sends all available information to AIVA).
LLM assigns scores + ranks. Computes extended metrics: Recall@50, Retention Rate, etc.

Usage:
  python benchmark_gene_ranking.py --limit 100
  python benchmark_gene_ranking.py --evaluate-only ranked_predictions.jsonl
"""

import sys
import os
import json
import math
import re
import asyncio
import argparse
import httpx
from pathlib import Path
from datetime import datetime, timezone
from collections import defaultdict

# Strip trailing "(ENSGXXX)" parentheticals from AIVA gene output
_ID_IN_PARENS = re.compile(r'^(.+?)\s*\([^)]+\)\s*$')

# Import reusable helpers from predict_genes.py
sys.path.insert(0, str(Path(__file__).parent))
from predict_genes import (
    load_disease_id_map,
    EnsemblResolver,
    load_patients,
    fmt_hpo,
    _variant_gene_symbol,
    call_aiva,
)

# ──────────────────────────────────────────────────────────────
# Paths and constants
# ──────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
# Configurable via env vars; defaults point at repo-bundled data.
DEFAULT_INPUT = Path(os.getenv("AIVA_BENCHMARK_INPUTS", REPO_ROOT / "inputs" / "benchmark_inputs.jsonl"))
DISEASE_LIST = Path(os.getenv("AIVA_DISEASE_LIST", REPO_ROOT / "inputs" / "simulated_diseases_list.txt"))
DEFAULT_RESULTS = REPO_ROOT / "results"

HPO_LABELS_CACHE = DEFAULT_RESULTS / "hp_labels.json"
# HPO OBO file download (we'll parse it to extract ID→label mapping)
HPO_OBO_URL = "https://purl.obolibrary.org/obo/hp.obo"

# ──────────────────────────────────────────────────────────────
# HPO Auto-Download
# ──────────────────────────────────────────────────────────────
def ensure_hpo_labels(cache_path: Path = HPO_LABELS_CACHE) -> dict[str, str]:
    """Return {HP:XXXXXXX: label} — download from OBO file once and cache."""
    if cache_path.exists():
        with open(cache_path, encoding="utf-8") as f:
            return json.load(f)

    import urllib.request
    print(f"Downloading HPO OBO file → {cache_path} …", file=sys.stderr)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        with urllib.request.urlopen(HPO_OBO_URL, timeout=120) as r:
            obo_text = r.read().decode("utf-8")

        # Parse OBO format: extract [Term] blocks with id: and name:
        mapping = {}
        lines = obo_text.split("\n")
        current_id = None
        for line in lines:
            line = line.strip()
            if line.startswith("id: HP:"):
                current_id = line.split(":", 1)[1].strip()
            elif line.startswith("name:") and current_id:
                name = line.split(":", 1)[1].strip()
                mapping[current_id] = name
                current_id = None

        with open(cache_path, "w", encoding="utf-8") as f:
            json.dump(mapping, f, indent=2)
        print(f"  Cached {len(mapping):,} HPO terms.", file=sys.stderr)
        return mapping
    except Exception as exc:
        print(f"Warning: HPO download failed: {exc}", file=sys.stderr)
        return {}


# ──────────────────────────────────────────────────────────────
# JSON schema and prompt templates
# ──────────────────────────────────────────────────────────────
_JSON_SCHEMA = (
    'Return ONLY a valid JSON object — no markdown, no explanation — in this exact format:\n'
    '{"predicted_genes": [{"gene": "GENE1", "score": 0.95}, {"gene": "GENE2", "score": 0.88}, ...]}\n'
    "List genes in order from most to least likely causal (rank = array order). "
    "Score each gene 0.0–1.0 (1.0 = most likely causal). "
    "Include all genes from the Candidate Genes list; put non-causal ones last with low scores."
)

_AIVA_TASK = (
    "\nGiven the above patient information, use appropriate tools to identify and rank likely genes based on the genetic variants provided.\n"
    "Prioritize literature search and don't web search unless absolutely necessary.\n"
)

_GENERIC_TASK = (
    "Identify and rank the most likely causal gene(s) for the patient described above."
)

PROMPT_STYLES = {
    "aiva":    _AIVA_TASK,
    "generic": _GENERIC_TASK,
}


def _age_block(age: str) -> str:
    if age and age.lower() not in ("unknown", ""):
        return f"\nAge of Onset: {age.replace('_', ' ')}\n"
    return ""


def _hpo_block(positive: list[str], negative: list[str], hpo_labels: dict) -> str:
    lines = []
    if positive:
        labels = [fmt_hpo(t, hpo_labels) for t in positive]
        lines.append(f"\nPositive Phenotypes:\n  " + "\n  ".join(labels))
    if negative:
        labels = [fmt_hpo(t, hpo_labels) for t in negative]
        lines.append(f"\nNegative Phenotypes:\n  " + "\n  ".join(labels))
    return "\n".join(lines)


def _fmt_variant_annot(v: dict) -> str:
    """Compact annotation string for variant."""
    parts = []
    csq = v.get("csq", "")
    if csq:
        parts.append(f"csq={csq}")
    cadd = v.get("cadd")
    if cadd is not None:
        parts.append(f"CADD={cadd:.1f}")
    pp2 = v.get("polyphen2")
    if pp2 is not None:
        parts.append(f"PP2={pp2:.3f}")
    if v.get("maf_absent"):
        parts.append("AF=rare")
    elif v.get("maf") is not None:
        parts.append(f"AF={v['maf']:.2e}")
    biotype = v.get("biotype", "")
    if biotype and biotype != "protein_coding":
        parts.append(f"bt={biotype}")
    return f" [{', '.join(parts)}]" if parts else ""


def _vcf_block_named(variants: list[dict]) -> str:
    """Variant lines WITH gene names and coordinates visible."""
    if not variants:
        return ""
    lines = ["\nCandidate Genes (with variant coordinates):"]
    for v in variants[:50]:  # show up to 50 genes
        gene = _variant_gene_symbol(v)
        chrom = v.get("chrom", "?")
        pos = v.get("pos", "?")
        ref = v.get("ref", "?")
        alt = v.get("alt", "?")
        coord = f"chr{chrom}:{pos} {ref}>{alt}"
        annot = _fmt_variant_annot(v)
        # _fmt_variant_annot returns " [csq=..., ...]" with brackets, so just concatenate
        lines.append(f"  {gene}: {coord}{annot}")
    if len(variants) > 50:
        lines.append(f"  ... and {len(variants) - 50} more")
    return "\n".join(lines)


def build_prompt_full_info(patient: dict, hpo_labels: dict, task_prompt: str = _AIVA_TASK, include_disease_id: bool = True) -> str:
    """Build prompt for hpo+vcf+genes+disease ablation (full information)."""
    parts = [_age_block(patient.get("age", ""))]

    # Disease ID (suspected diagnosis) - optional for comparison studies
    if include_disease_id:
        disease_id = patient.get("disease_id", "")
        disease_name = patient.get("disease_name", f"Disease_{disease_id}")
        if disease_id:
            parts.append(f"\nSuspected disease: {disease_name} (ORPHA:{disease_id})\n")

    # HPO terms (all tiers)
    parts.append(_hpo_block(patient.get("positive_hpo", []), patient.get("negative_hpo", []), hpo_labels))

    # Variants with gene names
    parts.append(_vcf_block_named(patient.get("variants", [])))

    # Task instruction + schema
    post_instruction = (
        f"\n{task_prompt}\n\n"
        "If gene symbols are listed in the `Candidate Genes` section, you MUST include every one in the returned `predicted_genes` array.\n"
        "If a gene is not supported as causal, list it last with a low score. Do not omit any gene from the list.\n\n"
        + _JSON_SCHEMA
    )
    parts.append(post_instruction)

    return "".join(parts)


# ──────────────────────────────────────────────────────────────
# Response parsing (extracts ranked list + LLM-assigned scores)
# ──────────────────────────────────────────────────────────────
def parse_predictions(response_text: str) -> tuple[list[str], dict[str, float]]:
    """Parse ranked+scored gene list. Rank = array order; score = LLM confidence."""
    genes: list[str] = []
    scores: dict[str, float] = {}

    decoder = json.JSONDecoder()
    idx = 0
    while idx < len(response_text):
        if response_text[idx] != "{":
            idx += 1
            continue
        try:
            parsed, span = decoder.raw_decode(response_text[idx:])
        except json.JSONDecodeError:
            idx += 1
            continue

        if isinstance(parsed, dict) and "predicted_genes" in parsed:
            pg = parsed["predicted_genes"]
            if isinstance(pg, list):
                for entry in pg:
                    if isinstance(entry, dict):
                        gene = entry.get("gene") or entry.get("symbol") or entry.get("name")
                        scr = entry.get("score")
                        if gene:
                            genes.append(str(gene))
                            if scr is not None:
                                try:
                                    scores[str(gene).upper()] = float(scr)
                                except (ValueError, TypeError):
                                    pass
                    elif isinstance(entry, str):
                        # graceful fallback if model returns plain strings
                        genes.append(entry)
            break
        idx += max(span, 1)

    # De-duplicate preserving array order (= rank); array order is authoritative
    seen: set[str] = set()
    deduped: list[str] = []
    for g in genes:
        if g.upper() not in seen:
            seen.add(g.upper())
            deduped.append(g)

    return deduped, scores


# ──────────────────────────────────────────────────────────────
# Evaluation metrics
# ──────────────────────────────────────────────────────────────
def rank_of_true_gene(predicted: list[str], true_genes: list[str]) -> int | None:
    """Return 1-based rank of the first true gene in the predicted list, or None.

    Normalises predicted strings so "TEK (ENSG00000120156)" matches "TEK".
    """
    true_upper = {g.upper() for g in true_genes}
    for i, g in enumerate(predicted, 1):
        m = _ID_IN_PARENS.match(g.strip())
        g_upper = (m.group(1).strip() if m else g.strip()).upper()
        if g_upper in true_upper:
            return i
    return None


def compute_metrics(results: list[dict], at_ks: tuple[int, ...] = (1, 5, 10, 20, 50)) -> dict:
    """Compute Recall@K, MRR, NDCG@10, Retention Rate, Filtered-out Rate."""
    found_at = {k: 0 for k in at_ks}
    retention_count = 0
    ranks: list[int] = []
    reciprocal_ranks: list[float] = []
    dcg_list: list[float] = []
    no_prediction = 0

    per_category: dict[str, list[int | None]] = defaultdict(list)

    for r in results:
        predicted = r.get("predicted", [])
        raw = r.get("predicted_before_force", predicted)
        true_genes = r.get("true_genes", [])
        cat = r.get("category", "unknown")

        if not predicted:
            no_prediction += 1
            per_category[cat].append(None)
            reciprocal_ranks.append(0.0)
            dcg_list.append(0.0)
            continue

        # Retention = true gene in raw model output (before force-inclusion)
        if rank_of_true_gene(raw, true_genes) is not None:
            retention_count += 1

        # Rank in final list (after force-inclusion)
        rk = rank_of_true_gene(predicted, true_genes)
        per_category[cat].append(rk)

        if rk is not None:
            ranks.append(rk)
            reciprocal_ranks.append(1.0 / rk)
            dcg_list.append(1.0 / math.log2(rk + 1))
            for k in at_ks:
                if rk <= k:
                    found_at[k] += 1
        else:
            reciprocal_ranks.append(0.0)
            dcg_list.append(0.0)

    n = len(results)
    if n == 0:
        return {}

    idcg = 1.0  # ideal DCG when true gene is at rank 1

    metrics: dict = {
        "n_patients": n,
        "no_prediction": no_prediction,
        "mrr": round(sum(reciprocal_ranks) / n, 4),
        "ndcg@10": round(sum(dcg_list) / (n * idcg), 4),
        "mean_rank": round(sum(ranks) / len(ranks), 2) if ranks else None,
        "median_rank": float(sorted(ranks)[len(ranks) // 2]) if ranks else None,
        "retention_rate": round(retention_count / n, 4),
        "filtered_out_rate": round(1.0 - retention_count / n, 4),
    }
    for k in at_ks:
        metrics[f"recall@{k}"] = round(found_at[k] / n, 4)

    # Per-category recall@K
    cat_metrics: dict[str, dict] = {}
    for cat, cat_ranks in per_category.items():
        n_cat = len(cat_ranks)
        cm: dict = {"n": n_cat}
        for k in at_ks:
            found_k = sum(1 for rk in cat_ranks if rk is not None and rk <= k)
            cm[f"recall@{k}"] = round(found_k / n_cat, 4)
        cat_metrics[cat] = cm
    metrics["per_category"] = cat_metrics

    return metrics


def print_metrics(metrics: dict):
    print(f"\n{'═' * 60}")
    print(f"  BENCHMARK RESULTS  —  hpo+vcf+genes+disease")
    print(f"{'═' * 60}")
    print(f"  Patients evaluated  : {metrics.get('n_patients')}")
    print(f"  No prediction       : {metrics.get('no_prediction')}")
    print(f"  MRR                 : {metrics.get('mrr')}")
    print(f"  NDCG@10             : {metrics.get('ndcg@10')}")
    print(f"  Mean rank           : {metrics.get('mean_rank')}")
    print(f"  Median rank         : {metrics.get('median_rank')}")
    print(f"  Retention rate      : {metrics.get('retention_rate'):.1%}")
    print(f"  Filtered-out rate   : {metrics.get('filtered_out_rate'):.1%}")
    for k in (1, 5, 10, 20, 50):
        key = f"recall@{k}"
        if key in metrics:
            print(f"  Recall@{k:<2}           : {metrics[key]:.1%}")
    print()

    cat_m = metrics.get("per_category", {})
    if cat_m:
        print(f"  Per-category Recall@5:")
        for cat, cm in sorted(cat_m.items()):
            print(f"    {cat:50s} : {cm.get('recall@5', 0):.1%} (n={cm.get('n', 0)})")
    print()


# ──────────────────────────────────────────────────────────────
# Patient processing loop
# ──────────────────────────────────────────────────────────────
async def process_patients(
    patients: list[dict],
    output_path: Path,
    hpo_labels: dict,
    dry_run: bool = False,
    task_prompt: str = _AIVA_TASK,
    save_debug: bool = False,
    include_disease_id: bool = True,
) -> list[dict]:
    """Process patients and write results to output_path."""
    results = []
    output_path.parent.mkdir(parents=True, exist_ok=True)

    async with httpx.AsyncClient() as client:
        with open(output_path, "w", encoding="utf-8") as out_f:
            for i, p in enumerate(patients, 1):
                pid = p["patient_id"]
                print(f"[{i}/{len(patients)}] Processing patient {pid} ...", end=" ", flush=True)

                # Build prompt
                prompt = build_prompt_full_info(p, hpo_labels, task_prompt, include_disease_id)

                if dry_run:
                    print("\n" + "─" * 60)
                    print(prompt)
                    print("─" * 60)
                    continue

                # Call AIVA
                try:
                    response_text = await call_aiva(client, prompt)
                except Exception as exc:
                    print(f"API error: {exc}")
                    response_text = ""
                # Normalize None -> empty string: call_aiva may return None on errors/status
                if response_text is None:
                    response_text = ""

                # Parse predictions
                predicted, scores = parse_predictions(response_text)
                predicted_before_force = list(predicted)

                # Force-inclusion: any variant gene not returned by model is appended
                variant_genes = []
                seen_vg = set()
                for v in p.get("variants", []):
                    gene = _variant_gene_symbol(v)
                    if gene and gene != "Unknown" and gene.upper() not in seen_vg:
                        seen_vg.add(gene.upper())
                        variant_genes.append(gene)

                predicted_upper = {g.upper() for g in predicted}
                for g in variant_genes:
                    if g.upper() not in predicted_upper:
                        predicted.append(g)
                        scores[g.upper()] = 0.0  # mark as model-absent

                # Rank of true gene
                rank = rank_of_true_gene(predicted, p["true_genes"])

                # Build output record
                rec = {
                    "patient_id": pid,
                    "disease_id": p["disease_id"],
                    "disease_name": p.get("disease_name", f"Disease_{p['disease_id']}"),
                    "true_genes": p["true_genes"],
                    "predicted": predicted,
                    "predicted_scores": {g: scores.get(g.upper(), 0.0) for g in predicted},
                    "predicted_before_force": predicted_before_force,
                    "rank_of_true_gene": rank,
                    "category": p.get("category", ""),
                    "in_gene_gene_kg": p.get("in_gene_gene_kg"),
                    "in_disease_gene_kg": p.get("in_disease_gene_kg"),
                    "timestamp": datetime.now(timezone.utc).isoformat(),
                }
                if save_debug:
                    rec["prompt"] = prompt
                    rec["raw_response"] = response_text or ""

                out_f.write(json.dumps(rec) + "\n")
                results.append(rec)

                print(f"rank={rank}")

        return results


# ──────────────────────────────────────────────────────────────
# Evaluate-only mode (load existing results and compute metrics)
# ──────────────────────────────────────────────────────────────
def evaluate_only(results_path: Path):
    """Load existing results and compute metrics."""
    print(f"Loading results from: {results_path}")
    results = []
    with open(results_path, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                results.append(json.loads(line))

    print(f"Loaded {len(results)} patient results.")
    metrics = compute_metrics(results)
    print_metrics(metrics)

    # Save evaluation summary
    summary_path = results_path.parent / (results_path.stem + "_evaluation_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({"metrics": metrics, "timestamp": datetime.now(timezone.utc).isoformat()}, f, indent=2)
    print(f"Saved evaluation summary to: {summary_path}")


# ──────────────────────────────────────────────────────────────
# CLI and main
# ──────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(description="Gene prioritization benchmark (full-info mode)")
    p.add_argument("--input", default=str(DEFAULT_INPUT),
                   help="Path to benchmark inputs JSONL (default: inputs/benchmark_inputs.jsonl)")
    p.add_argument("--output", default=str(DEFAULT_RESULTS / "ranked_predictions.jsonl"),
                   help="Output file for predictions (default: Results/ranked_predictions.jsonl)")
    p.add_argument("--hpo-labels", default=None,
                   help="Optional path to HPO label JSON. Auto-downloads from JAX if omitted.")
    p.add_argument("--prompt-style", choices=list(PROMPT_STYLES), default="aiva",
                   help="'aiva' = directs AIVA to query KG (default); 'generic' = minimal instruction.")
    p.add_argument("--limit", type=int, default=0,
                   help="Cap the number of patients to process (0 = all).")
    p.add_argument("--match-patients", metavar="FILE", default=None,
                   help="Run on the same patients (by ID) as in this existing predictions JSONL file. Ensures fair comparison.")
    p.add_argument("--no-disease-id", action="store_true",
                   help="Exclude disease ID from prompts (for comparison study: with vs without disease hint).")
    p.add_argument("--dry-run", action="store_true",
                   help="Print prompts without calling the API.")
    p.add_argument("--evaluate-only", metavar="FILE", default=None,
                   help="Skip inference; compute metrics from an existing results file.")
    p.add_argument("--save-debug", action="store_true",
                   help="Include 'prompt' and 'raw_response' fields in each output record.")
    return p.parse_args()


def main():
    args = parse_args()

    # ── Evaluate-only shortcut ──
    if args.evaluate_only:
        path = Path(args.evaluate_only)
        if not path.exists():
            sys.exit(f"File not found: {path}")
        evaluate_only(path)
        return

    # ── Load patients ──
    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Input file not found: {input_path}")

    id_map = load_disease_id_map(DISEASE_LIST)
    resolver = EnsemblResolver()
    hpo_labels = ensure_hpo_labels(Path(args.hpo_labels)) if args.hpo_labels else ensure_hpo_labels()

    print(f"Loading patients from: {input_path}")
    patients_all = load_patients(input_path, id_map, resolver, limit=0)

    # ── Filter by --match-patients if provided ──
    if args.match_patients:
        match_path = Path(args.match_patients)
        if not match_path.exists():
            sys.exit(f"--match-patients file not found: {match_path}")

        print(f"Matching patients from: {match_path}")
        target_ids = []
        with open(match_path, encoding="utf-8") as f:
            for line in f:
                if line.strip():
                    rec = json.loads(line)
                    target_ids.append(rec["patient_id"])

        # Build lookup dict for fast access
        patients_by_id = {p["patient_id"]: p for p in patients_all}
        patients = []
        for pid in target_ids:
            if pid in patients_by_id:
                patients.append(patients_by_id[pid])
            else:
                print(f"Warning: patient {pid} from {match_path} not found in input data.", file=sys.stderr)

        print(f"Matched {len(patients)}/{len(target_ids)} patients from comparison file.")
    elif args.limit > 0:
        patients = patients_all[:args.limit]
    else:
        patients = patients_all

    print(f"Processing {len(patients)} patients (total available: {len(patients_all)})")

    # ── Build prompts for all patients (always send all HPO terms + variants + disease_id) ──
    for p in patients:
        p["positive_hpo"] = sorted(p.get("positive_hpo", []))
        p["negative_hpo"] = sorted(p.get("negative_hpo", []))

    # ── Run predictions ──
    output_path = Path(args.output)
    task_prompt = PROMPT_STYLES[args.prompt_style]

    print(f"Output will be written to: {output_path}")

    results = asyncio.run(process_patients(
        patients,
        output_path,
        hpo_labels,
        dry_run=args.dry_run,
        task_prompt=task_prompt,
        save_debug=args.save_debug,
        include_disease_id=not args.no_disease_id,
    ))

    if args.dry_run:
        return

    # ── Evaluate ──
    metrics = compute_metrics(results)
    print_metrics(metrics)

    # Save evaluation summary
    summary_path = output_path.parent / (output_path.stem + "_evaluation_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({"metrics": metrics, "timestamp": datetime.now(timezone.utc).isoformat()}, f, indent=2)
    print(f"Saved evaluation summary to: {summary_path}")


if __name__ == "__main__":
    main()
