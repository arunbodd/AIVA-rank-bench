#!/usr/bin/env python3
"""
AIVA Gene Prioritisation Benchmark
====================================
Sends patient data to AIVA and evaluates gene ranking quality.

Ablation conditions (--ablation):
  hpo              : HPO terms + age only                         [floor: phenotype-only]
  hpo+vcf          : HPO + blinded variants (no gene names)       [PRIMARY benchmark]
  hpo+vcf+disease  : HPO + blinded variants + disease name        [disease-label uplift]
  hpo+vcf+genes    : HPO + variants WITH gene names (no disease)  [closed-list ranking]

"Blinded" means variant lines show only genomic coordinates + pathogenicity annotations;
gene symbols are withheld so AIVA must resolve gene identity from its own internal databases.
In hpo+vcf+genes each variant line is paired with its HGNC gene symbol, letting AIVA rank
the ~10 candidate genes purely by phenotype-gene fit without knowing the disease.
"Candidate Genes" from the benchmark construction are NEVER disclosed as a flat list.

Reads from:
  data/final_extracted.jsonl        (default, main benchmark set)
  data/generated_vcfs/simulated/    (per-patient VCFs with CADD/PolyPhen2/MAF)

Writes to:
  Results/<ablation>_predictions.jsonl        (raw per-patient predictions)
  Results/<ablation>_evaluation_summary.json  (aggregate metrics)

Usage:
  # HPO-only ablation, first 20 patients for smoke-test
  python scripts/predict_genes.py --ablation hpo --limit 20 --dry-run

  # Full HPO+VCF benchmark run
  python scripts/predict_genes.py --ablation hpo+vcf

  # Evaluate only (no API calls) — score a predictions file
  python scripts/predict_genes.py --evaluate-only Results/hpo+vcf_predictions.jsonl

Environment:
  AIVA_API_TOKEN   API bearer token (required)
  AIVA_API_URL     Base URL or full /agents/chat/stream path (optional)
  AIVA_MODEL       Model override (optional)
"""

import argparse
import asyncio
import csv
import json
import math
import os
import random
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Normalise predicted gene strings: strip trailing "(ENSGXXX)" or similar parentheticals
# that AIVA sometimes appends, e.g. "TEK (ENSG00000120156)" → "TEK"
_ID_IN_PARENS = re.compile(r'^(.+?)\s*\([^)]+\)\s*$')

HGNC_PATH = Path(__file__).parents[4] / "Lirical_benchmark" / "tools" / \
    "lirical-cli-2.2.1" / "data" / "hgnc_complete_set.txt"

try:
    import httpx
except ImportError:
    sys.exit("httpx not installed. Run: pip install httpx")

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    pass

# ──────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────
RAW_API_URL = os.getenv("AIVA_API_URL", "https://api.aivaportal.com/agents/chat/stream")
API_KEY     = os.getenv("AIVA_API_KEY") or os.getenv("AIVA_API_TOKEN", "")
MODEL       = os.getenv("AIVA_MODEL") or os.getenv("DEFAULT_MODEL", "google/gemini-3.1-pro-preview")

API_URL = RAW_API_URL.rstrip("/")
if "/agents/chat/stream" not in API_URL:
    API_URL = f"{API_URL}/agents/chat/stream"

SCRIPT_DIR   = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent                          # Simulated_Datasets/
BENCH_ROOT   = PROJECT_ROOT.parent                        # AIVA-var-bench/
DEFAULT_INPUT    = PROJECT_ROOT / "data" / "final_extracted.jsonl"
DEFAULT_RESULTS  = PROJECT_ROOT / "Results"
VCF_DIR          = PROJECT_ROOT / "data" / "generated_vcfs" / "simulated"
DISEASE_LIST     = PROJECT_ROOT / "data" / "simulated_diseases_list.txt"

REQUEST_TIMEOUT  = 120.0
CONCURRENCY      = int(os.getenv("CONCURRENCY_LIMIT", "3"))  # max simultaneous API calls

VALID_ABLATIONS  = ("hpo", "hpo+vcf", "hpo+vcf+disease", "hpo+vcf+genes", "hpo+vcf+genes+disease")

# Columns whose source tag begins with "init_phenotypes" = definitive diagnostic HPOs.
# We expose ALL sources to AIVA (as Positive / Negative) because the benchmark intentionally
# includes distractors and noise to simulate real clinical uncertainty.


# ──────────────────────────────────────────────────────────────
# VCF annotation loader  (CADD, PolyPhen2, MAF, CSQ, biotype)
# ──────────────────────────────────────────────────────────────
def load_vcf_annotations(patient_id: int) -> dict:
    """Return {(chrom, pos): {gene_symbol, ensembl_id, cadd, polyphen2,
                               maf, maf_absent, csq, biotype}}."""
    vcf_path = VCF_DIR / f"simulated_{patient_id}_hg19.vcf"
    annots: dict = {}
    if not vcf_path.exists():
        return annots

    with open(vcf_path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos_str = cols[0], cols[1]
            info: dict = {}
            for field in cols[7].split(";"):
                if "=" in field:
                    k, v = field.split("=", 1)
                    info[k] = v

            maf_raw  = float(info["MAF"]) if "MAF" in info else None
            pp2_raw  = float(info["PolyPhen2"]) if "PolyPhen2" in info else None

            # Some VCF records (e.g. mitochondrial genes) have the fields inverted:
            # ensembl_id contains the HGNC symbol and name contains the Ensembl ID.
            # Detect this by checking which field looks like an Ensembl accession.
            raw_name = info.get("name", "")
            raw_eid  = info.get("ensembl_id", "")
            if raw_name.startswith("ENSG") and not raw_eid.startswith("ENSG"):
                # Swap: name is the Ensembl ID, ensembl_id is the HGNC symbol
                gene_symbol = raw_eid
                ensembl_id  = raw_name
            else:
                gene_symbol = raw_name
                ensembl_id  = raw_eid

            annots[(chrom, int(pos_str))] = {
                "gene_symbol":  gene_symbol,
                "ensembl_id":   ensembl_id,
                "cadd":         float(info["CADD"]) if "CADD" in info else None,
                "polyphen2":    pp2_raw,
                "maf":          None if (maf_raw is None or maf_raw < 0) else maf_raw,
                "maf_absent":   maf_raw is not None and maf_raw < 0,
                "csq":          info.get("CSQ", ""),
                "biotype":      info.get("biotype", ""),
            }
    return annots


def build_ensembl_symbol_map(annots: dict) -> dict[str, str]:
    """Return {ensembl_id: gene_symbol} from a patient's VCF annotations."""
    return {
        v["ensembl_id"]: v["gene_symbol"]
        for v in annots.values()
        if v["ensembl_id"] and v["gene_symbol"]
    }


# ──────────────────────────────────────────────────────────────
# HPO label loader  (optional; loads obo-style flat map)
# ──────────────────────────────────────────────────────────────
def load_hpo_labels(path: Path | None) -> dict[str, str]:
    """Load a simple TSV/JSON mapping HP:XXXXXXX → label.

    Accepts:
      - JSON: {"HP:0001234": "Hypotonia", ...}
      - TSV: HP:0001234\\tHypotonia (no header required)

    Returns empty dict if path is None or file not found.
    """
    if not path or not path.exists():
        return {}
    try:
        if path.suffix == ".json":
            with open(path, encoding="utf-8") as f:
                return json.load(f)
        else:
            mapping: dict[str, str] = {}
            with open(path, encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t", 1)
                    if len(parts) == 2 and parts[0].startswith("HP:"):
                        mapping[parts[0]] = parts[1]
            return mapping
    except Exception as exc:
        print(f"Warning: could not load HPO labels from {path}: {exc}", file=sys.stderr)
        return {}


def fmt_hpo(term: str, hpo_labels: dict[str, str]) -> str:
    label = hpo_labels.get(term)
    return f"{term} ({label})" if label else term


# ──────────────────────────────────────────────────────────────
# Prompt builders — one per ablation
# ──────────────────────────────────────────────────────────────
# Human-readable annotations written into every evaluation summary JSON
# ──────────────────────────────────────────────────────────────
SUMMARY_METRICS_EXPLANATION = {
    "n_patients":    "Total number of patients evaluated in this run.",
    "no_prediction": "Patients where AIVA returned no parseable gene list (API failure or empty response). These count as missed for all Recall@K.",
    "retention_rate": "Fraction of patients where the model actually returned the true gene (BEFORE force-inclusion). This is the true model performance.",
    "filtered_out_rate": "Fraction of patients where the model missed the true gene and it was only added via force-inclusion.",

    "recall@1":      "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene was the model's TOP-ranked gene. The strictest measure.",
    "recall@5":      "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 5.",
    "recall@10":     "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 10.",
    "recall@20":     "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 20.",
    "recall@50":     "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 50.",

    "recall_after_force@1":  "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@5":  "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@10": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@20": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@50": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",

    "mrr":           "Mean Reciprocal Rank (only counts genes with score > 0.0). Genes force-included with score 0.0 are treated as rank = None. 1.0 = always rank 1; 0.5 = average rank 2.",
    "ndcg@10":       "Normalised Discounted Cumulative Gain at 10 (only counts genes with score > 0.0). Genes force-included with score 0.0 are treated as rank = None. 1.0 = perfect.",
    "mean_rank":     "Average position among genes with score > 0.0 (model-assigned ranks only). Ignores arbitrary rankings of force-included genes.",
    "median_rank":   "Median position among genes with score > 0.0 (model-assigned ranks only). Ignores arbitrary rankings of force-included genes.",
    "mean_rank_before_force": "Average position of the true gene in the model's ACTUAL output (before force-inclusion). Only includes cases where model found it.",
    "median_rank_before_force": "Median position of the true gene in the model's ACTUAL output (before force-inclusion). Only includes cases where model found it.",
}

SUMMARY_CATEGORY_DESCRIPTIONS = {
    "known_gene_disease_gene_in_kg":            "ESTABLISHED CASE (baseline). Both the causal gene and its link to this disease are already documented in AIVA's knowledge graph. AIVA has direct evidence to draw on.",
    "known_gene_known_disease_gene_in_kg":       "FULLY DOCUMENTED CASE. The gene and disease-gene association are both in the KG, and the gene is additionally connected via gene-gene relationships. Highest information density.",
    "known_gene_new_disease_gene_in_kg":         "NOVEL DISEASE ASSOCIATION. The causal gene is well-known, but this specific gene-disease pairing is new — not yet in the disease-gene KG. Tests whether AIVA can reason from gene biology alone.",
    "new_gene_known_disease_gene_in_kg":         "NOVEL GENE WITH KG SUPPORT. The causal gene is novel for this disease but does appear in the KG via some connection. Intermediate difficulty.",
    "new_gene_known_disease_gene_only_in_gg_kg": "NOVEL GENE (hardest). The causal gene is absent from AIVA's disease-gene KG and only reachable via gene-gene connections. Tests reasoning from indirect evidence and phenotype alone.",
}

# ──────────────────────────────────────────────────────────────
_JSON_SCHEMA = (
    'Return ONLY a valid JSON object — no markdown, no explanation — in this exact format:\n'
    '{"predicted_genes": [{"gene": "GENE1", "score": 0.95}, {"gene": "GENE2", "score": 0.88}, ...]}\n'
    "List genes in order from most to least likely causal (rank = array order). "
    "Score each gene 0.0–1.0 (1.0 = most likely causal). "
    "Include all genes from the Candidate Genes list; put non-causal ones last with low scores."
)


def build_prompt_full_info(patient: dict, hpo_labels: dict, task_prompt: str | None = None, include_disease_id: bool = True) -> str:
    """Build prompt for full-information ablation (hpo+vcf+genes+disease).

    Uses HPO labels to render human-readable phenotype names and lists
    candidate genes with coordinates and compact annotations.
    """
    # Allow callers to omit task_prompt; default to module-level _AIVA_TASK
    if task_prompt is None:
        task_prompt = _AIVA_TASK

    parts = [_age_block(patient.get("age", ""))]

    if include_disease_id:
        disease_id = patient.get("disease_id", "")
        disease_name = patient.get("disease_name", f"Disease_{disease_id}")
        if disease_id:
            parts.append(f"\nSuspected disease: {disease_name} (ORPHA:{disease_id})\n")

    parts.append(_hpo_block(patient.get("positive_hpo", []), patient.get("negative_hpo", []), hpo_labels))
    parts.append(_vcf_block_named(patient.get("variants", [])))

    post_instruction = (
        f"\n{task_prompt}\n\n"
        "If gene symbols are listed in the `Candidate Genes` section, you MUST include every one in the returned `predicted_genes` array.\n"
        "If a gene is not supported as causal, list it last with a low score. Do not omit any gene from the list.\n\n"
        + _JSON_SCHEMA
    )
    parts.append(post_instruction)
    return "".join(parts)

_AIVA_TASK = (
    "\n\nGiven the above patient information, use appropriate tools to identify and rank highly likely genes based on the genetic variants provided.\n"
    "Prioritize literature search and don't web search unless absolutely necessary\n"
)

_GENERIC_TASK = (
    "Identify and rank the most likely causal gene(s) for the patient described above."
)

PROMPT_STYLES = {
    "aiva":    _AIVA_TASK,    # directs AIVA to use internal DBs / KG
    "generic": _GENERIC_TASK, # minimal instruction; model decides how to reason
}


def _age_block(age: str) -> str:
    if age and age.lower() not in ("unknown", ""):
        return f"\nAge of Onset: {age.replace('_', ' ')}\n"
    return ""


def _hpo_block(positive: list[str], negative: list[str], hpo_labels: dict) -> str:
    lines = []
    if positive:
        labels = [fmt_hpo(t, hpo_labels) for t in positive]
        lines.append(f"\nPositive Phenotypes (present):\n  " + "\n  ".join(labels))
    if negative:
        labels = [fmt_hpo(t, hpo_labels) for t in negative]
        lines.append(f"\nNegative Phenotypes (absent):\n  " + "\n  ".join(labels))
    return "\n".join(lines)


def _fmt_variant_annot(v: dict) -> str:
    """Compact annotation string shared by both variant block builders."""
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
    if biotype and biotype != "protein_coding":   # skip the default; keep non-coding annotation
        parts.append(f"bt={biotype}")
    return f" [{', '.join(parts)}]" if parts else ""


def _vcf_block_blinded(variants: list[dict]) -> str:
    """Variant lines WITHOUT gene names — only coordinates + pathogenicity annotations."""
    if not variants:
        return ""
    lines = ["\nVariants (gene identity withheld; rank by phenotype match):"]
    for v in variants[:30]:
        chrom = v.get("chrom", "?")
        pos   = v.get("pos", "?")
        ref   = v.get("ref", "?")
        alt   = v.get("alt", "?")
        lines.append(f"  chr{chrom}:{pos} {ref}>{alt}{_fmt_variant_annot(v)}")
    if len(variants) > 30:
        lines.append(f"  ... and {len(variants) - 30} more")
    return "\n".join(lines)


def _variant_gene_symbol(v: dict) -> str:
    """Return the HGNC gene symbol from a variant dict, handling swapped gene_id/gene_symbol fields."""
    gs  = (v.get("gene_symbol") or "").strip()
    gid = (v.get("gene_id")     or "").strip()
    # If gene_symbol is an ENSG ID the fields are swapped — use gene_id instead
    if gs and not gs.upper().startswith("ENSG"):
        return gs
    if gid and not gid.upper().startswith("ENSG"):
        return gid
    return gs or gid or "Unknown"


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
        lines.append(f"  {gene}: {coord}{annot}")
    if len(variants) > 50:
        lines.append(f"  ... and {len(variants) - 50} more")
    return "\n".join(lines)



def build_prompt(patient: dict, ablation: str, hpo_labels: dict,
                 task_prompt: str = _AIVA_TASK) -> str:
    """Build a benchmark prompt for the given ablation condition.

    Ablation conditions:
      hpo             — HPO + age only               (floor: phenotype-only)
      hpo+vcf         — HPO + blinded variants        (PRIMARY benchmark)
      hpo+vcf+disease — HPO + blinded variants + disease name
      hpo+vcf+genes   — HPO + variants WITH gene names, no disease name
                        (closed-list phenotype ranking; AIVA ranks ~10 named genes by HPO fit)

    task_prompt: the opening instruction line (see PROMPT_STYLES).
    """
    if ablation not in VALID_ABLATIONS:
        raise ValueError(f"Unknown ablation '{ablation}'. Choose from: {VALID_ABLATIONS}")

    parts = []
    parts.append(_age_block(patient.get("age", "")))
    parts.append(_hpo_block(patient["positive_hpo"], patient["negative_hpo"], hpo_labels))

    if ablation == "hpo":
        pass  # phenotype + age only

    elif ablation == "hpo+vcf":
        parts.append(_vcf_block_blinded(patient.get("variants", [])))

    elif ablation == "hpo+vcf+disease":
        if patient.get("disease_name"):
            parts.append(f"\nSuspected Disease: {patient['disease_name']}")
        parts.append(_vcf_block_blinded(patient.get("variants", [])))

    elif ablation == "hpo+vcf+genes":
        # Each variant line is paired with its HGNC gene symbol.
        # No disease name is given — AIVA must rank by phenotype-gene fit alone.
        parts.append(_vcf_block_named(patient.get("variants", [])))

    elif ablation == "hpo+vcf+genes+disease":
        # Full information: show candidate genes AND include suspected disease name
        if patient.get("disease_name"):
            parts.append(f"\nSuspected Disease: {patient['disease_name']}")
        parts.append(_vcf_block_named(patient.get("variants", [])))

    # Post-case instruction: task + force-inclusion of all candidate genes in output JSON.
    post_instruction = (
        f"\n{task_prompt}\n\n"
        "If gene symbols are listed in the `Candidate Genes` section, you MUST include every one in the returned `predicted_genes` array.\n"
        "If a gene is not supported as causal, list it last. Do not omit any gene from the list.\n\n"
        + _JSON_SCHEMA
    )

    parts.append(post_instruction)
    return "\n".join(p for p in parts if p)


# ──────────────────────────────────────────────────────────────
# AIVA API caller
# ──────────────────────────────────────────────────────────────
async def call_aiva(client: httpx.AsyncClient, prompt: str) -> str | None:
    if not API_KEY:
        print("ERROR: AIVA_API_TOKEN / AIVA_API_KEY not set.", file=sys.stderr)
        return None

    headers = {"Authorization": f"Bearer {API_KEY}"}
    data    = {"message": prompt, "model": MODEL, "stream": "false"}

    try:
        async with client.stream("POST", API_URL, data=data, headers=headers,
                                 timeout=REQUEST_TIMEOUT) as resp:
            if resp.status_code != 200:
                body = (await resp.aread()).decode("utf-8", errors="replace")
                print(f"  API error {resp.status_code}: {body[:300]}", file=sys.stderr)
                return None
            text = ""
            async for line in resp.aiter_lines():
                line = line.strip()
                if not line or not line.startswith("data: "):
                    continue
                payload = line[6:]
                if payload == "[DONE]":
                    break
                try:
                    ev = json.loads(payload)
                    if ev.get("type") in ("text_delta", "final_output"):
                        text += ev.get("content", "")
                except json.JSONDecodeError:
                    pass
            if not text:
                print("  Warning: API returned 200 but streamed no content.", file=sys.stderr)
            return text
    except Exception as exc:
        print(f"  Request failed: {exc}", file=sys.stderr)
        return None


# ──────────────────────────────────────────────────────────────
# Response parser
# ──────────────────────────────────────────────────────────────
def parse_predictions(response_text: str) -> tuple[list[str], dict]:
    """Parse ranked+scored gene list. Rank = array order; score = LLM confidence.

    Returns (deduped_gene_list, {GENE_UPPER: score}).
    """
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
                            gene = str(gene)
                            m = _ID_IN_PARENS.match(gene)
                            if m:
                                gene = m.group(1).strip()
                            genes.append(gene)
                            if scr is not None:
                                try:
                                    scores[gene.upper()] = float(scr)
                                except (ValueError, TypeError):
                                    pass
                    elif isinstance(entry, str):
                        genes.append(entry)
            break
        idx += max(span, 1)

    # De-duplicate preserving array order (= rank)
    seen: set[str] = set()
    deduped: list[str] = []
    for g in genes:
        if g.upper() not in seen:
            seen.add(g.upper())
            deduped.append(g)

    return deduped, scores


# ──────────────────────────────────────────────────────────────
# Ensembl → HGNC symbol resolver
# ──────────────────────────────────────────────────────────────
class EnsemblResolver:
    """Accumulates Ensembl→symbol mappings seen in loaded VCFs/JSONL.

    On construction, pre-populates from the HGNC complete set so that Ensembl
    IDs absent from the current run's variants can still be resolved.
    """

    def __init__(self):
        self._map: dict[str, str] = {}
        self._load_hgnc()

    def _load_hgnc(self):
        """Pre-populate from HGNC complete set (Ensembl gene ID → symbol)."""
        if not HGNC_PATH.exists():
            return
        with open(HGNC_PATH, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                sym = row.get("symbol", "").strip()
                ens = row.get("ensembl_gene_id", "").strip()
                if ens and sym:
                    self._map.setdefault(ens, sym)  # don't override VCF-observed mappings

    def update(self, annots: dict):
        for v in annots.values():
            eid = v.get("ensembl_id", "")
            sym = v.get("gene_symbol", "")
            if eid and sym:
                self._map[eid] = sym

    def update_from_variants(self, variants: list[dict]):
        for v in variants:
            eid = v.get("gene_id") or ""
            sym = v.get("gene_symbol") or ""
            # Some entries have the fields swapped (gene_symbol contains ENSG, gene_id contains symbol)
            if sym.startswith("ENSG") and eid and not eid.startswith("ENSG"):
                eid, sym = sym, eid
            if eid.startswith("ENSG") and sym and not sym.startswith("ENSG"):
                self._map[eid] = sym

    def resolve(self, ensembl_id: str) -> str:
        """Return HGNC symbol or the original ID if not found."""
        return self._map.get(ensembl_id, ensembl_id)

    def resolve_list(self, ids: list[str]) -> list[str]:
        return [self.resolve(i) for i in ids]


# ──────────────────────────────────────────────────────────────
# Patient input loader
# ──────────────────────────────────────────────────────────────
def load_disease_id_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    with open(path) as f:
        f.readline()  # skip header
        for line in f:
            parts = line.strip().split("\t", 1)
            if len(parts) == 2:
                mapping[parts[0].strip()] = parts[1].strip()
    return mapping


def load_patients(jsonl_path: Path, id_map: dict[str, str],
                  resolver: EnsemblResolver,
                  limit: int = 0) -> list[dict]:
    """Load patients from final_extracted.jsonl.

    Each patient dict has:
      patient_id, disease_id, disease_name, positive_hpo, negative_hpo,
      age, variants (enriched from VCF), true_genes (HGNC symbols),
      category, broad_category

    limit: if > 0, stop after reading this many patients (avoids loading
           all 42k VCF files when only a subset is needed).
    """
    patients: list[dict] = []
    with open(jsonl_path, encoding="utf-8") as f:
        for line in f:
            if limit and len(patients) >= limit:
                break
            line = line.strip()
            if not line:
                continue
            entry = json.loads(line)
            patient_id = entry["id"]

            # Enrich variants from VCF
            vcf_annots = load_vcf_annotations(patient_id)
            resolver.update(vcf_annots)

            variants = entry.get("variants", [])
            resolver.update_from_variants(variants)
            for v in variants:
                key = (str(v.get("chrom", "")), int(v.get("pos", 0)))
                if key in vcf_annots:
                    a = vcf_annots[key]
                    v["cadd"]       = a["cadd"]
                    v["polyphen2"]  = a["polyphen2"]
                    v["maf"]        = a["maf"]
                    v["maf_absent"] = a["maf_absent"]
                    v["csq"]        = a["csq"]
                    v["biotype"]    = a["biotype"]
                    if a["gene_symbol"]:
                        v["gene_symbol"] = a["gene_symbol"]

            did = str(entry["disease_id"])
            patients.append({
                "patient_id":    patient_id,
                "disease_id":    did,
                "disease_name":  id_map.get(did, f"Disease_{did}"),
                "positive_hpo":  sorted(entry.get("positive_phenotypes", {}).keys()),
                "negative_hpo":  sorted(entry.get("negative_phenotypes", {}).keys()),
                "age":           entry.get("age", "Unknown"),
                "variants":      variants,
                # true_genes as Ensembl IDs — resolved to HGNC symbols after all patients loaded
                "true_genes_raw": entry.get("true_genes", []),
                "category":       entry.get("category", ""),
                "broad_category": entry.get("broad_category", ""),
                "in_gene_gene_kg":    entry.get("in_gene_gene_kg", None),
                "in_disease_gene_kg": entry.get("in_disease_gene_kg", None),
            })

    # Resolve true_genes Ensembl IDs to HGNC symbols for all patients at once
    for p in patients:
        p["true_genes"] = resolver.resolve_list(p["true_genes_raw"])

    return patients


# ──────────────────────────────────────────────────────────────
# Evaluation metrics
# ──────────────────────────────────────────────────────────────
def rank_of_true_gene(predicted: list[str], true_genes: list[str]) -> int | None:
    """Return 1-based rank of the first true gene in the predicted list, or None.

    Normalises predicted gene strings before comparison so that AIVA output
    like "TEK (ENSG00000120156)" matches the true gene symbol "TEK".
    """
    true_upper = {g.upper() for g in true_genes}
    for i, g in enumerate(predicted, 1):
        norm = _ID_IN_PARENS.match(g.strip())
        g_upper = (norm.group(1).strip() if norm else g.strip()).upper()
        if g_upper in true_upper:
            return i
    return None


def compute_metrics(results: list[dict], at_ks: tuple[int, ...] = (1, 5, 10, 20, 50)) -> dict:
    """Compute aggregate benchmarking metrics over a list of result records.

    Each record must have:
      true_genes  : list[str]   (HGNC symbols)
      predicted   : list[str]   (HGNC symbols, ranked - AFTER force-inclusion)
      predicted_before_force : list[str] (BEFORE force-inclusion - TRUE model output)
      category    : str

    Recall@K metrics use predicted_before_force (true model performance).
    Ranking metrics (MRR, NDCG, mean/median rank) use predicted (post-force-inclusion).
    """
    ranks: list[int] = []
    ranks_before_force: list[int] = []
    found_at: dict[int, int] = {k: 0 for k in at_ks}
    found_at_after_force: dict[int, int] = {k: 0 for k in at_ks}
    reciprocal_ranks: list[float] = []
    dcg_list: list[float] = []
    no_prediction = 0
    retention_count = 0

    per_category: dict[str, list[int | None]] = defaultdict(list)
    per_category_before_force: dict[str, list[int | None]] = defaultdict(list)

    for r in results:
        predicted  = r.get("predicted", [])
        raw        = r.get("predicted_before_force", predicted)
        true_genes = r.get("true_genes", [])
        cat        = r.get("category", "unknown")

        if not predicted:
            no_prediction += 1
            per_category[cat].append(None)
            per_category_before_force[cat].append(None)
            reciprocal_ranks.append(0.0)
            dcg_list.append(0.0)
            for k in at_ks:
                pass  # not found
            continue

        # TRUE RECALL: rank in model output BEFORE force-inclusion
        rk_before_force = rank_of_true_gene(raw, true_genes)
        per_category_before_force[cat].append(rk_before_force)

        # Retention: whether the true gene appeared in the model output before we forced inclusion
        if rk_before_force is not None:
            retention_count += 1
            ranks_before_force.append(rk_before_force)
            for k in at_ks:
                if rk_before_force <= k:
                    found_at[k] += 1

        # Ranking metrics: use AFTER force-inclusion, but ONLY if gene has non-zero score
        # If score = 0.0, the gene was force-included and ranking is arbitrary (based on VCF order)
        rk = rank_of_true_gene(predicted, true_genes)
        scores = r.get("predicted_scores", {})

        # Check if true gene has a real score (not force-included with 0.0)
        true_gene_has_real_score = False
        if rk is not None and true_genes:
            for tg in true_genes:
                gene_score = scores.get(tg, 0.0)
                if gene_score > 0.0:
                    true_gene_has_real_score = True
                    break

        # Only count rank if gene has a real score from the model
        if rk is not None and true_gene_has_real_score:
            ranks.append(rk)
            reciprocal_ranks.append(1.0 / rk)
            dcg_list.append(1.0 / math.log2(rk + 1))
            per_category[cat].append(rk)
            for k in at_ks:
                if rk <= k:
                    found_at_after_force[k] += 1
        else:
            # Gene was force-included (score 0.0) or not found - treat as null rank
            reciprocal_ranks.append(0.0)
            dcg_list.append(0.0)
            per_category[cat].append(None)

    n = len(results)
    if n == 0:
        return {}

    idcg = 1.0  # ideal DCG when true gene is at rank 1

    metrics: dict = {
        "n_patients":       n,
        "no_prediction":    no_prediction,
        "mrr":              round(sum(reciprocal_ranks) / n, 4),
        "ndcg@10":          round(sum(dcg_list) / (n * idcg), 4),
        "mean_rank":        round(sum(ranks) / len(ranks), 2) if ranks else None,
        "median_rank":      float(sorted(ranks)[len(ranks) // 2]) if ranks else None,
        "mean_rank_before_force": round(sum(ranks_before_force) / len(ranks_before_force), 2) if ranks_before_force else None,
        "median_rank_before_force": float(sorted(ranks_before_force)[len(ranks_before_force) // 2]) if ranks_before_force else None,
        "retention_rate":   round(retention_count / n, 4),
        "filtered_out_rate": round(1.0 - retention_count / n, 4),
    }
    # TRUE Recall@K (before force-inclusion) - the real performance metric
    for k in at_ks:
        metrics[f"recall@{k}"] = round(found_at[k] / n, 4)

    # Recall@K after force-inclusion (for debugging/comparison)
    for k in at_ks:
        metrics[f"recall_after_force@{k}"] = round(found_at_after_force[k] / n, 4)

    # Per-category recall@K (BEFORE force-inclusion - true model performance)
    cat_metrics: dict[str, dict] = {}
    for cat, cat_ranks in per_category_before_force.items():
        n_cat = len(cat_ranks)
        cm: dict = {"n": n_cat}
        for k in at_ks:
            found_k = sum(1 for rk in cat_ranks if rk is not None and rk <= k)
            cm[f"recall@{k}"] = round(found_k / n_cat, 4)
        cat_metrics[cat] = cm
    metrics["per_category"] = cat_metrics

    # Per-category recall@K AFTER force-inclusion (for comparison)
    cat_metrics_after: dict[str, dict] = {}
    for cat, cat_ranks in per_category.items():
        n_cat = len(cat_ranks)
        cm: dict = {"n": n_cat}
        for k in at_ks:
            found_k = sum(1 for rk in cat_ranks if rk is not None and rk <= k)
            cm[f"recall_after_force@{k}"] = round(found_k / n_cat, 4)
        cat_metrics_after[cat] = cm
    metrics["per_category_after_force"] = cat_metrics_after

    return metrics


def print_metrics(metrics: dict, ablation: str):
    print(f"\n{'═' * 80}")
    print(f"  BENCHMARK RESULTS  —  ablation: {ablation}")
    print(f"{'═' * 80}")
    print(f"  Patients evaluated  : {metrics.get('n_patients')}")
    print(f"  No prediction       : {metrics.get('no_prediction')}")
    print(f"  Retention rate      : {metrics.get('retention_rate', 0):.1%}  (model found gene before force-inclusion)")
    print(f"  Filtered-out rate   : {metrics.get('filtered_out_rate', 0):.1%}  (model missed gene, force-added)")
    print()
    print(f"  TRUE MODEL PERFORMANCE (before force-inclusion):")
    print(f"  ─────────────────────────────────────────────────")
    for k in (1, 5, 10, 20, 50):
        key = f"recall@{k}"
        if key in metrics:
            print(f"    Recall@{k:<2}         : {metrics[key]:.1%}")
    mean_bf = metrics.get('mean_rank_before_force')
    median_bf = metrics.get('median_rank_before_force')
    if mean_bf or median_bf:
        print(f"    Mean rank          : {mean_bf}")
        print(f"    Median rank        : {median_bf}")

    print()
    print(f"  AFTER FORCE-INCLUSION (inflated by missing genes appended):")
    print(f"  ─────────────────────────────────────────────────────────────")
    for k in (1, 5, 10, 20, 50):
        key = f"recall_after_force@{k}"
        if key in metrics:
            print(f"    Recall@{k:<2}         : {metrics[key]:.1%}")
    print(f"    MRR                : {metrics.get('mrr')}")
    print(f"    NDCG@10            : {metrics.get('ndcg@10')}")
    print(f"    Mean rank          : {metrics.get('mean_rank')}")
    print(f"    Median rank        : {metrics.get('median_rank')}")

    print()
    cat_m = metrics.get("per_category", {})
    if cat_m:
        recall_ks = [k for k in (1, 5, 10, 20, 50) if f"recall@{k}" in next(iter(cat_m.values()), {})]
        header_str = " / ".join(f"@{k}" for k in recall_ks)
        print(f"  Per-category Recall {header_str} (BEFORE force-inclusion):")
        for cat, cm in sorted(cat_m.items()):
            recall_str = "  ".join(f"R@{k}={cm[f'recall@{k}']:.1%}" for k in recall_ks)
            print(f"    {cat[:45]:<45}  {recall_str}  (n={cm['n']})")
    print(f"{'═' * 80}\n")


# ──────────────────────────────────────────────────────────────
# Evaluate-only mode  (score an existing predictions file)
# ──────────────────────────────────────────────────────────────
def evaluate_only(predictions_path: Path, ablation: str):
    results: list[dict] = []
    with open(predictions_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                results.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    metrics = compute_metrics(results)
    print_metrics(metrics, ablation)

    summary_path = predictions_path.with_name(
        predictions_path.stem + "_evaluation_summary.json"
    )
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump({
                    "ablation":              ablation,
                    "metrics_explanation":   SUMMARY_METRICS_EXPLANATION,
                    "category_descriptions": SUMMARY_CATEGORY_DESCRIPTIONS,
                    "metrics":               metrics,
                    "timestamp":             datetime.now(timezone.utc).isoformat(),
                  }, f, indent=2, ensure_ascii=False)
    print(f"Evaluation summary saved to: {summary_path}")


# ──────────────────────────────────────────────────────────────
# Main processing loop
# ──────────────────────────────────────────────────────────────
async def process_patients(
    patients:      list[dict],
    output_path:   Path,
    ablation:      str,
    hpo_labels:    dict,
    dry_run:       bool = False,
    task_prompt:   str  = _AIVA_TASK,
    save_debug:    bool = False,
) -> list[dict]:
    """Submit each patient to AIVA, parse predictions, write results, return records."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Resume support — skip already-written patient IDs
    done_ids: set = set()
    results:  list[dict] = []
    if output_path.exists():
        with open(output_path, encoding="utf-8") as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    done_ids.add(str(rec["patient_id"]))
                    results.append(rec)
                except (json.JSONDecodeError, KeyError):
                    pass
        if done_ids:
            print(f"Resuming — {len(done_ids)} patients already done.")

    remaining = [p for p in patients if str(p["patient_id"]) not in done_ids]
    print(f"Patients to process: {len(remaining)} / {len(patients)}")

    sem     = asyncio.Semaphore(CONCURRENCY)
    counter = {"done": 0}
    total   = len(remaining)

    print(f"Concurrency: {CONCURRENCY} simultaneous requests")

    async def process_one(p: dict) -> dict | None:
        pid   = p["patient_id"]
        dname = p["disease_name"]
        prompt = build_prompt(p, ablation, hpo_labels, task_prompt=task_prompt)

        if dry_run:
            print(f"\n── Patient {pid}  ({dname}) ──")
            print(prompt)
            print("─" * 60)
            return None

        async with sem:                          # cap to CONCURRENCY parallel calls
            response_text = await call_aiva(client, prompt)

        predicted, scores = parse_predictions(response_text or "")
        predicted_before_force = list(predicted)
        true_genes = p.get("true_genes", [])

        # ── Force inclusion of all variant genes ──
        # Even if the model omits them, append them at the end so the true gene
        # (always in the variants list) is never "missing" from the ranked list.
        variant_genes = []
        seen_vg = set()
        for v in p.get("variants", []):
            gene = _variant_gene_symbol(v)
            if gene and gene != "Unknown" and gene.upper() not in seen_vg:
                seen_vg.add(gene.upper())
                variant_genes.append(gene)

        # Append missing genes (ranked last) and set score 0.0 for model-absent genes
        predicted_upper = {g.upper() for g in predicted}
        for g in variant_genes:
            if g.upper() not in predicted_upper:
                predicted.append(g)
                try:
                    scores[g.upper()] = 0.0
                except Exception:
                    scores[g.upper()] = 0.0

        rank = rank_of_true_gene(predicted, true_genes)

        rec = {
            "patient_id":        pid,
            "disease_id":        p["disease_id"],
            "disease_name":      dname,
            "ablation":          ablation,
            "true_genes":        true_genes,
            "predicted":         predicted,
            "predicted_before_force": predicted_before_force,
            "predicted_scores":  {g: scores.get(g.upper(), 0.0) for g in predicted},
            "rank_of_true_gene": rank,
            "category":          p.get("category", ""),
            "in_gene_gene_kg":   p.get("in_gene_gene_kg"),
            "in_disease_gene_kg":p.get("in_disease_gene_kg"),
            "timestamp":         datetime.now(timezone.utc).isoformat(),
        }
        if save_debug:
            rec["prompt"]       = prompt
            rec["raw_response"] = response_text or ""

        counter["done"] += 1

        # Build detailed status line showing top gene, true gene, and scores
        top_gene = predicted[0] if predicted else "NONE"
        top_score = scores.get(top_gene.upper(), 0.0) if predicted else 0.0

        true_gene_str = true_genes[0] if true_genes else "UNKNOWN"
        true_score = scores.get(true_gene_str.upper(), 0.0) if true_genes else 0.0

        # Determine status
        if rank == 1:
            status = "✓ CORRECT"
        elif rank is not None and rank <= 5:
            status = f"RANK_{rank}"
        elif rank is not None:
            status = f"RANK_{rank}"
        else:
            status = "✗ MISSED"

        # Check if gene was in raw prediction (before force)
        in_raw = "✓" if rank_of_true_gene(predicted_before_force, true_genes) is not None else "✗"

        print(f"  [{counter['done']}/{total}] Patient {pid}: "
              f"Top={top_gene}(s={top_score:.2f}) | "
              f"True={true_gene_str}(rank={rank if rank else 'None'}, s={true_score:.2f}, in_model={in_raw}) | "
              f"{status}")

        fout.write(json.dumps(rec, ensure_ascii=False) + "\n")
        fout.flush()
        return rec

    async with httpx.AsyncClient() as client:
        with open(output_path, "a", encoding="utf-8") as fout:
            if dry_run:
                for p in remaining:
                    await process_one(p)
            else:
                task_results = await asyncio.gather(
                    *[process_one(p) for p in remaining],
                    return_exceptions=True,
                )
                for r in task_results:
                    if isinstance(r, Exception):
                        print(f"  Task error: {r}", file=sys.stderr)
                    elif isinstance(r, dict):
                        results.append(r)

    return results


# ──────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="AIVA gene prioritisation benchmark.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--input", default=str(DEFAULT_INPUT),
                   help="Path to final_extracted.jsonl benchmark file.")
    p.add_argument("--output", default=None,
                   help="Output .jsonl file. Default: Results/<ablation>_predictions.jsonl")
    p.add_argument("--ablation", choices=VALID_ABLATIONS, default="hpo+vcf",
                   help=(
                       "Ablation condition: "
                       "'hpo' = HPO only (floor); "
                       "'hpo+vcf' = HPO + blinded variants (PRIMARY); "
                       "'hpo+vcf+disease' = HPO + blinded variants + disease name; "
                       "'hpo+vcf+genes' = HPO + named variants, no disease (closed-list ranking); "
                       "'hpo+vcf+genes+disease' = HPO + named variants + disease name (full-info)."
                   ))
    p.add_argument("--prompt-style", choices=list(PROMPT_STYLES), default="aiva",
                   help=(
                       "Opening task instruction sent to AIVA: "
                       "'aiva' = directs AIVA to query its internal databases and KG (default); "
                       "'generic' = minimal instruction, lets the model decide how to reason."
                   ))
    p.add_argument("--hpo-labels", default=None,
                   help="Optional path to HPO label TSV/JSON (HP:XXXXX→label). "
                        "Enriches prompts with human-readable phenotype names.")
    p.add_argument("--limit", type=int, default=0,
                   help="Cap the number of patients to process (0 = all).")
    p.add_argument("--batch-size", type=int, default=0,
                   help="Process only a batch of N random patients from remaining (0 = all remaining).")
    p.add_argument("--progress-file", type=str, default=str(DEFAULT_RESULTS / "progress.json"),
                   help="File to track processed patient IDs for batching.")
    p.add_argument("--dry-run", action="store_true",
                   help="Print prompts without calling the API.")
    p.add_argument("--evaluate-only",  metavar="PREDICTIONS_FILE", default=None,
                   help="Skip inference; compute metrics from an existing .jsonl results file.")
    p.add_argument("--no-eval", action="store_true",
                   help="Skip automatic evaluation step after inference.")
    p.add_argument("--save-debug", action="store_true",
                   help="Include 'prompt' and 'raw_response' fields in each output record (useful for debugging).")
    return p.parse_args()


def main():
    args = parse_args()

    # ── Evaluate-only shortcut ──
    if args.evaluate_only:
        path = Path(args.evaluate_only)
        if not path.exists():
            sys.exit(f"File not found: {path}")
        evaluate_only(path, ablation=args.ablation)
        return

    # ── Load patients ──
    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"Input file not found: {input_path}")

    id_map   = load_disease_id_map(DISEASE_LIST)
    resolver = EnsemblResolver()
    hpo_labels = load_hpo_labels(Path(args.hpo_labels) if args.hpo_labels else None)

    print(f"Loading patients from: {input_path}")
    
    # ── Handle batch processing vs simple limit ──
    if args.batch_size > 0:
        # For batch processing, need to load all patients first to track progress consistently
        all_patients = load_patients(input_path, id_map, resolver, limit=0)
        
        # Apply limit after loading if specified (for testing batch processing on subsets)
        if args.limit > 0:
            all_patients = all_patients[:args.limit]
            print(f"Limited to first {len(all_patients)} patients for batch processing")
        
        # Load progress
        progress_file = Path(args.progress_file)
        done_ids = set()
        if progress_file.exists():
            try:
                with open(progress_file, "r") as f:
                    done_ids = set(json.load(f))
            except Exception as e:
                print(f"Warning: Could not load progress file {progress_file}: {e}")
        
        print(f"Total patients: {len(all_patients)}, Already processed: {len(done_ids)}")
        
        # Filter remaining patients
        remaining_patients = [p for p in all_patients if str(p["patient_id"]) not in done_ids]
        
        if not remaining_patients:
            print("All patients have been processed.")
            return
        
        # Select batch
        if args.batch_size < len(remaining_patients):
            selected_patients = random.sample(remaining_patients, args.batch_size)
        else:
            selected_patients = remaining_patients
        
        print(f"Selected batch of {len(selected_patients)} patients to process.")
        patients = selected_patients
    else:
        # For simple runs, load only what we need (more efficient)
        load_limit = args.limit if args.limit > 0 else 0
        patients = load_patients(input_path, id_map, resolver, limit=load_limit)
        
    print(f"Processing {len(patients)} patients.")

    task_prompt = PROMPT_STYLES[args.prompt_style]

    # ── Output path ──
    # Include prompt-style in the default filename so aiva vs generic runs
    # never overwrite each other (e.g. hpo+vcf_aiva_predictions.jsonl).
    if args.output:
        output_path = Path(args.output)
    else:
        DEFAULT_RESULTS.mkdir(parents=True, exist_ok=True)
        output_path = DEFAULT_RESULTS / f"{args.ablation}_{args.prompt_style}_predictions.jsonl"

    print(f"Ablation condition : {args.ablation}")
    print(f"Prompt style       : {args.prompt_style}")
    print(f"Output             : {output_path}")
    print()

    # ── Run inference ──
    results = asyncio.run(
        process_patients(
            patients    = patients,
            output_path = output_path,
            ablation    = args.ablation,
            hpo_labels  = hpo_labels,
            dry_run     = args.dry_run,
            task_prompt = task_prompt,
            save_debug  = args.save_debug,
        )
    )

    if args.dry_run or args.no_eval:
        return

    # ── Update progress file for batch processing ──
    if args.batch_size > 0 and results:
        progress_file = Path(args.progress_file)
        done_ids = set()
        if progress_file.exists():
            try:
                with open(progress_file, "r") as f:
                    done_ids = set(json.load(f))
            except Exception:
                pass
        
        # Add newly processed IDs
        newly_done = {str(r.get("patient_id", r.get("id", ""))) for r in results if r.get("patient_id") or r.get("id")}
        all_done = done_ids | newly_done
        
        try:
            progress_file.parent.mkdir(parents=True, exist_ok=True)
            with open(progress_file, "w") as f:
                json.dump(list(all_done), f)
            print(f"Updated progress file: {len(all_done)} patients processed total.")
        except Exception as e:
            print(f"Warning: Could not save progress file {progress_file}: {e}")

    # ── Auto-evaluate ──
    if results:
        metrics = compute_metrics(results)
        print_metrics(metrics, args.ablation)
        summary_path = output_path.with_name(output_path.stem + "_evaluation_summary.json")
        with open(summary_path, "w", encoding="utf-8") as f:
            json.dump(
                {
                    "ablation":              args.ablation,
                    "input":                 str(input_path),
                    "metrics_explanation":   SUMMARY_METRICS_EXPLANATION,
                    "category_descriptions": SUMMARY_CATEGORY_DESCRIPTIONS,
                    "metrics":               metrics,
                    "timestamp":             datetime.now(timezone.utc).isoformat(),
                },
                f, indent=2, ensure_ascii=False,
            )
        print(f"Evaluation summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
