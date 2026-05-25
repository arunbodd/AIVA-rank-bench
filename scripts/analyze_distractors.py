"""
Hallucinated Gene Analysis — HGNC-aware edition
================================================
Classifies the top-1 prediction for every non-rank-1 patient into one of:

  format_correct        — AIVA returned "GENE (ENSGXXX)"; after stripping the
                          parenthetical it matches the true gene symbol exactly.
  alias_correct         — AIVA returned a known HGNC alias / previous symbol for
                          the true gene; not truly hallucinated.
  truly_hallucinated    — confirmed different gene (true gene resolvable, AIVA wrong).
  unresolvable_true     — true gene ID is a numeric Entrez ID not found in HGNC or
                          the mygene cache; cannot determine correctness.
  unlabelled            — gene IS in the VCF/candidate pool but was shown only as a
                          raw Ensembl ID in the prompt (no symbol label).

The older "hallucinated" bucket mapped to format_correct + alias_correct +
truly_hallucinated + unresolvable_true.  Separating them clarifies what fraction
of the "8.6% LLM-confabulated" claim in the manuscript represents genuine model
errors vs. benchmark artefacts.
"""

import csv
import json
import re
from pathlib import Path

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE          = Path("/data/AIVA_Benchmarking")
V2_PATH       = BASE / "AIVA-var-bench/Simulated_Datasets/Results/Final_Broad_Cat_noDiseaseName_v2.jsonl"
INPUT_PATH    = BASE / "AIVA-var-bench/Simulated_Datasets/data/final_extracted.jsonl"
HGNC_PATH     = BASE / "Lirical_benchmark/tools/lirical-cli-2.2.1/data/hgnc_complete_set.txt"
MYGENE_CACHE  = BASE / "AIVA-var-bench/data/ensg_symbol_cache.json"

# Normalise predicted gene strings: strip trailing "(ENSGXXX)" etc.
_ID_IN_PARENS = re.compile(r'^(.+?)\s*\([^)]+\)\s*$')


def normalize_pred(g: str) -> str:
    m = _ID_IN_PARENS.match(g.strip())
    return (m.group(1).strip() if m else g.strip()).upper()


# ── Load HGNC ──────────────────────────────────────────────────────────────────
def load_hgnc(path: Path):
    """Return (ensembl_map, entrez_map, alias_to_canonical).

    ensembl_map      : ENSGXXX  → canonical symbol (upper)
    entrez_map       : entrez_id (str) → canonical symbol (upper)
    alias_to_canonical: any alias / prev symbol (upper) → canonical symbol (upper)
    """
    ensembl_map: dict[str, str] = {}
    entrez_map:  dict[str, str] = {}
    alias_to_canonical: dict[str, str] = {}

    if not path.exists():
        print(f"WARNING: HGNC file not found at {path}")
        return ensembl_map, entrez_map, alias_to_canonical

    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            sym     = row.get("symbol",          "").strip().upper()
            ens     = row.get("ensembl_gene_id", "").strip()
            entrez  = row.get("entrez_id",       "").strip()
            aliases = row.get("alias_symbol",    "").strip().replace('"', '')
            prev    = row.get("prev_symbol",     "").strip().replace('"', '')

            if ens and sym:
                ensembl_map[ens] = sym
            if entrez and sym:
                entrez_map[entrez] = sym
            if sym:
                alias_to_canonical[sym] = sym  # canonical maps to itself
                for a in aliases.split("|"):
                    a = a.strip().upper()
                    if a:
                        alias_to_canonical.setdefault(a, sym)
                for p in prev.split("|"):
                    p = p.strip().upper()
                    if p:
                        alias_to_canonical.setdefault(p, sym)

    return ensembl_map, entrez_map, alias_to_canonical


# ── Load mygene cache ──────────────────────────────────────────────────────────
def load_mygene_cache(path: Path) -> dict[str, str]:
    if path.exists():
        with open(path, encoding="utf-8") as f:
            raw = json.load(f)
        return {k: v.upper() for k, v in raw.items() if v}
    return {}


# ── Build combined resolver ────────────────────────────────────────────────────
def build_resolver(variant_ensembl_map: dict[str, str],
                   hgnc_ens: dict[str, str],
                   hgnc_entrez: dict[str, str],
                   mygene: dict[str, str]) -> dict[str, str]:
    """Merge all sources: variant-observed > HGNC > mygene."""
    resolved: dict[str, str] = {}
    for src in (mygene, hgnc_ens, hgnc_entrez, variant_ensembl_map):
        for k, v in src.items():
            resolved.setdefault(k.upper(), v.upper())
    return resolved


def resolve_gene_id(gid: str, resolver: dict[str, str]) -> str | None:
    """Resolve an Ensembl or Entrez gene ID to a canonical symbol."""
    gid = gid.strip()
    if not gid:
        return None
    u = gid.upper()
    # Already a symbol if it doesn't look like a numeric/ENSG ID
    if not u.startswith("ENSG") and not gid.isdigit():
        return u
    return resolver.get(u) or resolver.get(gid)


# ── Main ───────────────────────────────────────────────────────────────────────
def main():
    print("Loading HGNC …")
    hgnc_ens, hgnc_entrez, alias_to_canonical = load_hgnc(HGNC_PATH)
    print(f"  {len(hgnc_ens):,} Ensembl→symbol entries")
    print(f"  {len(alias_to_canonical):,} alias entries")

    mygene_cache = load_mygene_cache(MYGENE_CACHE)
    print(f"  {len(mygene_cache):,} mygene cache entries")

    # 1. Load true genes from input
    true_genes_map: dict[int, list[str]] = {}
    with open(INPUT_PATH, encoding="utf-8") as f:
        for line in f:
            if line.strip():
                e = json.loads(line)
                true_genes_map[int(e["id"])] = e.get("true_genes", [])

    # 2. Build variant-observed Ensembl→symbol map from results file
    variant_ensembl_map: dict[str, str] = {}
    records: list[dict] = []
    with open(V2_PATH, encoding="utf-8") as f:
        for line in f:
            if not line.strip():
                continue
            rec = json.loads(line)
            records.append(rec)
            for v in rec.get("variants", []):
                eid = v.get("gene_id") or ""
                sym = v.get("gene_symbol") or ""
                if sym.startswith("ENSG") and eid and not eid.startswith("ENSG"):
                    eid, sym = sym, eid
                if eid.startswith("ENSG") and sym and not sym.startswith("ENSG"):
                    variant_ensembl_map[eid.upper()] = sym.upper()

    resolver = build_resolver(variant_ensembl_map, hgnc_ens, hgnc_entrez, mygene_cache)
    print(f"\nTotal result records: {len(records)}")

    # 3. Identify non-rank-1 cases using normalised matching
    rank1_count = 0
    non_rank1: list[dict] = []

    for r in records:
        pid = int(r["patient_id"])
        true_raw = true_genes_map.get(pid, [])
        true_syms: set[str] = set()
        for gid in true_raw:
            sym = resolve_gene_id(gid, resolver)
            if sym:
                true_syms.add(sym)
            else:
                true_syms.add(gid.upper())

        preds = r.get("predicted_genes", [])
        rank = None
        for i, pg in enumerate(preds, 1):
            norm = normalize_pred(pg["gene"])
            if norm in true_syms:
                rank = i
                break
            # Also check HGNC alias: pred may be a previous/alias symbol for true gene
            canon = alias_to_canonical.get(norm)
            if canon and canon in true_syms:
                rank = i
                break

        if rank == 1:
            rank1_count += 1
        else:
            r["_true_syms"] = true_syms
            r["_true_raw"]  = true_raw
            non_rank1.append(r)

    print(f"Rank-1 (correct) after normalisation: {rank1_count} ({rank1_count/len(records):.1%})")
    print(f"Non-rank-1 cases: {len(non_rank1)}")

    # 4. Classify non-rank-1 top predictions
    counters = {
        "format_correct":     0,
        "alias_correct":      0,
        "truly_hallucinated": 0,
        "unresolvable_true":  0,
        "unlabelled":         0,
        "other":              0,
        "no_prediction":      0,
    }
    examples: dict[str, list] = {k: [] for k in counters}

    for r in non_rank1:
        preds = r.get("predicted_genes", [])
        if not preds:
            counters["no_prediction"] += 1
            continue

        top_gene     = preds[0]["gene"]
        top_norm     = normalize_pred(top_gene)   # stripped, upper
        top_raw_upper = top_gene.strip().upper()  # un-stripped, upper
        true_syms    = r["_true_syms"]
        true_raw     = r["_true_raw"]
        pid          = int(r["patient_id"])

        # Build candidate set for this patient
        cand_syms: set[str] = set()
        cand_ensg: set[str] = set()
        for v in r.get("variants", []):
            sym = v.get("gene_symbol", "")
            eid = v.get("gene_id", "")
            if sym: cand_syms.add(sym.upper())
            if eid: cand_ensg.add(eid.upper())
        for cg in r.get("all_candidate_genes", []):
            cand_ensg.add(cg.upper())
            resolved = resolver.get(cg.upper())
            if resolved:
                cand_syms.add(resolved)

        in_candidates = (top_norm in cand_syms or top_norm in cand_ensg
                         or top_raw_upper in cand_syms or top_raw_upper in cand_ensg)

        if in_candidates:
            # Check if it was labelled in the prompt
            prompt = r.get("prompt", "")
            is_labeled = False
            in_cand_section = False
            for line in prompt.split("\n"):
                if "Candidate Genes" in line:
                    in_cand_section = True
                    continue
                if in_cand_section:
                    parts = line.strip().split()
                    if parts and parts[0].upper() == top_norm:
                        is_labeled = True
                        break
            if not is_labeled:
                counters["unlabelled"] += 1
                examples["unlabelled"].append((pid, top_gene))
            else:
                counters["other"] += 1
        else:
            # Not in candidate pool — determine why
            # (a) Format-only issue: top_norm matches true gene, but raw didn't
            if top_norm in true_syms and top_raw_upper not in true_syms:
                counters["format_correct"] += 1
                examples["format_correct"].append((pid, top_gene, list(true_syms)))
                continue

            # (b) Alias / previous symbol match
            canon_of_pred = alias_to_canonical.get(top_norm)
            if canon_of_pred and canon_of_pred in true_syms:
                counters["alias_correct"] += 1
                examples["alias_correct"].append((pid, top_gene, canon_of_pred))
                continue

            # (c) Check if any true gene ID is unresolvable
            all_unresolvable = all(
                gid.isdigit() and gid not in mygene_cache and gid not in hgnc_entrez
                for gid in true_raw
                if gid.isdigit()
            ) and any(gid.isdigit() for gid in true_raw)

            if all_unresolvable:
                counters["unresolvable_true"] += 1
                examples["unresolvable_true"].append((pid, top_gene, true_raw))
            else:
                counters["truly_hallucinated"] += 1
                examples["truly_hallucinated"].append((pid, top_gene, list(true_syms)))

    # 5. Report
    total = len(records)
    print("\n" + "=" * 70)
    print("HALLUCINATION ANALYSIS RESULTS")
    print("=" * 70)
    print(f"{'Category':<30} {'Count':>7}  {'% of total':>10}  {'% of non-rank-1':>16}")
    print("-" * 70)
    nr1 = len(non_rank1)
    for cat, cnt in counters.items():
        pct_total = cnt / total * 100
        pct_nr1   = cnt / nr1   * 100 if nr1 else 0
        print(f"  {cat:<28} {cnt:>7}  {pct_total:>9.1f}%  {pct_nr1:>15.1f}%")
    print("-" * 70)
    old_hallucinated = (counters["format_correct"] + counters["alias_correct"]
                        + counters["truly_hallucinated"] + counters["unresolvable_true"])
    print(f"  {'OLD hallucinated total':<28} {old_hallucinated:>7}  "
          f"{old_hallucinated/total*100:>9.1f}%  {old_hallucinated/nr1*100:>15.1f}%")
    print("=" * 70)

    print("\nNotes:")
    print(f"  format_correct   : AIVA was RIGHT — prediction stripped of parenthetical ID matches true gene")
    print(f"  alias_correct    : AIVA was RIGHT — prediction is a known HGNC alias/prev symbol")
    print(f"  truly_hallucinated: confirmed wrong gene (true gene resolvable)")
    print(f"  unresolvable_true: true gene ID is a numeric Entrez ID with no HGNC/mygene mapping")
    print(f"  unlabelled       : gene in VCF pool but shown as Ensembl ID only in prompt")

    if examples["format_correct"]:
        print("\nformat_correct examples:")
        for pid, pred, true in examples["format_correct"]:
            print(f"  pid={pid}: AIVA='{pred}' → true={true}")

    if examples["alias_correct"]:
        print("\nalias_correct examples:")
        for pid, pred, canon in examples["alias_correct"]:
            print(f"  pid={pid}: AIVA='{pred}' (alias of {canon})")

    if examples["truly_hallucinated"][:5]:
        print("\ntruly_hallucinated examples (first 5):")
        for pid, pred, true in examples["truly_hallucinated"][:5]:
            print(f"  pid={pid}: AIVA='{pred}', true={true[:3]}")

    if examples["unresolvable_true"][:5]:
        print("\nunresolvable_true examples (first 5):")
        for pid, pred, raw_ids in examples["unresolvable_true"][:5]:
            print(f"  pid={pid}: AIVA='{pred}', true_raw={raw_ids}")

    if examples["unlabelled"][:5]:
        print("\nunlabelled examples (first 5):")
        for pid, gene in examples["unlabelled"][:5]:
            print(f"  pid={pid}: '{gene}'")


if __name__ == "__main__":
    main()
