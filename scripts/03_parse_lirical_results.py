#!/usr/bin/env python3
"""
03_parse_lirical_results.py
Parse LIRICAL JSON output files and extract per-case gene rankings.

LIRICAL JSON structure (v2.2.1):
  {
    "analysisData": {...},
    "analysisMetadata": {...},
    "analysisResults": [
      {
        "diseaseId": "OMIM:123456",
        "posttestProbability": 0.85,
        "genotypeLR": {
          "geneId": {"id": "NCBIGene:xxxx", "symbol": "GENE_SYMBOL"},
          "lr": 10.0
        },
        "compositeLR": 5.0,
        ...
      },
      ...
    ]
  }

Output: results/lirical_rankings.jsonl
  Each record: patient_id, disease_id, true_gene, broad_category,
               lirical_rank (1-based, None if not found),
               lirical_top_genes (list of {gene, prob})

Usage:
    python scripts/03_parse_lirical_results.py [--verbose] [--delete-after]
    python scripts/03_parse_lirical_results.py --delete-after   # parse + free disk space
"""

import argparse
import json
import sys
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_PATH = REPO_ROOT / "inputs" / "manifest.jsonl"
OUTPUTS_DIR = REPO_ROOT / "outputs"
RANKINGS_PATH = REPO_ROOT / "results" / "lirical_rankings.jsonl"


# ── LIRICAL JSON parsing ──────────────────────────────────────────────────────

def parse_lirical_json(json_path: Path) -> list:
    """
    Parse a LIRICAL output JSON file.
    Returns a list of {gene_symbol, prob} sorted by prob descending (deduplicated by gene).
    Gene is taken from each analysisResult's genotypeLR.geneId.symbol.
    """
    with open(json_path) as f:
        data = json.load(f)

    raw_results = data.get("analysisResults", [])
    if not raw_results:
        return []

    # Collect best probability per gene symbol
    gene_best_prob: dict = {}
    for res in raw_results:
        prob = res.get("posttestProbability", 0.0)
        genotype_lr = res.get("genotypeLR")
        if genotype_lr is None:
            # Phenotype-only result — no gene associated; skip for gene ranking
            continue
        gene_id_obj = genotype_lr.get("geneId", {})
        symbol = gene_id_obj.get("symbol", "")
        if not symbol:
            # Try alt key names used in some LIRICAL versions
            symbol = gene_id_obj.get("name", "") or gene_id_obj.get("value", "")
        if not symbol:
            continue
        if symbol not in gene_best_prob or prob > gene_best_prob[symbol]:
            gene_best_prob[symbol] = prob

    # Sort by probability descending
    ranked = sorted(gene_best_prob.items(), key=lambda x: x[1], reverse=True)
    return [{"gene": sym, "prob": round(prob, 6)} for sym, prob in ranked]


def find_rank(ranked_genes: list, true_symbol: str) -> int | None:
    """Return 1-based rank of true_symbol in ranked_genes list, or None if not found."""
    if not true_symbol:
        return None
    # Case-insensitive match as fallback
    true_upper = true_symbol.upper()
    for i, entry in enumerate(ranked_genes):
        if entry["gene"] == true_symbol or entry["gene"].upper() == true_upper:
            return i + 1
    return None


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Parse LIRICAL results into rankings")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Print per-case details")
    parser.add_argument("--top-n", type=int, default=20,
                        help="Max top genes to save per case (default: 20)")
    parser.add_argument("--delete-after", action="store_true",
                        help="Delete each large JSON file after parsing (reclaims disk space)")
    args = parser.parse_args()

    if not MANIFEST_PATH.exists():
        print(f"ERROR: Manifest not found: {MANIFEST_PATH}")
        print("  → Run: python scripts/01_prepare_inputs.py")
        sys.exit(1)

    # Load manifest
    manifest = {}
    with open(MANIFEST_PATH) as f:
        for line in f:
            entry = json.loads(line)
            manifest[entry["patient_id"]] = entry

    total = len(manifest)
    stats = {"parsed": 0, "missing": 0, "parse_error": 0, "gene_found": 0, "gene_not_found": 0}
    rankings = []

    print(f"Parsing LIRICAL results for {total} cases...")

    for pid, entry in manifest.items():
        json_path = OUTPUTS_DIR / str(pid) / f"{pid}.json"

        if not json_path.exists():
            stats["missing"] += 1
            if args.verbose:
                print(f"  [MISSING] patient {pid}")
            continue

        try:
            ranked_genes = parse_lirical_json(json_path)
        except Exception as e:
            stats["parse_error"] += 1
            print(f"  [PARSE ERROR] patient {pid}: {e}", file=sys.stderr)
            continue
        finally:
            if args.delete_after and json_path.exists():
                try:
                    json_path.unlink()
                except OSError:
                    pass

        stats["parsed"] += 1
        true_sym = entry.get("true_gene_symbol")
        rank = find_rank(ranked_genes, true_sym)

        if rank is not None:
            stats["gene_found"] += 1
        else:
            stats["gene_not_found"] += 1
            if args.verbose:
                print(f"  [NOT FOUND] patient {pid}, true={true_sym}, "
                      f"top3={[g['gene'] for g in ranked_genes[:3]]}")

        record = {
            "patient_id": pid,
            "disease_id": entry.get("disease_id"),
            "broad_category": entry.get("broad_category"),
            "category": entry.get("category"),
            "true_gene_symbol": true_sym,
            "true_gene_ensg": entry.get("true_gene_ensg"),
            "lirical_rank": rank,
            "lirical_n_genes": len(ranked_genes),
            "lirical_top_genes": ranked_genes[: args.top_n],
        }
        rankings.append(record)

    # Write output
    RANKINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(RANKINGS_PATH, "w") as out:
        for rec in rankings:
            out.write(json.dumps(rec) + "\n")

    print(f"\n=== Parsing complete ===")
    print(f"  Parsed OK    : {stats['parsed']}")
    print(f"  Missing      : {stats['missing']}")
    print(f"  Parse errors : {stats['parse_error']}")
    print(f"  Gene found   : {stats['gene_found']} ({100*stats['gene_found']/max(stats['parsed'],1):.1f}%)")
    print(f"  Gene missing : {stats['gene_not_found']}")
    print(f"\nRankings written: {RANKINGS_PATH} ({len(rankings)} records)")
    print("Next step: python scripts/04_compare_metrics.py")


if __name__ == "__main__":
    main()
