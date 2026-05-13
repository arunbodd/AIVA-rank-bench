#!/usr/bin/env python3
"""
07_validate_runs.py
Audit all tool outputs before publication table generation.

Checks performed:
  1. Duplicate patient_ids in lirical_rankings.jsonl (FAIL if any found)
  2. All manifest PIDs present in lirical_rankings.jsonl
  3. Distribution of lirical_n_genes, lirical_rank, posttestProbability
  4. Cross-reference lirical_errors.jsonl with rankings
     (are any error PIDs also in rankings? → confirms retry success)
  5. AIVA predicted_genes non-empty check
  6. Optional: Exomiser rankings audit (if --include-exomiser)

Output:
  - results/validation_report.json  (machine-readable summary)
  - Console summary

Usage:
    python scripts/07_validate_runs.py
    python scripts/07_validate_runs.py --include-exomiser
"""

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
AIVA_ROOT = Path("/data/AIVA_Benchmarking/AIVA-var-bench/Simulated_Datasets")
AIVA_RESULTS = AIVA_ROOT / "Results" / "Final_Broad_Cat_noDiseaseName_v2.jsonl"
MANIFEST_PATH = REPO_ROOT / "inputs" / "manifest.jsonl"
LIRICAL_RANKINGS = REPO_ROOT / "results" / "lirical_rankings.jsonl"
LIRICAL_ERRORS = REPO_ROOT / "results" / "lirical_errors.jsonl"
EXOMISER_RANKINGS  = REPO_ROOT / "results" / "exomiser_rankings.jsonl"
EXOMISER_ERRORS    = REPO_ROOT / "results" / "exomiser_errors.jsonl"
CLINPRIOR_RANKINGS = REPO_ROOT / "results" / "clinprior_rankings.jsonl"
CLINPRIOR_ERRORS   = REPO_ROOT / "results" / "clinprior_errors.jsonl"
VALIDATION_REPORT = REPO_ROOT / "results" / "validation_report.json"


# ── helpers ───────────────────────────────────────────────────────────────────

def load_jsonl(path: Path) -> list:
    """Load all records from a JSONL file."""
    records = []
    with open(path) as f:
        for i, line in enumerate(f, 1):
            line = line.strip()
            if not line:
                continue
            try:
                records.append(json.loads(line))
            except json.JSONDecodeError as e:
                print(f"  [WARN] {path.name} line {i}: JSON decode error — {e}")
    return records


def percentile(sorted_vals: list, p: float) -> float:
    """Return p-th percentile (0–100) of a sorted list."""
    if not sorted_vals:
        return float("nan")
    idx = (p / 100.0) * (len(sorted_vals) - 1)
    lo = int(idx)
    hi = min(lo + 1, len(sorted_vals) - 1)
    return sorted_vals[lo] + (sorted_vals[hi] - sorted_vals[lo]) * (idx - lo)


def summarize_numeric(vals: list, name: str) -> dict:
    """Return summary statistics for a list of numeric values."""
    if not vals:
        return {"n": 0, "mean": None, "median": None, "min": None, "max": None,
                "p25": None, "p75": None}
    sv = sorted(vals)
    n = len(sv)
    return {
        "n": n,
        "mean": round(sum(sv) / n, 4),
        "median": round(percentile(sv, 50), 4),
        "min": sv[0],
        "max": sv[-1],
        "p25": round(percentile(sv, 25), 4),
        "p75": round(percentile(sv, 75), 4),
    }


# ── LIRICAL audit ─────────────────────────────────────────────────────────────

def audit_lirical(report: dict) -> bool:
    """
    Full audit of lirical_rankings.jsonl.

    Returns True if audit passes (no critical issues), False otherwise.
    """
    section = {}
    ok = True

    print(f"\n{'='*65}")
    print("  LIRICAL Rankings Audit")
    print(f"{'='*65}")

    if not LIRICAL_RANKINGS.exists():
        print(f"  ERROR: {LIRICAL_RANKINGS} not found")
        report["lirical"] = {"error": "rankings file not found"}
        return False

    rankings = load_jsonl(LIRICAL_RANKINGS)
    section["total_lines"] = len(rankings)
    print(f"  Total lines parsed        : {len(rankings)}")

    # 1. Duplicate patient_ids
    pids = [r.get("patient_id") for r in rankings]
    pid_counts = Counter(pids)
    duplicates = {pid: cnt for pid, cnt in pid_counts.items() if cnt > 1}
    section["duplicate_pids"] = duplicates
    if duplicates:
        print(f"  [FAIL] Duplicate patient_ids found: {list(duplicates.keys())[:10]}")
        ok = False
    else:
        print(f"  [PASS] No duplicate patient_ids")

    # 2. Cross-check manifest PIDs
    if not MANIFEST_PATH.exists():
        print(f"  [WARN] Manifest not found: {MANIFEST_PATH}")
        section["manifest_check"] = "skipped"
    else:
        manifest = load_jsonl(MANIFEST_PATH)
        manifest_pids = {r["patient_id"] for r in manifest}
        ranking_pids = set(pids)
        missing_from_rankings = manifest_pids - ranking_pids
        extra_in_rankings = ranking_pids - manifest_pids
        section["manifest_pids"] = len(manifest_pids)
        section["ranking_pids"] = len(ranking_pids)
        section["missing_from_rankings"] = sorted(missing_from_rankings)
        section["extra_in_rankings"] = sorted(extra_in_rankings)

        if missing_from_rankings:
            print(f"  [FAIL] {len(missing_from_rankings)} manifest PIDs not in rankings:")
            for pid in sorted(missing_from_rankings)[:10]:
                print(f"         {pid}")
            if len(missing_from_rankings) > 10:
                print(f"         ... and {len(missing_from_rankings) - 10} more")
            ok = False
        else:
            print(f"  [PASS] All {len(manifest_pids)} manifest PIDs present in rankings")

        if extra_in_rankings:
            print(f"  [WARN] {len(extra_in_rankings)} extra PIDs in rankings (not in manifest)")

    # 3. Rank distribution
    ranks = [r.get("lirical_rank") for r in rankings]
    rank_none = sum(1 for r in ranks if r is None)
    rank_found = [r for r in ranks if r is not None]
    section["rank_absent"] = rank_none
    section["rank_absent_pct"] = round(100 * rank_none / len(ranks), 2) if ranks else None
    section["rank_distribution"] = summarize_numeric(rank_found, "lirical_rank")

    print(f"\n  Rank distribution:")
    print(f"    Absent (rank=None)   : {rank_none}  ({section['rank_absent_pct']:.1f}%)")
    if rank_found:
        rs = summarize_numeric(rank_found, "lirical_rank")
        print(f"    Found cases          : {rs['n']}")
        print(f"    Min / Median / Max   : {rs['min']} / {rs['median']} / {rs['max']}")
        print(f"    Mean                 : {rs['mean']}")

    # Rank bins
    bins = {1: 0, 3: 0, 5: 0, 10: 0, 20: 0, 50: 0, ">50": 0}
    for r in rank_found:
        if r <= 1:
            bins[1] += 1
        elif r <= 3:
            bins[3] += 1
        elif r <= 5:
            bins[5] += 1
        elif r <= 10:
            bins[10] += 1
        elif r <= 20:
            bins[20] += 1
        elif r <= 50:
            bins[50] += 1
        else:
            bins[">50"] += 1
    section["rank_bins"] = bins
    print(f"    Rank bins: {bins}")

    # 4. n_genes distribution
    n_genes_vals = [r.get("lirical_n_genes", 0) for r in rankings]
    section["n_genes_distribution"] = summarize_numeric(n_genes_vals, "lirical_n_genes")
    ng = summarize_numeric(n_genes_vals, "lirical_n_genes")
    print(f"\n  n_genes distribution (LIRICAL output genes per case):")
    print(f"    Min / Median / Max   : {ng['min']} / {ng['median']} / {ng['max']}")
    print(f"    Mean                 : {ng['mean']}")
    zero_gene_cases = sum(1 for v in n_genes_vals if v == 0)
    section["zero_gene_cases"] = zero_gene_cases
    if zero_gene_cases:
        print(f"    [WARN] {zero_gene_cases} cases with 0 output genes")

    # 5. posttestProbability distribution (from top_genes[0])
    probs = []
    zero_prob_cases = 0
    for r in rankings:
        top_genes = r.get("lirical_top_genes", [])
        if top_genes:
            p = top_genes[0].get("prob", 0.0)
            probs.append(p)
            if p == 0.0:
                zero_prob_cases += 1
    section["top_prob_distribution"] = summarize_numeric(probs, "top_prob")
    section["zero_top_prob_cases"] = zero_prob_cases
    nonzero = sum(1 for p in probs if p > 0)
    print(f"\n  Top gene posttestProbability (stored in lirical_top_genes[0].prob):")
    print(f"    Cases with top gene   : {len(probs)}")
    print(f"    Non-zero prob cases   : {nonzero}  ({100*nonzero/len(probs):.1f}%)")
    print(f"    Zero prob cases       : {zero_prob_cases}  ({100*zero_prob_cases/len(probs):.1f}%)")
    if probs:
        nonzero_probs = sorted(p for p in probs if p > 0)
        if nonzero_probs:
            pp = summarize_numeric(nonzero_probs, "nonzero_probs")
            print(f"    Non-zero min/max     : {pp['min']:.2e} / {pp['max']:.4f}")
            print(f"    (Note: values <5e-7 are normal for --global mode; see plan DOC note)")

    # 6. Cross-reference error log
    if LIRICAL_ERRORS.exists():
        errors = load_jsonl(LIRICAL_ERRORS)
        error_pids = {e.get("patient_id") for e in errors}
        ranking_pid_set = set(pids)
        errors_also_ranked = error_pids & ranking_pid_set
        errors_not_ranked = error_pids - ranking_pid_set
        section["error_log"] = {
            "total_errors": len(errors),
            "error_pids_in_rankings": len(errors_also_ranked),
            "error_pids_missing_from_rankings": len(errors_not_ranked),
            "error_causes": Counter(
                e.get("error", "")[:60] for e in errors
            ).most_common(5),
        }
        print(f"\n  Error log ({LIRICAL_ERRORS.name}): {len(errors)} records")
        print(f"    Error PIDs also in rankings  : {len(errors_also_ranked)}  (successfully retried)")
        print(f"    Error PIDs missing from rankings: {len(errors_not_ranked)}")
        if errors_not_ranked:
            print(f"    [WARN] These PIDs have no rankings entry: {sorted(errors_not_ranked)[:10]}")
        if errors:
            print("    Top error causes:")
            for cause, cnt in section["error_log"]["error_causes"]:
                print(f"      ({cnt}×) {cause}")
    else:
        section["error_log"] = {"note": "lirical_errors.jsonl not found"}
        print(f"\n  [INFO] No error log found at {LIRICAL_ERRORS}")

    report["lirical"] = section
    return ok


# ── AIVA audit ────────────────────────────────────────────────────────────────

def audit_aiva(report: dict) -> bool:
    """Audit AIVA results for non-empty predicted_genes and sort order."""
    section = {}
    ok = True

    print(f"\n{'='*65}")
    print("  AIVA Results Audit")
    print(f"{'='*65}")

    if not AIVA_RESULTS.exists():
        print(f"  ERROR: {AIVA_RESULTS} not found")
        report["aiva"] = {"error": "results file not found"}
        return False

    aiva_recs = load_jsonl(AIVA_RESULTS)
    section["total_cases"] = len(aiva_recs)
    print(f"  Total cases               : {len(aiva_recs)}")

    # Check for empty predicted_genes
    empty_pg = [r.get("patient_id") for r in aiva_recs
                if not r.get("predicted_genes")]
    section["empty_predicted_genes"] = len(empty_pg)
    if empty_pg:
        print(f"  [WARN] {len(empty_pg)} cases with empty predicted_genes")
        for pid in empty_pg[:5]:
            print(f"         {pid}")
    else:
        print(f"  [PASS] All cases have non-empty predicted_genes")

    # Check sort order of predicted_genes
    unsorted_count = 0
    n_genes_vals = []
    for rec in aiva_recs:
        genes = rec.get("predicted_genes", [])
        n_genes_vals.append(len(genes))
        if len(genes) >= 2:
            scores = [g.get("score", 0) for g in genes]
            if scores != sorted(scores, reverse=True):
                unsorted_count += 1
    section["unsorted_predicted_genes"] = unsorted_count
    ng = summarize_numeric(n_genes_vals, "n_genes")
    section["n_genes_distribution"] = ng

    if unsorted_count:
        print(f"  [WARN] {unsorted_count} cases with non-sorted predicted_genes")
        print(f"         (04_compare_metrics.py applies defensive sort before rank lookup)")
    else:
        print(f"  [PASS] All predicted_genes lists are sorted by score DESC")

    print(f"\n  Predicted genes per case:")
    print(f"    Min / Median / Max   : {ng['min']} / {ng['median']} / {ng['max']}")
    print(f"    Mean                 : {ng['mean']}")

    report["aiva"] = section
    return ok


# ── Exomiser audit ────────────────────────────────────────────────────────────

def audit_exomiser(report: dict) -> bool:
    """Audit Exomiser rankings."""
    section = {}
    ok = True

    print(f"\n{'='*65}")
    print("  Exomiser Rankings Audit")
    print(f"{'='*65}")

    if not EXOMISER_RANKINGS.exists():
        print(f"  [SKIP] {EXOMISER_RANKINGS} not found (run 06_run_exomiser.py first)")
        report["exomiser"] = {"skipped": "rankings file not found"}
        return True  # Not an error — Exomiser may not have run yet

    rankings = load_jsonl(EXOMISER_RANKINGS)
    section["total_lines"] = len(rankings)
    print(f"  Total lines parsed        : {len(rankings)}")

    pids = [r.get("patient_id") for r in rankings]
    pid_counts = Counter(pids)
    duplicates = {pid: cnt for pid, cnt in pid_counts.items() if cnt > 1}
    section["duplicate_pids"] = duplicates
    if duplicates:
        print(f"  [FAIL] Duplicate patient_ids found: {list(duplicates.keys())[:10]}")
        ok = False
    else:
        print(f"  [PASS] No duplicate patient_ids")

    if MANIFEST_PATH.exists():
        manifest = load_jsonl(MANIFEST_PATH)
        manifest_pids = {r["patient_id"] for r in manifest}
        missing = manifest_pids - set(pids)
        section["missing_from_rankings"] = len(missing)
        if missing:
            print(f"  [WARN] {len(missing)} manifest PIDs not in Exomiser rankings")
        else:
            print(f"  [PASS] All manifest PIDs present in Exomiser rankings")

    ranks = [r.get("exomiser_rank") for r in rankings]
    rank_none = sum(1 for r in ranks if r is None)
    rank_found = [r for r in ranks if r is not None]
    section["rank_absent"] = rank_none
    section["rank_distribution"] = summarize_numeric(rank_found, "exomiser_rank")

    print(f"\n  Rank distribution:")
    print(f"    Absent (rank=None)   : {rank_none}  ({100*rank_none/len(ranks):.1f}%)")
    if rank_found:
        rs = summarize_numeric(rank_found, "exomiser_rank")
        print(f"    Min / Median / Max   : {rs['min']} / {rs['median']} / {rs['max']}")
        print(f"    Mean                 : {rs['mean']}")

    n_genes_vals = [r.get("exomiser_n_genes", 0) for r in rankings]
    section["n_genes_distribution"] = summarize_numeric(n_genes_vals, "exomiser_n_genes")
    ng = summarize_numeric(n_genes_vals, "exomiser_n_genes")
    print(f"\n  n_genes distribution:")
    print(f"    Min / Median / Max   : {ng['min']} / {ng['median']} / {ng['max']}")

    if EXOMISER_ERRORS.exists():
        errors = load_jsonl(EXOMISER_ERRORS)
        section["error_log"] = {"total_errors": len(errors)}
        print(f"\n  Error log: {len(errors)} records")
    else:
        section["error_log"] = {"note": "not found"}

    report["exomiser"] = section
    return ok


# ── ClinPrior audit ──────────────────────────────────────────────────────────

def audit_clinprior(report: dict) -> bool:
    """Audit ClinPrior rankings (phenotype-only mode)."""
    section = {}
    ok = True

    print(f"\n{'='*65}")
    print("  ClinPrior Rankings Audit (phenotype-only mode)")
    print(f"{'='*65}")

    if not CLINPRIOR_RANKINGS.exists():
        print(f"  [SKIP] {CLINPRIOR_RANKINGS} not found (run 08_run_clinprior.py first)")
        report["clinprior"] = {"skipped": "rankings file not found"}
        return True

    rankings = load_jsonl(CLINPRIOR_RANKINGS)
    section["total_lines"] = len(rankings)
    print(f"  Total lines parsed        : {len(rankings)}")

    pids = [r.get("patient_id") for r in rankings]
    pid_counts = Counter(pids)
    duplicates = {pid: cnt for pid, cnt in pid_counts.items() if cnt > 1}
    section["duplicate_pids"] = duplicates
    if duplicates:
        print(f"  [FAIL] Duplicate patient_ids: {list(duplicates.keys())[:10]}")
        ok = False
    else:
        print(f"  [PASS] No duplicate patient_ids")

    if MANIFEST_PATH.exists():
        manifest = load_jsonl(MANIFEST_PATH)
        manifest_pids = {r["patient_id"] for r in manifest}
        missing = manifest_pids - set(pids)
        section["missing_from_rankings"] = len(missing)
        if missing:
            print(f"  [WARN] {len(missing)} manifest PIDs not in ClinPrior rankings")
        else:
            print(f"  [PASS] All manifest PIDs present in ClinPrior rankings")

    ranks = [r.get("clinprior_rank") for r in rankings]
    rank_none = sum(1 for r in ranks if r is None)
    rank_found = [r for r in ranks if r is not None]
    section["rank_absent"] = rank_none
    section["rank_distribution"] = summarize_numeric(rank_found, "clinprior_rank")

    print(f"\n  Rank distribution:")
    print(f"    Absent (rank=None)   : {rank_none}  ({100*rank_none/len(ranks):.1f}%)")
    if rank_found:
        rs = summarize_numeric(rank_found, "clinprior_rank")
        print(f"    Min / Median / Max   : {rs['min']} / {rs['median']} / {rs['max']}")
        print(f"    Mean                 : {rs['mean']}")

    n_genes_vals = [r.get("clinprior_n_genes", 0) for r in rankings]
    section["n_genes_distribution"] = summarize_numeric(n_genes_vals, "clinprior_n_genes")
    ng = summarize_numeric(n_genes_vals, "clinprior_n_genes")
    print(f"\n  n_genes distribution:")
    print(f"    Min / Median / Max   : {ng['min']} / {ng['median']} / {ng['max']}")

    if CLINPRIOR_ERRORS.exists():
        errors = load_jsonl(CLINPRIOR_ERRORS)
        section["error_log"] = {"total_errors": len(errors)}
        print(f"\n  Error log: {len(errors)} records")
    else:
        section["error_log"] = {"note": "not found"}

    report["clinprior"] = section
    return ok


# ── cross-tool consistency check ──────────────────────────────────────────────

def audit_cross_tool(report: dict) -> None:
    """Cross-tool consistency: confirm shared PIDs and category counts."""
    section = {}

    print(f"\n{'='*65}")
    print("  Cross-Tool Consistency")
    print(f"{'='*65}")

    pid_sets = {}

    if MANIFEST_PATH.exists():
        manifest = load_jsonl(MANIFEST_PATH)
        pid_sets["manifest"] = {r["patient_id"] for r in manifest}
        print(f"  Manifest PIDs         : {len(pid_sets['manifest'])}")

    if LIRICAL_RANKINGS.exists():
        rankings = load_jsonl(LIRICAL_RANKINGS)
        pid_sets["lirical"] = {r["patient_id"] for r in rankings}
        print(f"  LIRICAL ranking PIDs  : {len(pid_sets['lirical'])}")

    if AIVA_RESULTS.exists():
        recs = load_jsonl(AIVA_RESULTS)
        pid_sets["aiva"] = {r["patient_id"] for r in recs}
        print(f"  AIVA PIDs             : {len(pid_sets['aiva'])}")

    if EXOMISER_RANKINGS.exists():
        rankings = load_jsonl(EXOMISER_RANKINGS)
        pid_sets["exomiser"] = {r["patient_id"] for r in rankings}
        print(f"  Exomiser ranking PIDs : {len(pid_sets['exomiser'])}")

    if CLINPRIOR_RANKINGS.exists():
        rankings = load_jsonl(CLINPRIOR_RANKINGS)
        pid_sets["clinprior"] = {r["patient_id"] for r in rankings}
        print(f"  ClinPrior ranking PIDs: {len(pid_sets['clinprior'])}")

    # Compute intersections
    active_sets = {k: v for k, v in pid_sets.items() if k != "manifest"}
    if len(active_sets) >= 2:
        shared = set.intersection(*active_sets.values())
        section["shared_pids_all_tools"] = len(shared)
        print(f"\n  Shared across all active tools: {len(shared)}")

    # Broad category breakdown (from LIRICAL rankings which has broad_category)
    if LIRICAL_RANKINGS.exists():
        rankings = load_jsonl(LIRICAL_RANKINGS)
        cat_counts = Counter(r.get("broad_category", "unknown") for r in rankings)
        section["category_counts"] = dict(cat_counts)
        print(f"\n  Broad category distribution (from LIRICAL rankings):")
        for cat in sorted(cat_counts.keys()):
            print(f"    {cat:<35} : {cat_counts[cat]}")

    report["cross_tool"] = section


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit all tool outputs before publication table generation"
    )
    parser.add_argument(
        "--include-exomiser", action="store_true",
        help="Include Exomiser rankings audit (requires exomiser_rankings.jsonl)",
    )
    parser.add_argument(
        "--include-clinprior", action="store_true",
        help="Include ClinPrior rankings audit (requires clinprior_rankings.jsonl)",
    )
    args = parser.parse_args()

    print("LIRICAL Benchmark — Validation Report")
    print(f"Output: {VALIDATION_REPORT}")

    report = {}
    all_ok = True

    all_ok &= audit_lirical(report)
    all_ok &= audit_aiva(report)
    if args.include_exomiser:
        all_ok &= audit_exomiser(report)
    if args.include_clinprior:
        all_ok &= audit_clinprior(report)
    audit_cross_tool(report)

    # ── write validation report ────────────────────────────────────────────
    report["validation_passed"] = all_ok
    VALIDATION_REPORT.parent.mkdir(parents=True, exist_ok=True)
    with open(VALIDATION_REPORT, "w") as f:
        json.dump(report, f, indent=2)

    print(f"\n{'='*65}")
    if all_ok:
        print("  VALIDATION PASSED — safe to proceed with publication metrics")
    else:
        print("  VALIDATION FAILED — see issues above before generating metrics")
    print(f"{'='*65}")
    print(f"\nReport saved: {VALIDATION_REPORT}")

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
