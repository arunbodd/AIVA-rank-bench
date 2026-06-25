#!/usr/bin/env python3
"""
04_compare_metrics.py
Compute and compare LIRICAL vs Exomiser vs ClinPrior vs AIVA performance on the benchmark.

Metrics per tool:
  - Top-1, Top-3, Top-5, Top-10, Top-20 accuracy (% of cases where true gene ranked ≤ k)
  - Mean Reciprocal Rank (MRR) — mean(1/rank); rank=None counts as 0; denominator = n_total
  - Gene Absent % (true gene not returned in output)
  - NDCG@5, NDCG@10, NDCG@20 (normalized discounted cumulative gain; IDCG=1 for single causal gene)
  - Median rank (found cases only)
  - Bootstrap 95% CIs on Top-1 and MRR (n=1000 iterations, seed=42)
  - McNemar's test p-values for pairwise Top-1/Top-5/Top-10 comparisons

Breakdown by broad_category.

Inputs:
  - results/lirical_rankings.jsonl      (from 02_run_lirical.py)
  - results/exomiser_rankings.jsonl     (from 06_run_exomiser.py)
  - results/clinprior_rankings.jsonl    (from 08_run_clinprior.py)
  - AIVA: Final_Broad_Cat_noDiseaseName_v2.jsonl
  - results/manifest.jsonl              (true gene symbols + categories)

Outputs:
  - results/comparison_metrics.csv
  - Printed summary table with CIs and McNemar p-values

Usage:
    python scripts/04_compare_metrics.py [--output-csv results/comparison_metrics.csv]
    python scripts/04_compare_metrics.py --no-exomiser     # exclude Exomiser
    python scripts/04_compare_metrics.py --no-clinprior    # exclude ClinPrior
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict
from pathlib import Path

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
# Configurable via env vars; defaults point at repo-bundled data. AIVA_RESULTS
# is the AIVA prediction JSONL you generate with scripts/predict_genes.py
# (bring your own API key) — see README.
AIVA_RESULTS = Path(os.getenv("AIVA_RESULTS", REPO_ROOT / "results" / "aiva_rankings.jsonl"))
INPUT_JSONL = Path(os.getenv("AIVA_BENCHMARK_INPUTS", REPO_ROOT / "inputs" / "benchmark_inputs.jsonl"))
LIRICAL_RANKINGS   = REPO_ROOT / "results" / "lirical_rankings.jsonl"
EXOMISER_RANKINGS  = REPO_ROOT / "results" / "exomiser_rankings.jsonl"
CLINPRIOR_RANKINGS = REPO_ROOT / "results" / "clinprior_rankings.jsonl"
MANIFEST_PATH      = REPO_ROOT / "inputs" / "manifest.jsonl"
DEFAULT_CSV        = REPO_ROOT / "results" / "comparison_metrics.csv"

TOP_K = [1, 3, 5, 10, 20]
TOOLS = ["LIRICAL", "Exomiser", "ClinPrior", "AIVA"]


# ── core metric helpers ───────────────────────────────────────────────────────

def compute_metrics(ranks: list) -> dict:
    """
    Compute metrics from a list of integer ranks (1-based) or None values.
    None = true gene absent from output.

    MRR denominator is n_total (absent genes contribute 0), following the
    conservative convention that penalises tools for not ranking the true gene.
    """
    n = len(ranks)
    if n == 0:
        return {}
    not_found = sum(1 for r in ranks if r is None)
    found_ranks = [r for r in ranks if r is not None]

    metrics = {
        "n": n,
        "gene_absent": not_found,
        "gene_absent_pct": round(100 * not_found / n, 2),
        "mrr": round(sum(1.0 / r for r in found_ranks) / n, 4) if found_ranks else 0.0,
    }
    for k in TOP_K:
        count = sum(1 for r in found_ranks if r <= k)
        metrics[f"top{k}"] = count
        metrics[f"top{k}_pct"] = round(100 * count / n, 2)
    return metrics


# ── publication-grade metric helpers ─────────────────────────────────────────

def compute_ndcg(ranks: list, k: int) -> float:
    """
    Mean NDCG@k over all cases.

    For single-label cases (one causal gene per patient):
      IDCG@k = 1.0  (ideal: true gene at rank 1 → 1/log2(2) = 1.0)
      Per-case DCG@k = 1/log2(rank+1) if rank ≤ k else 0
      Mean NDCG@k = mean(per-case DCG@k) over n_total cases

    MRR denominator convention applied: absent genes (rank=None) contribute 0.
    """
    n = len(ranks)
    if n == 0:
        return 0.0
    dcg_sum = sum(
        1.0 / math.log2(r + 1) for r in ranks if r is not None and r <= k
    )
    return round(dcg_sum / n, 4)


def compute_median_rank(ranks: list) -> float:
    """Median rank among found cases (rank is not None)."""
    found = sorted(r for r in ranks if r is not None)
    if not found:
        return float("nan")
    n = len(found)
    mid = n // 2
    return float(found[mid]) if n % 2 == 1 else (found[mid - 1] + found[mid]) / 2.0


def _mrr_fn(ranks: list) -> float:
    """MRR function for bootstrap sampling."""
    n = len(ranks)
    return sum(1.0 / r for r in ranks if r is not None) / n if n > 0 else 0.0


def _topk_pct_fn(k: int):
    """Returns a Top-k% function for bootstrap sampling."""
    def fn(ranks: list) -> float:
        n = len(ranks)
        return 100.0 * sum(1 for r in ranks if r is not None and r <= k) / n if n > 0 else 0.0
    return fn


def compute_bootstrap_ci(
    ranks: list,
    metric_fn,
    n_iter: int = 1000,
    seed: int = 42,
    ci: float = 0.95,
) -> tuple:
    """
    Bootstrap confidence interval for metric_fn applied to ranks list.

    Args:
        ranks:     List of integer ranks (1-based) or None values.
        metric_fn: Function mapping list → float (e.g. _mrr_fn).
        n_iter:    Number of bootstrap iterations (default 1000).
        seed:      Random seed for reproducibility (default 42).
        ci:        Confidence level (default 0.95 → 95% CI).

    Returns:
        (lo, hi) tuple rounded to 4 decimal places.
    """
    if not _HAS_NUMPY:
        return (float("nan"), float("nan"))
    rng = np.random.default_rng(seed)
    n = len(ranks)
    ranks_arr = np.array(ranks, dtype=object)
    boot_vals = []
    for _ in range(n_iter):
        idx = rng.integers(0, n, size=n)
        sample = [ranks_arr[i] for i in idx]
        boot_vals.append(metric_fn(sample))
    boot_vals.sort()
    # Use max(0, ...) and min(n_iter-1, ...) to guard against edge cases.
    # Subtract 1 from lo_idx so boot_vals[lo_idx] is the 2.5th percentile
    # (not the 2.6th, which int(0.025 * 1000) = 25 would give).
    lo_idx = max(0, round((1.0 - ci) / 2.0 * n_iter) - 1)
    hi_idx = min(n_iter - 1, round((1.0 + ci) / 2.0 * n_iter))
    return round(boot_vals[lo_idx], 4), round(boot_vals[hi_idx], 4)


def pairwise_mcnemar(ranks_a: list, ranks_b: list, k: int) -> float:
    """
    McNemar's test p-value for pairwise Top-k comparison between two tools.

    Uses chi-squared with continuity correction (Yates' correction).
    Discordant pairs: b_10 = A correct / B wrong; b_01 = B correct / A wrong.
    p-value computed from chi-squared distribution with df=1 via math.erfc.

    Bonferroni correction should be applied to the output p-values when
    comparing multiple k-values and tool pairs (see console output).
    """
    a_correct = [r is not None and r <= k for r in ranks_a]
    b_correct = [r is not None and r <= k for r in ranks_b]
    b_01 = sum(1 for a, b in zip(a_correct, b_correct) if not a and b)
    b_10 = sum(1 for a, b in zip(a_correct, b_correct) if a and not b)
    n_discordant = b_01 + b_10
    if n_discordant == 0:
        return 1.0
    # McNemar's chi-squared with continuity correction
    chi2_stat = (abs(b_10 - b_01) - 1.0) ** 2 / n_discordant
    # P(chi2 > x | df=1) = erfc(sqrt(x/2))
    p_val = math.erfc(math.sqrt(chi2_stat / 2.0))
    return round(float(p_val), 4)


# ── display helpers ───────────────────────────────────────────────────────────

def print_metrics_table(results: dict, title: str = "", cis: dict = None) -> None:
    """
    Pretty-print a comparison table for all tools.

    Args:
        results: {tool: metrics_dict}
        cis:     Optional {tool: {"top1_ci": (lo, hi), "mrr_ci": (lo, hi)}}
    """
    tools = [t for t in TOOLS if t in results and results[t]]
    if not tools:
        return
    cols = ["top1_pct", "top3_pct", "top5_pct", "top10_pct", "top20_pct",
            "mrr", "gene_absent_pct", "n"]
    col_names = ["Top-1 %", "Top-3 %", "Top-5 %", "Top-10 %", "Top-20 %",
                 "MRR", "Absent %", "N"]
    width = 11

    if title:
        print(f"\n{'─'*70}")
        print(f"  {title}")
        print(f"{'─'*70}")

    print(f"{'Tool':<12}", end="")
    for name in col_names:
        print(f"{name:>{width}}", end="")
    print()
    print("-" * (12 + width * len(cols)))

    for tool in tools:
        m = results[tool]
        print(f"{tool:<12}", end="")
        for col in cols:
            val = m.get(col, "-")
            print(f"{val:>{width}}", end="")
        print()
        # Print CI row if available
        if cis and tool in cis:
            tool_ci = cis[tool]
            top1_lo, top1_hi = tool_ci.get("top1_ci", (None, None))
            mrr_lo, mrr_hi = tool_ci.get("mrr_ci", (None, None))
            ci_parts = []
            if top1_lo is not None:
                ci_parts.append(f"  Top-1 95%CI [{top1_lo:.2f}, {top1_hi:.2f}]")
            if mrr_lo is not None:
                ci_parts.append(f"  MRR 95%CI [{mrr_lo:.4f}, {mrr_hi:.4f}]")
            if ci_parts:
                print(f"{'':12}{'  '.join(ci_parts)}")


def print_mcnemar_table(overall: dict, tools: list) -> None:
    """Print McNemar's test p-values for pairwise Top-k comparisons."""
    if len(tools) < 2:
        return

    k_vals = [1, 5, 10]
    # Number of comparisons for Bonferroni: C(n_tools, 2) × len(k_vals)
    n_pairs = len(tools) * (len(tools) - 1) // 2
    n_comparisons = n_pairs * len(k_vals)
    bonf_threshold = 0.05 / n_comparisons

    print(f"\n{'─'*70}")
    print(f"  McNemar's Test — Pairwise Top-k Comparisons (raw p-values)")
    print(f"  Bonferroni threshold: 0.05 / {n_comparisons} = {bonf_threshold:.4f}")
    print(f"  * = significant after Bonferroni correction")
    print(f"{'─'*70}")

    header = f"{'Comparison':<28}"
    for k in k_vals:
        header += f"{'Top-' + str(k) + ' p':>12}"
    print(header)
    print("-" * (28 + 12 * len(k_vals)))

    from itertools import combinations
    for tool_a, tool_b in combinations(tools, 2):
        row = f"{tool_a + ' vs ' + tool_b:<28}"
        for k in k_vals:
            p = pairwise_mcnemar(overall[tool_a], overall[tool_b], k)
            sig = "*" if p < bonf_threshold else " "
            row += f"{p:>11.4f}{sig}"
        print(row)


# ── data loaders ──────────────────────────────────────────────────────────────

def load_manifest(path: Path) -> dict:
    manifest = {}
    with open(path) as f:
        for line in f:
            rec = json.loads(line)
            manifest[rec["patient_id"]] = rec
    return manifest


def load_rankings(path: Path, rank_key: str) -> dict:
    """Load a rankings JSONL keyed by patient_id."""
    rankings = {}
    with open(path) as f:
        for line in f:
            rec = json.loads(line)
            rankings[rec["patient_id"]] = rec
    return rankings


def load_input_categories(path: Path) -> dict:
    cats = {}
    with open(path) as f:
        for line in f:
            rec = json.loads(line)
            cats[rec["id"]] = rec.get("broad_category", "unknown")
    return cats


def get_rank_from_list(gene_list: list, true_symbol: str, gene_key: str = "gene") -> int | None:
    """Find 1-based rank of true_symbol in a list of {gene_key: ..., score: ...} dicts."""
    if not true_symbol:
        return None
    true_upper = true_symbol.upper()
    for i, entry in enumerate(gene_list):
        gene = entry.get(gene_key, "")
        if gene == true_symbol or gene.upper() == true_upper:
            return i + 1
    return None


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare LIRICAL vs Exomiser vs ClinPrior vs AIVA metrics"
    )
    parser.add_argument(
        "--output-csv", type=Path, default=DEFAULT_CSV,
        help=f"Output CSV path (default: {DEFAULT_CSV})",
    )
    parser.add_argument(
        "--no-exomiser", action="store_true",
        help="Exclude Exomiser from comparison",
    )
    parser.add_argument(
        "--no-clinprior", action="store_true",
        help="Exclude ClinPrior from comparison",
    )
    parser.add_argument(
        "--bootstrap-iter", type=int, default=1000,
        help="Bootstrap iterations for CIs (default: 1000)",
    )
    args = parser.parse_args()

    # ── validate required files ───────────────────────────────────────────────
    missing = []
    if not LIRICAL_RANKINGS.exists():
        missing.append(f"{LIRICAL_RANKINGS}  → Run: python scripts/02_run_lirical.py")
    if not args.no_exomiser and not EXOMISER_RANKINGS.exists():
        missing.append(
            f"{EXOMISER_RANKINGS}  → Run: python scripts/06_run_exomiser.py\n"
            f"    (or use --no-exomiser to exclude Exomiser)"
        )
    if not args.no_clinprior and not CLINPRIOR_RANKINGS.exists():
        missing.append(
            f"{CLINPRIOR_RANKINGS}  → Run: python scripts/08_run_clinprior.py\n"
            f"    (or use --no-clinprior to exclude ClinPrior)"
        )
    if not AIVA_RESULTS.exists():
        missing.append(str(AIVA_RESULTS))
    if not MANIFEST_PATH.exists():
        missing.append(f"{MANIFEST_PATH}  → Run: python scripts/01_prepare_inputs.py")
    if missing:
        print("Missing required files:")
        for m in missing:
            print(f"  {m}")
        sys.exit(1)

    if not _HAS_NUMPY:
        print("[WARNING] numpy not found — bootstrap CIs will be skipped.")
        print("          Install: pip install numpy")

    active_tools = [
        t for t in TOOLS
        if not (t == "Exomiser" and args.no_exomiser)
        and not (t == "ClinPrior" and args.no_clinprior)
    ]

    print("Loading data...")
    manifest = load_manifest(MANIFEST_PATH)
    lirical = load_rankings(LIRICAL_RANKINGS, "lirical_rank")
    input_cats = load_input_categories(INPUT_JSONL)

    aiva_data = {}
    with open(AIVA_RESULTS) as f:
        for line in f:
            rec = json.loads(line)
            aiva_data[rec["patient_id"]] = rec

    exomiser_data = {}
    if not args.no_exomiser:
        exomiser_data = load_rankings(EXOMISER_RANKINGS, "exomiser_rank")

    clinprior_data = {}
    if not args.no_clinprior:
        clinprior_data = load_rankings(CLINPRIOR_RANKINGS, "clinprior_rank")

    # ── shared patient IDs ────────────────────────────────────────────────────
    shared_pids = set(lirical.keys()) & set(aiva_data.keys())
    if not args.no_exomiser:
        shared_pids &= set(exomiser_data.keys())
    if not args.no_clinprior:
        shared_pids &= set(clinprior_data.keys())

    print(f"  LIRICAL cases    : {len(lirical)}")
    if not args.no_exomiser:
        print(f"  Exomiser cases   : {len(exomiser_data)}")
    if not args.no_clinprior:
        print(f"  ClinPrior cases  : {len(clinprior_data)}")
    print(f"  AIVA cases       : {len(aiva_data)}")
    print(f"  Shared cases     : {len(shared_pids)}")

    if not shared_pids:
        print("ERROR: No shared cases across all tools.")
        sys.exit(1)

    # ── collect per-case ranks ────────────────────────────────────────────────
    overall: dict = {t: [] for t in active_tools}
    all_ranks: dict = defaultdict(lambda: {t: [] for t in active_tools})

    for pid in shared_pids:
        man = manifest.get(pid, {})
        true_sym = man.get("true_gene_symbol")
        category = input_cats.get(pid, "unknown")

        # LIRICAL rank (stored directly in record)
        lirical_rank = lirical[pid].get("lirical_rank")
        overall["LIRICAL"].append(lirical_rank)
        all_ranks[category]["LIRICAL"].append(lirical_rank)

        # Exomiser rank
        if not args.no_exomiser:
            exomiser_rank = exomiser_data[pid].get("exomiser_rank")
            # Fallback: recompute from top_genes list if rank is missing
            if exomiser_rank is None and true_sym:
                exomiser_rank = get_rank_from_list(
                    exomiser_data[pid].get("exomiser_top_genes", []), true_sym
                )
            overall["Exomiser"].append(exomiser_rank)
            all_ranks[category]["Exomiser"].append(exomiser_rank)

        # ClinPrior rank (stored directly; MT gene mapping already applied in R worker)
        if not args.no_clinprior:
            clinprior_rank = clinprior_data[pid].get("clinprior_rank")
            overall["ClinPrior"].append(clinprior_rank)
            all_ranks[category]["ClinPrior"].append(clinprior_rank)

        # AIVA rank — sort predicted_genes defensively before lookup
        # (1/1396 cases observed with non-sorted predicted_genes)
        sorted_aiva_genes = sorted(
            aiva_data[pid].get("predicted_genes", []),
            key=lambda x: x.get("score", 0), reverse=True,
        )
        aiva_rank = get_rank_from_list(sorted_aiva_genes, true_sym)
        overall["AIVA"].append(aiva_rank)
        all_ranks[category]["AIVA"].append(aiva_rank)

    # ── overall metrics ───────────────────────────────────────────────────────
    overall_metrics = {t: compute_metrics(overall[t]) for t in active_tools}

    # Bootstrap 95% CIs on Top-1 and MRR (overall only; per-category CIs
    # are less reliable at small n and are omitted from the table)
    overall_cis: dict = {}
    if _HAS_NUMPY:
        print(f"\nComputing bootstrap CIs (n={args.bootstrap_iter} iterations, seed=42)...")
        for tool in active_tools:
            top1_ci = compute_bootstrap_ci(
                overall[tool], _topk_pct_fn(1), n_iter=args.bootstrap_iter
            )
            mrr_ci = compute_bootstrap_ci(
                overall[tool], _mrr_fn, n_iter=args.bootstrap_iter
            )
            overall_cis[tool] = {"top1_ci": top1_ci, "mrr_ci": mrr_ci}

    print_metrics_table(overall_metrics,
                        f"Overall ({len(shared_pids)} shared cases)",
                        cis=overall_cis)

    # ── extended overall metrics (NDCG, median rank) ─────────────────────────
    print(f"\n{'─'*70}")
    print("  Extended Metrics (Overall)")
    print(f"{'─'*70}")
    ndcg_ks = [5, 10, 20]
    header = f"{'Tool':<12}"
    for k in ndcg_ks:
        header += f"{'NDCG@' + str(k):>12}"
    header += f"{'Med.Rank':>12}"
    print(header)
    print("-" * (12 + 12 * (len(ndcg_ks) + 1)))
    for tool in active_tools:
        row = f"{tool:<12}"
        for k in ndcg_ks:
            row += f"{compute_ndcg(overall[tool], k):>12.4f}"
        med = compute_median_rank(overall[tool])
        med_str = f"{med:.1f}" if not math.isnan(med) else "N/A"
        row += f"{med_str:>12}"
        print(row)

    # ── per-category metrics ──────────────────────────────────────────────────
    category_rows = []
    for cat in sorted(all_ranks.keys()):
        cat_metrics = {t: compute_metrics(all_ranks[cat][t]) for t in active_tools}
        n = cat_metrics.get("LIRICAL", {}).get("n", 0)
        print_metrics_table(cat_metrics, f"Category: {cat}  (n={n})")
        for tool, m in cat_metrics.items():
            row = {"category": cat, "tool": tool}
            row.update(m)
            # NDCG and median rank for category rows
            for k in ndcg_ks:
                row[f"ndcg{k}"] = compute_ndcg(all_ranks[cat][tool], k)
            med = compute_median_rank(all_ranks[cat][tool])
            row["median_rank"] = med if not math.isnan(med) else None
            # CIs for category rows (computed but not printed in table)
            if _HAS_NUMPY:
                top1_ci = compute_bootstrap_ci(
                    all_ranks[cat][tool], _topk_pct_fn(1), n_iter=args.bootstrap_iter
                )
                mrr_ci = compute_bootstrap_ci(
                    all_ranks[cat][tool], _mrr_fn, n_iter=args.bootstrap_iter
                )
                row["top1_ci_lo"], row["top1_ci_hi"] = top1_ci
                row["mrr_ci_lo"], row["mrr_ci_hi"] = mrr_ci
            category_rows.append(row)

    # ── McNemar's test (overall only) ────────────────────────────────────────
    if len(active_tools) >= 2:
        print_mcnemar_table(overall, active_tools)

    # ── write CSV ─────────────────────────────────────────────────────────────
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "category", "tool", "n",
        "top1", "top1_pct", "top3", "top3_pct",
        "top5", "top5_pct", "top10", "top10_pct", "top20", "top20_pct",
        "mrr", "gene_absent", "gene_absent_pct",
        "ndcg5", "ndcg10", "ndcg20", "median_rank",
        "top1_ci_lo", "top1_ci_hi", "mrr_ci_lo", "mrr_ci_hi",
    ]

    all_rows = []
    for tool in active_tools:
        m = overall_metrics[tool]
        row = {"category": "OVERALL", "tool": tool}
        row.update(m)
        for k in ndcg_ks:
            row[f"ndcg{k}"] = compute_ndcg(overall[tool], k)
        med = compute_median_rank(overall[tool])
        row["median_rank"] = med if not math.isnan(med) else None
        if tool in overall_cis:
            row["top1_ci_lo"], row["top1_ci_hi"] = overall_cis[tool]["top1_ci"]
            row["mrr_ci_lo"], row["mrr_ci_hi"] = overall_cis[tool]["mrr_ci"]
        all_rows.append(row)
    all_rows.extend(category_rows)

    with open(args.output_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    print(f"\nMetrics saved: {args.output_csv}")
    print("Next step: python scripts/05_plot_comparison.py")


if __name__ == "__main__":
    main()
