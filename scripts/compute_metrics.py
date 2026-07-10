#!/usr/bin/env python3
"""
AIVA Benchmark Metrics Computation
===================================
Reads prediction JSONL files and computes comprehensive evaluation metrics.

This script:
- Loads predictions from predict_genes_simple.py output
- Performs force-inclusion (appends missing candidate genes with score 0.0)
- Computes recall@K metrics (before and after force-inclusion)
- Computes ranking metrics (MRR, NDCG, mean/median rank)
- Handles per-category breakdowns
- Generates evaluation summary JSON

Key Metrics:
- retention_rate: Fraction where model returned true gene
- recall@K (before force): TRUE model performance
- recall@K (after force): Only counts genes with score > 0.0
- MRR, NDCG: Ranking quality (excludes score 0.0 genes)

Usage:
  python scripts/compute_metrics.py Results/hpo+vcf+genes+disease_predictions.jsonl
"""

import argparse
import json
import math
import re
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

# Strip trailing "(ENSGXXX)" parentheticals from AIVA gene output, e.g.
# "TEK (ENSG00000120156)" → "TEK"
_ID_IN_PARENS = re.compile(r'^(.+?)\s*\([^)]+\)\s*$')


# ──────────────────────────────────────────────────────────────
# Metric explanations
# ──────────────────────────────────────────────────────────────
SUMMARY_METRICS_EXPLANATION = {
    "n_patients": "Total number of patients evaluated in this run.",
    "no_prediction": "Patients where AIVA returned no parseable gene list (API failure or empty response). These count as missed for all Recall@K.",
    "retention_rate": "Fraction of patients where the model actually returned the true gene (BEFORE force-inclusion). This is the true model performance.",
    "filtered_out_rate": "Fraction of patients where the model missed the true gene and it was only added via force-inclusion.",
    "recall@1": "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene was the model's TOP-ranked gene. The strictest measure.",
    "recall@5": "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 5.",
    "recall@10": "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 10.",
    "recall@20": "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 20.",
    "recall@50": "TRUE MODEL PERFORMANCE (before force-inclusion): Fraction where the true gene appeared in the model's top 50.",
    "recall_after_force@1": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@5": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@10": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@20": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "recall_after_force@50": "SECONDARY METRIC (after force-inclusion): Only counts genes with score > 0.0 (model actually returned them). Genes with score = 0.0 are ignored.",
    "mrr": "Mean Reciprocal Rank (only counts genes with score > 0.0). Genes force-included with score 0.0 are treated as rank = None. 1.0 = always rank 1; 0.5 = average rank 2.",
    "ndcg@10": "Normalised Discounted Cumulative Gain at 10 (only counts genes with score > 0.0). Genes force-included with score 0.0 are treated as rank = None. 1.0 = perfect.",
    "mean_rank": "Average position among genes with score > 0.0 (model-assigned ranks only). Ignores arbitrary rankings of force-included genes.",
    "median_rank": "Median position among genes with score > 0.0 (model-assigned ranks only). Ignores arbitrary rankings of force-included genes.",
    "mean_rank_before_force": "Average position of the true gene in the model's ACTUAL output (before force-inclusion). Only includes cases where model found it.",
    "median_rank_before_force": "Median position of the true gene in the model's ACTUAL output (before force-inclusion). Only includes cases where model found it.",
}

SUMMARY_CATEGORY_DESCRIPTIONS = {
    "known_gene_disease_gene_in_kg": "ESTABLISHED CASE (baseline). Both the causal gene and its link to this disease are already documented in AIVA's knowledge graph. AIVA has direct evidence to draw on.",
    "known_gene_known_disease_gene_in_kg": "FULLY DOCUMENTED CASE. The gene and disease-gene association are both in the KG, and the gene is additionally connected via gene-gene relationships. Highest information density.",
    "known_gene_new_disease_gene_in_kg": "NOVEL DISEASE ASSOCIATION. The causal gene is well-known, but this specific gene-disease pairing is new — not yet in the disease-gene KG. Tests whether AIVA can reason from gene biology alone.",
    "new_gene_known_disease_gene_in_kg": "NOVEL GENE WITH KG SUPPORT. The causal gene is novel for this disease but does appear in the KG via some connection. Intermediate difficulty.",
    "new_gene_known_disease_gene_only_in_gg_kg": "NOVEL GENE (hardest). The causal gene is absent from AIVA's disease-gene KG and only reachable via gene-gene connections. Tests reasoning from indirect evidence and phenotype alone.",
}


# ──────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────
def rank_of_true_gene(predicted: list[str], true_genes: list[str]) -> int | None:
    """Return 1-based rank of first true gene in predicted list, or None.

    Normalises predicted strings so "TEK (ENSG00000120156)" matches "TEK".
    """
    if not predicted or not true_genes:
        return None

    true_upper = [g.upper() for g in true_genes]

    for i, pred in enumerate(predicted, start=1):
        m = _ID_IN_PARENS.match(pred.strip())
        pred_upper = (m.group(1).strip() if m else pred.strip()).upper()
        if pred_upper in true_upper:
            return i
    return None


def force_include_candidates(
    predicted: list[str],
    scores: dict[str, float],
    candidate_genes: list[str]
) -> tuple[list[str], dict[str, float]]:
    """Append missing candidate genes with score 0.0.

    Returns:
        (updated_predicted, updated_scores)
    """
    predicted_copy = list(predicted)
    scores_copy = dict(scores)

    predicted_upper = {g.upper() for g in predicted}

    for gene in candidate_genes:
        if gene.upper() not in predicted_upper:
            predicted_copy.append(gene)
            scores_copy[gene.upper()] = 0.0

    return predicted_copy, scores_copy


# ──────────────────────────────────────────────────────────────
# Metrics computation
# ──────────────────────────────────────────────────────────────
def compute_metrics(results: list[dict], at_ks: tuple[int, ...] = (1, 5, 10, 20, 50)) -> dict:
    """Compute comprehensive evaluation metrics.

    Args:
        results: List of prediction records from JSONL
        at_ks: Recall@K thresholds

    Returns:
        Dictionary with all metrics
    """
    n = len(results)
    if n == 0:
        return {}

    # Counters
    no_prediction = 0
    retention_count = 0
    found_at = {k: 0 for k in at_ks}
    found_at_after_force = {k: 0 for k in at_ks}

    # Ranking metrics (only for genes with score > 0.0)
    ranks = []
    ranks_before_force = []
    reciprocal_ranks = []

    # Per-category tracking
    per_category = defaultdict(list)
    per_category_before_force = defaultdict(list)

    for r in results:
        true_genes = r.get("true_genes", [])
        if not true_genes:
            continue

        # Original predictions (before force-inclusion)
        predicted_orig = r.get("predicted_genes", [])
        scores_orig = r.get("predicted_scores", {})

        # Force-include candidate genes
        candidate_genes = r.get("candidate_genes", [])
        predicted_force, scores_force = force_include_candidates(
            predicted_orig, scores_orig, candidate_genes
        )

        cat = r.get("category", "unknown")

        # ══════════════════════════════════════════════════════════
        # BEFORE FORCE-INCLUSION: True model performance
        # ══════════════════════════════════════════════════════════
        if not predicted_orig:
            no_prediction += 1
            rk_before = None
        else:
            rk_before = rank_of_true_gene(predicted_orig, true_genes)

        per_category_before_force[cat].append(rk_before)

        if rk_before is not None:
            retention_count += 1
            ranks_before_force.append(rk_before)
            for k in at_ks:
                if rk_before <= k:
                    found_at[k] += 1

        # ══════════════════════════════════════════════════════════
        # AFTER FORCE-INCLUSION: Ranking metrics (score > 0.0 only)
        # ══════════════════════════════════════════════════════════
        rk_after = rank_of_true_gene(predicted_force, true_genes)

        # Check if true gene has a real score (> 0.0)
        true_gene_has_real_score = False
        if rk_after is not None and true_genes:
            for tg in true_genes:
                gene_score = scores_force.get(tg.upper(), 0.0)
                if gene_score > 0.0:
                    true_gene_has_real_score = True
                    break

        per_category[cat].append(rk_after if true_gene_has_real_score else None)

        if rk_after is not None and true_gene_has_real_score:
            ranks.append(rk_after)
            reciprocal_ranks.append(1.0 / rk_after)
            for k in at_ks:
                if rk_after <= k:
                    found_at_after_force[k] += 1
        else:
            reciprocal_ranks.append(0.0)

    # ──────────────────────────────────────────────────────────────
    # Aggregate metrics
    # ──────────────────────────────────────────────────────────────
    retention_rate = retention_count / n if n > 0 else 0.0
    filtered_out_rate = 1.0 - retention_rate

    mrr = sum(reciprocal_ranks) / n if n > 0 else 0.0

    # NDCG@10
    ndcg_scores = []
    for rk in ranks:
        if rk <= 10:
            ndcg_scores.append(1.0 / math.log2(rk + 1))
        else:
            ndcg_scores.append(0.0)
    # Add zeros for cases where gene wasn't ranked with score > 0
    ndcg_scores.extend([0.0] * (n - len(ranks)))
    ndcg_10 = sum(ndcg_scores) / n if n > 0 else 0.0

    mean_rank = sum(ranks) / len(ranks) if ranks else None
    median_rank = sorted(ranks)[len(ranks) // 2] if ranks else None

    mean_rank_bf = sum(ranks_before_force) / len(ranks_before_force) if ranks_before_force else None
    median_rank_bf = sorted(ranks_before_force)[len(ranks_before_force) // 2] if ranks_before_force else None

    # ──────────────────────────────────────────────────────────────
    # Per-category metrics
    # ──────────────────────────────────────────────────────────────
    def category_recall(ranks_list: list, at_ks: tuple) -> dict:
        """Compute recall@K for a category."""
        n_cat = len(ranks_list)
        if n_cat == 0:
            return {}

        metrics = {"n": n_cat}
        for k in at_ks:
            count = sum(1 for r in ranks_list if r is not None and r <= k)
            metrics[f"recall@{k}"] = count / n_cat
        return metrics

    category_metrics = {}
    for cat, ranks_list in per_category_before_force.items():
        category_metrics[cat] = category_recall(ranks_list, at_ks)

    category_metrics_after = {}
    for cat, ranks_list in per_category.items():
        metrics = {"n": len(ranks_list)}
        for k in at_ks:
            count = sum(1 for r in ranks_list if r is not None and r <= k)
            metrics[f"recall_after_force@{k}"] = count / len(ranks_list) if len(ranks_list) > 0 else 0.0
        category_metrics_after[cat] = metrics

    # ──────────────────────────────────────────────────────────────
    # Return comprehensive metrics
    # ──────────────────────────────────────────────────────────────
    return {
        "n_patients": n,
        "no_prediction": no_prediction,
        "mrr": round(mrr, 4),
        "ndcg@10": round(ndcg_10, 4),
        "mean_rank": round(mean_rank, 2) if mean_rank is not None else None,
        "median_rank": median_rank,
        "mean_rank_before_force": round(mean_rank_bf, 2) if mean_rank_bf is not None else None,
        "median_rank_before_force": median_rank_bf,
        "retention_rate": round(retention_rate, 4),
        "filtered_out_rate": round(filtered_out_rate, 4),
        **{f"recall@{k}": round(found_at[k] / n, 4) for k in at_ks},
        **{f"recall_after_force@{k}": round(found_at_after_force[k] / n, 4) for k in at_ks},
        "per_category": category_metrics,
        "per_category_after_force": category_metrics_after,
    }


def print_metrics(metrics: dict, ablation: str = ""):
    """Pretty-print metrics to console."""
    n = metrics.get("n_patients", 0)
    if n == 0:
        print("No results to evaluate.")
        return

    print(f"\n{'═' * 80}")
    print(f"  EVALUATION SUMMARY: {ablation}")
    print(f"{'═' * 80}")
    print(f"  Total patients     : {n}")
    print(f"  No prediction      : {metrics.get('no_prediction', 0)} ({metrics.get('no_prediction', 0) / n:.1%})")
    print(f"  Retention rate     : {metrics.get('retention_rate', 0):.1%}  ← % where model returned true gene")
    print(f"  Filtered out       : {metrics.get('filtered_out_rate', 0):.1%}  ← % requiring force-inclusion")
    print()

    print(f"  BEFORE FORCE-INCLUSION (true model performance):")
    print(f"  ─────────────────────────────────────────────────────────────")
    for k in (1, 5, 10, 20, 50):
        key = f"recall@{k}"
        if key in metrics:
            print(f"    Recall@{k:<2}         : {metrics[key]:.1%}")

    mean_bf = metrics.get("mean_rank_before_force")
    median_bf = metrics.get("median_rank_before_force")
    if mean_bf or median_bf:
        print(f"    Mean rank          : {mean_bf}")
        print(f"    Median rank        : {median_bf}")

    print()
    print(f"  AFTER FORCE-INCLUSION (score > 0.0 only):")
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
# Main
# ──────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Compute metrics from AIVA prediction results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("predictions_file", type=str,
                   help="Path to predictions JSONL file")
    p.add_argument("--output", type=str, default=None,
                   help="Output summary JSON path (default: same name with _evaluation_summary.json suffix)")
    p.add_argument("--ablation", type=str, default="hpo+vcf+genes+disease",
                   help="Ablation label for summary")
    return p.parse_args()


def main():
    args = parse_args()

    predictions_path = Path(args.predictions_file)
    if not predictions_path.exists():
        sys.exit(f"File not found: {predictions_path}")

    # Load results
    print(f"Loading predictions from: {predictions_path}")
    results = []
    with open(predictions_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            try:
                results.append(json.loads(line))
            except json.JSONDecodeError:
                continue

    print(f"Loaded {len(results)} prediction records\n")

    # Compute metrics
    metrics = compute_metrics(results)

    # Print to console
    print_metrics(metrics, args.ablation)

    # Save summary JSON
    if args.output:
        summary_path = Path(args.output)
    else:
        summary_path = predictions_path.with_name(
            predictions_path.stem + "_evaluation_summary.json"
        )

    summary = {
        "ablation": args.ablation,
        "input": str(predictions_path),
        "metrics_explanation": SUMMARY_METRICS_EXPLANATION,
        "category_descriptions": SUMMARY_CATEGORY_DESCRIPTIONS,
        "metrics": metrics,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }

    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    print(f"Evaluation summary saved to: {summary_path}")


if __name__ == "__main__":
    main()
