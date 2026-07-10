#!/usr/bin/env python3
"""
05_plot_comparison.py
Generate comparison plots: LIRICAL vs Exomiser vs ClinPrior vs AIVA.
"""

import argparse
import json
import math
import os
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
except ImportError:
    print("ERROR: matplotlib and numpy required.")
    print("  Install: pip install matplotlib numpy")
    sys.exit(1)


# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
# Configurable via env vars; defaults point at repo-bundled data. AIVA_RESULTS
# is the AIVA prediction JSONL you generate with scripts/predict_genes.py.
AIVA_RESULTS = Path(os.getenv("AIVA_RESULTS", REPO_ROOT / "results" / "aiva_rankings.jsonl"))
INPUT_JSONL = Path(os.getenv("AIVA_BENCHMARK_INPUTS", REPO_ROOT / "inputs" / "benchmark_inputs.jsonl"))
LIRICAL_RANKINGS   = REPO_ROOT / "results" / "lirical_rankings.jsonl"
EXOMISER_RANKINGS  = REPO_ROOT / "results" / "exomiser_rankings.jsonl"
CLINPRIOR_RANKINGS = REPO_ROOT / "results" / "clinprior_rankings.jsonl"
MANIFEST_PATH      = REPO_ROOT / "inputs" / "manifest.jsonl"
FIGURES_DIR        = REPO_ROOT / "results" / "figures"

COLORS = {
    "LIRICAL":   "#0602FF",
    "Exomiser":  "#E07B00",
    "ClinPrior": "#2CA02C",
    "AIVA":      "#FF0D00",
}
TOP_K = [1, 3, 5, 10, 20]

# Explicit display labels for broad-category keys -> human-friendly strings
CATEGORY_LABEL_MAP = {
    "known_gene_disease": "Known Gene Disease",
    "known_gene_known_disease": "Known Gene + Known Disease",
    "known_gene_new_disease": "Known Gene + New Disease",
    "new_gene_known_disease": "New Gene + Known Disease",
    "new_gene_new_disease": "New Gene + New Disease",
}

# Fixed order for the 5 categories in the summary table
CATEGORY_ORDER = [
    "known_gene_disease",
    "known_gene_known_disease",
    "known_gene_new_disease",
    "new_gene_known_disease",
    "new_gene_new_disease",
]


# ── data loading ──────────────────────────────────────────────────────────────

def load_all_ranks(include_exomiser: bool = True, include_clinprior: bool = True):
    manifest = {}
    with open(MANIFEST_PATH) as f:
        for line in f:
            rec = json.loads(line)
            manifest[rec["patient_id"]] = rec

    lirical = {}
    with open(LIRICAL_RANKINGS) as f:
        for line in f:
            rec = json.loads(line)
            lirical[rec["patient_id"]] = rec

    exomiser = {}
    if include_exomiser and EXOMISER_RANKINGS.exists():
        with open(EXOMISER_RANKINGS) as f:
            for line in f:
                rec = json.loads(line)
                exomiser[rec["patient_id"]] = rec

    clinprior = {}
    if include_clinprior and CLINPRIOR_RANKINGS.exists():
        with open(CLINPRIOR_RANKINGS) as f:
            for line in f:
                rec = json.loads(line)
                clinprior[rec["patient_id"]] = rec

    input_cats = {}
    with open(INPUT_JSONL) as f:
        for line in f:
            rec = json.loads(line)
            input_cats[rec["id"]] = rec.get("broad_category", "unknown")

    aiva_data = {}
    with open(AIVA_RESULTS) as f:
        for line in f:
            rec = json.loads(line)
            aiva_data[rec["patient_id"]] = rec

    shared_pids = set(lirical.keys()) & set(aiva_data.keys())
    if include_exomiser and exomiser:
        shared_pids &= set(exomiser.keys())
    if include_clinprior and clinprior:
        shared_pids &= set(clinprior.keys())

    tools = ["LIRICAL"]
    if include_exomiser and exomiser:
        tools.append("Exomiser")
    if include_clinprior and clinprior:
        tools.append("ClinPrior")
    tools.append("AIVA")

    overall = {t: [] for t in tools}
    by_cat = defaultdict(lambda: {t: [] for t in tools})

    for pid in shared_pids:
        man = manifest.get(pid, {})
        true_sym = man.get("true_gene_symbol")
        category = input_cats.get(pid, "unknown")
        true_upper = true_sym.upper() if true_sym else ""

        # LIRICAL rank
        l_rank = lirical[pid].get("lirical_rank")
        overall["LIRICAL"].append(l_rank)
        by_cat[category]["LIRICAL"].append(l_rank)

        # Exomiser rank
        if "Exomiser" in tools:
            e_rec = exomiser.get(pid, {})
            e_rank = e_rec.get("exomiser_rank")
            if e_rank is None and true_sym:
                for i, entry in enumerate(e_rec.get("exomiser_top_genes", [])):
                    g = entry.get("gene", "")
                    if g == true_sym or g.upper() == true_upper:
                        e_rank = i + 1
                        break
            overall["Exomiser"].append(e_rank)
            by_cat[category]["Exomiser"].append(e_rank)

        # ClinPrior rank (already resolved for MT genes in R worker)
        if "ClinPrior" in tools:
            cp_rank = clinprior.get(pid, {}).get("clinprior_rank")
            overall["ClinPrior"].append(cp_rank)
            by_cat[category]["ClinPrior"].append(cp_rank)

        # AIVA rank
        a_rank = None
        if true_sym:
            sorted_aiva_genes = sorted(
                aiva_data[pid].get("predicted_genes", []),
                key=lambda x: x.get("score", 0), reverse=True,
            )
            for i, entry in enumerate(sorted_aiva_genes):
                g = entry.get("gene", "")
                if g == true_sym or g.upper() == true_upper:
                    a_rank = i + 1
                    break
        overall["AIVA"].append(a_rank)
        by_cat[category]["AIVA"].append(a_rank)

    return overall, dict(by_cat), tools


# ── metric helpers ────────────────────────────────────────────────────────────

def topk_pct(ranks, k):
    n = len(ranks)
    return 0.0 if n == 0 else 100 * sum(1 for r in ranks if r is not None and r <= k) / n


def topk_count(ranks, k):
    return sum(1 for r in ranks if r is not None and r <= k)


def mrr(ranks):
    n = len(ranks)
    return 0.0 if n == 0 else sum(1.0 / r for r in ranks if r is not None) / n


def compute_ndcg(ranks, k):
    n = len(ranks)
    if n == 0:
        return 0.0
    return sum(1.0 / math.log2(r + 1) for r in ranks if r is not None and r <= k) / n


def bootstrap_topk_ci(ranks, k, n_iter=500, seed=42):
    rng = np.random.default_rng(seed)
    n = len(ranks)
    ranks_arr = np.array(ranks, dtype=object)
    boot_vals = []
    for _ in range(n_iter):
        idx = rng.integers(0, n, size=n)
        sample = [ranks_arr[i] for i in idx]
        boot_vals.append(topk_pct(sample, k))
    boot_vals.sort()
    lo = boot_vals[max(0, round(0.025 * n_iter) - 1)]
    hi = boot_vals[min(n_iter - 1, round(0.975 * n_iter))]
    return lo, hi


def mcnemar_p(ranks_a, ranks_b, k):
    a_c = [r is not None and r <= k for r in ranks_a]
    b_c = [r is not None and r <= k for r in ranks_b]
    b_01 = sum(1 for a, b in zip(a_c, b_c) if not a and b)
    b_10 = sum(1 for a, b in zip(a_c, b_c) if a and not b)
    nd = b_01 + b_10
    if nd == 0:
        return 1.0
    chi2 = (abs(b_10 - b_01) - 1.0) ** 2 / nd
    return math.erfc(math.sqrt(chi2 / 2.0))


# ── plot 1: overall Top-k grouped bar ────────────────────────────────────────

def plot_overall_topk(overall: dict, tools: list, out_path: Path) -> None:
    n_total = len(overall[tools[0]])
    x = np.arange(len(TOP_K))
    x_labels = [f"Top-{k}\n(n={n_total})" for k in TOP_K]
    bar_width = 0.22
    fig, ax = plt.subplots(figsize=(11, 6))
    for i, tool in enumerate(tools):
        vals = [topk_pct(overall[tool], k) for k in TOP_K]
        counts = [topk_count(overall[tool], k) for k in TOP_K]
        ax.bar(
            x + i * bar_width - (bar_width * (len(tools) - 1) / 2),
            vals,
            width=bar_width,
            color=COLORS[tool],
            label=tool,
            alpha=0.85,
            edgecolor="black",
        )
        for xi, val, cnt in zip(x, vals, counts):
            ax.text(
                xi + i * bar_width - (bar_width * (len(tools) - 1) / 2),
                val + 2,
                f"{val:.1f}%\n({cnt}/{n_total})",
                ha="center",
                va="bottom",
                fontsize=8,
                color=COLORS[tool],
            )
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_ylabel("Cases with True Gene Ranked ≤ K (%)", fontsize=11)
    ax.set_title(
        f"{' vs '.join(tools)} — Overall Gene Ranking Performance\n"
        f"(n={n_total} cases)",
        fontsize=12,
    )
    ax.set_ylim(0, 110)
    handles, labels = ax.get_legend_handles_labels()
    label_order = [t for t in COLORS if t in tools]
    sorted_handles_labels = sorted(zip(labels, handles), key=lambda x: label_order.index(x[0]))
    sorted_labels, sorted_handles = zip(*sorted_handles_labels)
    ax.legend(sorted_handles, sorted_labels, fontsize=10, loc="lower right")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    png_path = out_path.with_suffix(".png")
    fig.savefig(png_path, dpi=150, format="png")
    plt.close(fig)


# ── plot 2: per-category grouped bars ────────────────────────────────────────

def plot_category_topk(by_cat: dict, tools: list, out_path: Path) -> None:
    tool_list = [t for t in ["LIRICAL", "Exomiser", "ClinPrior", "AIVA"] if t in tools]
    categories = sorted(by_cat.keys())
    n_tools = len(tool_list)
    bar_width = 0.22
    col_labels = [f"Top-{k}" for k in TOP_K]

    fig, axes = plt.subplots(1, len(TOP_K), figsize=(3.2 * len(TOP_K), 6), sharey=True)
    if len(TOP_K) == 1:
        axes = [axes]

    for j, (ax, k) in enumerate(zip(axes, TOP_K)):
        acc_by_tool_cat = {tool: [topk_pct(by_cat[cat][tool], k) for cat in categories] for tool in tool_list}
        sort_tool = tool_list[0]
        sorted_idx = np.argsort(acc_by_tool_cat[sort_tool])[::-1]
        sorted_categories = [categories[i] for i in sorted_idx]
        x = np.arange(len(sorted_categories))
        for i, tool in enumerate(tool_list):
            vals = [topk_pct(by_cat[cat][tool], k) for cat in sorted_categories]
            ax.bar(
                x + i * bar_width - (bar_width * (n_tools - 1) / 2),
                vals,
                width=bar_width,
                color=COLORS[tool],
                label=tool if j == 0 else None,
                alpha=0.85,
                edgecolor="black",
            )
        def format_cat_label(cat):
            return CATEGORY_LABEL_MAP.get(cat, cat.replace("_", " ").title())
        ax.set_xticks(x)
        ax.set_xticklabels([format_cat_label(cat) for cat in sorted_categories],
                           rotation=30, ha="right", fontsize=9)
        if j == 0:
            ax.set_ylabel("Top-K Accuracy (%)", fontsize=11)
        ax.set_title(col_labels[j], fontsize=11, pad=8)
        ax.set_ylim(0, 110)
        ax.grid(axis="y", linestyle="--", alpha=0.3)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=n_tools,
               fontsize=10, bbox_to_anchor=(0.5, 1.08))
    fig.suptitle("Top-K Accuracy by Broad Category", fontsize=13, y=1.12)
    fig.tight_layout(rect=[0, 0, 1, 1])
    fig.savefig(out_path, dpi=150, bbox_inches="tight", format="png")
    plt.close(fig)


# ── plot 3: rank distribution violin ─────────────────────────────────────────

def plot_rank_distribution(overall: dict, tools: list, out_path: Path) -> None:
    clip_at = 20
    positions = list(range(1, len(tools) + 1))
    fig, ax = plt.subplots(figsize=(3 * len(tools) + 1, 5))

    data = [[min(r, clip_at) for r in overall[t] if r is not None] for t in tools]
    parts = ax.violinplot(data, positions=positions, showmedians=True,
                          showextrema=True, widths=0.6)
    for pc, tool in zip(parts["bodies"], tools):
        pc.set_facecolor(COLORS[tool])
        pc.set_alpha(0.65)
    for key in ("cmedians", "cmins", "cmaxes", "cbars"):
        parts[key].set_color("black")
    parts["cmedians"].set_linewidth(2)

    rng = np.random.default_rng(42)
    for pos, tool, d in zip(positions, tools, data):
        jitter = rng.uniform(-0.08, 0.08, size=len(d))
        ax.scatter(np.full(len(d), pos) + jitter, d,
                   color=COLORS[tool], alpha=0.22, s=5, zorder=2)

    x_labels = [f"{t}\nfound: {len(d)}/{len(overall[t])}" for t, d in zip(tools, data)]
    ax.set_xticks(positions)
    ax.set_xticklabels(x_labels, fontsize=10)
    ax.set_ylabel(f"Rank of True Gene (capped at {clip_at})", fontsize=11)
    ax.set_title("Rank Distribution — Found Cases", fontsize=12)
    ax.set_ylim(0, clip_at + 1)
    ax.set_yticks(range(0, clip_at + 1, 5))
    ax.grid(axis="y", linestyle="--", alpha=0.4)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for pos, d in zip(positions, data):
        if d:
            med = float(np.median(d))
            mean = float(np.mean(d))
            ax.annotate(
                f"med={med:.0f}",
                xy=(pos, med), xytext=(pos + 0.35, med),
                fontsize=8.0, va="center",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.8),
            )
            ax.annotate(
                f"μ={mean:.1f}",
                xy=(pos, mean), xytext=(pos + 0.35, mean + 0.6),
                fontsize=7.5, va="center", color="dimgray",
                arrowprops=dict(arrowstyle="-", color="lightgray", lw=0.7),
            )
            ax.hlines(mean, pos - 0.28, pos + 0.28,
                      colors=["dimgray"], linewidths=1.5,
                      linestyles="dashed", zorder=4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, format="pdf")
    plt.close(fig)


# ── plot 4: MRR lollipop ─────────────────────────────────────────────────────

def plot_mrr(overall: dict, by_cat: dict, tools: list, out_path: Path) -> None:
    categories = ["OVERALL"] + sorted(by_cat.keys())
    y_labels = []
    for cat in categories:
        if cat == "OVERALL":
            n = len(overall[tools[0]])
            y_labels.append(f"OVERALL\n(n={n})")
        else:
            n = len(by_cat[cat].get(tools[0], []))
            y_labels.append(f"{cat.replace('_', ' ').title()}\n(n={n})")

    n_cats = len(categories)
    y_pos = np.arange(n_cats)
    n_tools = len(tools)
    offsets = np.linspace(0.22, -0.22, n_tools)

    fig, ax = plt.subplots(figsize=(11, 5))
    for tool, offset in zip(tools, offsets):
        vals = [
            mrr(overall[tool]) if cat == "OVERALL" else mrr(by_cat[cat].get(tool, []))
            for cat in categories
        ]
        y = y_pos + offset
        ax.hlines(y, 0, vals, color=COLORS[tool], linewidth=2, alpha=0.7)
        ax.scatter(vals, y, color=COLORS[tool], s=90, zorder=5, label=tool)
        for val, yi in zip(vals, y):
            ax.text(val + 0.012, yi, f"{val:.3f}",
                    va="center", ha="left", fontsize=7.5, color=COLORS[tool])

    ax.axhline(n_cats - 1.5, color="gray", linewidth=0.7, linestyle="--", alpha=0.4)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Mean Reciprocal Rank (MRR)", fontsize=11)
    tool_str = " vs ".join(tools)
    ax.set_title(f"{tool_str} — MRR by Category", fontsize=12)
    ax.set_xlim(0, 1.05)
    ax.legend(fontsize=10, loc="lower right")
    ax.grid(axis="x", linestyle="--", alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight", format="pdf")
    plt.close(fig)


# ── plot 5: NDCG grouped bars ────────────────────────────────────────────────

def plot_ndcg(overall: dict, tools: list, out_path: Path) -> None:
    ndcg_ks = [5, 10, 20]
    n_ks = len(ndcg_ks)
    n_tools = len(tools)

    ndcg_vals = {tool: [compute_ndcg(overall[tool], k) for k in ndcg_ks] for tool in tools}
    x = np.arange(n_ks)
    bar_width = 0.22
    offsets = np.linspace(-(n_tools - 1) / 2, (n_tools - 1) / 2, n_tools) * bar_width

    fig, ax = plt.subplots(figsize=(8, 5))
    for tool, offset in zip(tools, offsets):
        vals = ndcg_vals[tool]
        bars = ax.bar(
            x + offset,
            vals,
            width=bar_width * 0.9,
            label=tool,
            color=COLORS[tool],
            alpha=0.85,
            zorder=3,
        )
        for bar, val in zip(bars, vals):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.003,
                f"{val:.3f}",
                ha="center",
                va="bottom",
                fontsize=8.5,
                color=COLORS[tool],
                fontweight="bold",
            )

    ax.set_xticks(x)
    ax.set_xticklabels([f"NDCG@{k}" for k in ndcg_ks], fontsize=12)
    ax.set_ylabel("Mean NDCG (higher is better)", fontsize=11)
    tool_str = " vs ".join(tools)
    ax.set_title(
        f"{tool_str} — NDCG Comparison\n"
        f"(IDCG=1; single causal gene per case)",
        fontsize=12,
    )
    ax.set_ylim(0, min(1.0, max(v for vals in ndcg_vals.values() for v in vals) * 1.25))
    ax.legend(fontsize=10, loc="upper right")
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, format="pdf")
    plt.close(fig)


# ── plot 6+7: combined summary table ─────────────────────────────────────────

BAR_ORDER = ["AIVA", "LIRICAL", "ClinPrior", "Exomiser"]

def plot_combined_summary_table(
    overall: dict, by_cat: dict, tools: list, out_path: Path
) -> None:
    """
    Publication-style summary table.

    Fixes:
      - Category labels forced to the exact desired strings.
      - AIVA column highlight no longer overlaps the title.
      - Extra bottom padding so the footer is not cramped.
    """
    ordered = [t for t in BAR_ORDER if t in tools]
    n_total = len(overall[ordered[0]])

    # Overall rows (Top-1,5,10,20)
    overall_rows = [
        {
            "label": f"Top-{k}",
            "n": n_total,
            **{t: topk_pct(overall[t], k) for t in ordered},
        }
        for k in TOP_K
        if k != 3
    ]

    # Category rows in fixed order, using explicit display labels
    cat_rows = []
    for cat in CATEGORY_ORDER:
        if cat not in by_cat:
            continue
        label = CATEGORY_LABEL_MAP[cat]
        n_cat = len(by_cat[cat].get(ordered[0], []))
        row = {
            "label": label,
            "n": n_cat,
            **{t: topk_pct(by_cat[cat].get(t, []), 20) for t in ordered},
        }
        cat_rows.append(row)

    fig, ax = plt.subplots(figsize=(8.5 + 0.9 * max(0, len(ordered) - 3), 11))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    lx = 0.04
    # Dynamic column positions: space tools evenly in [0.54, 0.96]
    col_positions = np.linspace(0.54, 0.96, len(ordered))
    col_x = {tool: float(x) for tool, x in zip(ordered, col_positions)}
    row_h = 0.060
    cat_h = 0.078
    y = 0.955

    # AIVA column highlight (first in BAR_ORDER = highest priority)
    aiva_cx = col_x[ordered[0]]
    highlight_width = 0.13
    ax.add_patch(
        mpatches.FancyBboxPatch(
            (aiva_cx - highlight_width / 2, 0.090),  # keep same vertical offset
            highlight_width,
            0.85,
            boxstyle="round,pad=0.006",
            facecolor="#DBEAFE",
            edgecolor="#BFDBFE",
            linewidth=0.8,
            zorder=1,
            alpha=0.95,
        )
    )

    # Title
    ax.text(
        0.5,
        y,
        "Gene Prioritization Benchmark",
        ha="center",
        va="center",
        fontsize=13,
        color="#0F172A",
        fontweight="bold",
        zorder=2,
    )
    y -= row_h * 0.85

    # Column headers
    ax.text(
        lx,
        y,
        "Benchmark / Category",
        ha="left",
        va="center",
        fontsize=9.5,
        color="#334155",
        fontweight="bold",
        zorder=2,
    )
    for tool in ordered:
        hdr_color = "#1A56DB" if tool == "AIVA" else "#0F172A"
        ax.text(
            col_x[tool],
            y,
            tool,
            ha="center",
            va="center",
            fontsize=11,
            color=hdr_color,
            fontweight="bold",
            zorder=2,
        )
    y -= row_h * 0.42
    ax.axhline(y, xmin=0.02, xmax=0.98, color="#94A3B8", linewidth=0.9)
    y -= row_h * 0.58

    def zebra(y_pos: float, h: float, idx: int) -> None:
        if idx % 2 == 0:
            ax.add_patch(
                mpatches.Rectangle(
                    (0.025, y_pos - h * 0.52),
                    0.95,
                    h,
                    facecolor="#F8FAFC",
                    edgecolor="none",
                    zorder=0,
                )
            )

    def draw_overall_row(row_data: dict, y_pos: float, idx: int) -> None:
        zebra(y_pos, row_h, idx)
        vals = {t: row_data[t] for t in ordered}
        best_v = max(vals.values())

        ax.text(
            lx + 0.015,
            y_pos,
            row_data["label"],
            ha="left",
            va="center",
            fontsize=10,
            color="#1E293B",
            fontweight="bold",        # ← bold label
            zorder=2,
        )
        ax.text(
            lx + 0.12,
            y_pos,
            f"(n={row_data['n']:,})",
            ha="left",
            va="center",
            fontsize=8,
            color="#1E293B",          # darker
            fontweight="bold",        # ← bold n
            zorder=2,
        )
        for tool in ordered:
            v = vals[tool]
            ax.text(
                col_x[tool],
                y_pos,
                f"{v:.1f}%",
                ha="center",
                va="center",
                fontsize=10,
                color="#1E293B",
                fontweight="bold",    # ← always bold percentages
                zorder=2,
            )

    def draw_cat_row(row_data: dict, y_pos: float, idx: int) -> None:
        zebra(y_pos, cat_h, idx)
        vals = {t: row_data[t] for t in ordered}
        best_v = max(vals.values())

        # Label + n both bold on same row
        label_text = f"{row_data['label']}  (n={row_data['n']:,})"
        ax.text(
            lx + 0.015,
            y_pos,
            label_text,
            ha="left",
            va="center",
            fontsize=10,
            color="#1E293B",
            fontweight="bold",        # ← bold label + n
            zorder=2,
        )

        for tool in ordered:
            v = vals[tool]
            ax.text(
                col_x[tool],
                y_pos,
                f"{v:.1f}%",
                ha="center",
                va="center",
                fontsize=10,
                color="#1E293B",
                fontweight="bold",    # ← bold percentages
                zorder=2,
            )
    # Section 1: Overall
    ax.text(
        lx,
        y,
        f"Overall  ·  n = {n_total:,}",
        ha="left",
        va="center",
        fontsize=9,
        color="#64748B",
        fontweight="bold",
        fontstyle="italic",
        zorder=2,
    )
    y -= row_h
    for i, row in enumerate(overall_rows):
        draw_overall_row(row, y, i)
        y -= row_h

    y -= row_h * 0.25
    ax.axhline(y, xmin=0.02, xmax=0.98, color="#E2E8F0", linewidth=0.7)
    y -= row_h * 0.55

    # Section 2: Categories
    ax.text(
        lx,
        y,
        "Category Breakdown  ·  Top-20",
        ha="left",
        va="center",
        fontsize=9,
        color="#64748B",
        fontweight="bold",
        fontstyle="italic",
        zorder=2,
    )
    y -= row_h * 0.85
    for i, row in enumerate(cat_rows):
        draw_cat_row(row, y, i)
        y -= cat_h

    # Footer with extra padding
    y -= row_h * 0.25
    ax.axhline(y, xmin=0.02, xmax=0.98, color="#CBD5E1", linewidth=0.5)
    y -= row_h * 0.55
    ax.text(
        lx,
        y,
        "Category breakdown: Top-20 accuracy.  Bold = best per row.",
        ha="left",
        va="center",
        fontsize=7.5,
        color="#94A3B8",
        fontstyle="italic",
        zorder=2,
    )

    fig.savefig(out_path, dpi=150, bbox_inches="tight", format="pdf")
    plt.close(fig)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(description="Generate comparison plots")
    parser.add_argument(
        "--no-exomiser", action="store_true",
        help="Exclude Exomiser from plots",
    )
    parser.add_argument(
        "--no-clinprior", action="store_true",
        help="Exclude ClinPrior from plots",
    )
    args = parser.parse_args()

    for path, label in [
        (LIRICAL_RANKINGS, "LIRICAL rankings"),
        (AIVA_RESULTS, "AIVA results"),
        (MANIFEST_PATH, "Manifest"),
        (INPUT_JSONL, "Input JSONL"),
    ]:
        if not path.exists():
            print(f"ERROR: {label} not found: {path}")
            sys.exit(1)

    if not args.no_exomiser and not EXOMISER_RANKINGS.exists():
        print(f"WARNING: Exomiser rankings not found. Excluding Exomiser.")
        args.no_exomiser = True

    if not args.no_clinprior and not CLINPRIOR_RANKINGS.exists():
        print(f"WARNING: ClinPrior rankings not found. Excluding ClinPrior.")
        args.no_clinprior = True

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

    print("Loading data...")
    overall, by_cat, tools = load_all_ranks(
        include_exomiser=not args.no_exomiser,
        include_clinprior=not args.no_clinprior,
    )
    print(f"  {len(overall[tools[0]])} shared cases  |  tools: {', '.join(tools)}")

    print("\nGenerating plots...")
    plot_overall_topk(overall, tools, FIGURES_DIR / "overall_topk_line.pdf")
    plot_category_topk(by_cat, tools, FIGURES_DIR / "category_topk_heatmap.pdf")
    plot_rank_distribution(overall, tools, FIGURES_DIR / "rank_distribution.pdf")
    plot_mrr(overall, by_cat, tools, FIGURES_DIR / "mrr_comparison.pdf")
    plot_ndcg(overall, tools, FIGURES_DIR / "ndcg_comparison.pdf")
    plot_combined_summary_table(overall, by_cat, tools, FIGURES_DIR / "combined_summary_table.pdf")

    print(f"\nAll figures saved to: {FIGURES_DIR}/")


if __name__ == "__main__":
    main()