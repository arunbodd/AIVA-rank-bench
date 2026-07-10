# AIVA-rank-bench

**AIVA** (AI Variant Analysis) is a clinical variant-interpretation tool developed by [Mamidi Health (MHSPL)](https://github.com/MHSPL). This repository runs **AIVA's gene-prioritization** over a cohort of 1,396 simulated rare-disease cases and scores it with standard ranking metrics. Precomputed rankings from three established tools — **LIRICAL**, **Exomiser**, and **ClinPrior** — ship alongside as reference baselines for comparison.

> **The main workflow is running AIVA** (steps below). The baseline tools are optional and only needed if you want to regenerate their rankings from scratch — the repo already includes their outputs.

---

## Quick start — run AIVA

```bash
git clone https://github.com/MHSPL/AIVA-rank-bech.git
cd AIVA-rank-bech
pip install -r requirements.txt

# Bring your own AIVA API key
cp .env.example .env        # then edit .env and set AIVA_API_KEY
#   (or: export AIVA_API_KEY=...)

# 1) Generate AIVA gene predictions over the 1,396-case benchmark
python scripts/predict_genes.py --ablation hpo+vcf
#    → results/hpo+vcf_aiva_predictions.jsonl

# 2) Score AIVA against the bundled baselines
export AIVA_RESULTS=results/hpo+vcf_aiva_predictions.jsonl
python scripts/04_compare_metrics.py     # → results/comparison_metrics.csv
python scripts/05_plot_comparison.py     # → results/figures/
```

That's the whole loop: **predict → score**. Inputs (HPO terms + per-case VCFs) and the baseline rankings are already in the repo, so AIVA is the only thing you need to run.

### Running AIVA in more detail

`scripts/predict_genes.py` queries the AIVA API for each case and writes a ranked gene list per patient.

```bash
python scripts/predict_genes.py \
    --ablation hpo+vcf \        # which evidence to send: hpo, vcf, hpo+vcf, ...
    --prompt-style aiva \       # prompt template (default: aiva)
    --input inputs/benchmark_inputs.jsonl   # default; the 1,396-case set
# Output: results/<ablation>_<prompt_style>_predictions.jsonl
```

Useful options: re-score an existing run without re-calling the API via
`--evaluate-only results/<file>.jsonl`; tune throughput with `CONCURRENCY_LIMIT`.

AIVA configuration (all via env / `.env`):

| Env var | Default | Purpose |
|---|---|---|
| `AIVA_API_KEY` | _(required)_ | AIVA API bearer token |
| `AIVA_API_URL` | `https://api.aivaportal.com/agents/chat/stream` | API endpoint |
| `AIVA_MODEL` | `google/gemini-3.1-pro-preview` | model id |
| `CONCURRENCY_LIMIT` | `3` | max simultaneous API calls |

---

## What's in this repo

A clone has everything needed to run AIVA and score it — only your API key is missing.

```
AIVA-rank-bench/
├── inputs/
│   ├── benchmark_inputs.jsonl      # 1,396-case benchmark (HPO terms, true gene) — AIVA input
│   ├── manifest.jsonl              # per-case HPO + true gene + VCF path (repo-relative)
│   └── vcfs/                       # 1,396 per-case VCFs (simulated_<id>_hg19.vcf)
├── results/
│   ├── lirical_rankings.jsonl      # baseline rankings (precomputed, for comparison)
│   ├── exomiser_rankings.jsonl
│   ├── clinprior_rankings.jsonl
│   ├── comparison_metrics.csv      # Top-k / MRR / nDCG summary
│   └── figures/                    # comparison plots
├── scripts/
├── requirements.txt
├── .env.example                    # AIVA API config template
└── README.md
```

**Bring your own:** an `AIVA_API_KEY` (the only thing required to run AIVA). **Not bundled:** the full 42k simulated set (~180 MB; only the 1,396-case subset ships) and the LIRICAL/Exomiser JARs + reference data (~6 GB; fetched on demand — see below).

---

## Requirements

To run AIVA you only need **Python ≥ 3.9** and an **AIVA API key**:

```bash
pip install -r requirements.txt   # httpx, numpy, matplotlib (+ optional tqdm, python-dotenv)
```

Java 17+ and R ≥ 4.2 (`ClinPrior`, `jsonlite`) are needed **only** if you choose to regenerate the baseline rankings.

---

## Optional — reproduce the baseline rankings

The baseline outputs are already committed. Regenerate them only if you want to rebuild the comparison from scratch (needs ~6 GB of tool data, Java, and — for ClinPrior — R).

```bash
bash scripts/00_download_setup.sh --setup all    # download LIRICAL + Exomiser data

python scripts/02_run_lirical.py                 # → results/lirical_rankings.jsonl
python scripts/03_parse_lirical_results.py
python scripts/06_run_exomiser.py                # → results/exomiser_rankings.jsonl
python scripts/08_run_clinprior.py               # → results/clinprior_rankings.jsonl
python scripts/07_validate_runs.py               # sanity-check the runs
```

Each runner supports `--limit N` / `--test` for a quick subset and resumes from existing output. Paths are env-overridable (`AIVA_BENCHMARK_INPUTS`, `AIVA_VCF_DIR`, `AIVA_RESULTS`, `AIVA_HGNC_PATH`); defaults point at the bundled data, so no configuration is needed.

---

## Scripts

| Script | Role |
|---|---|
| **`predict_genes.py`** | **Run AIVA — query the API for ranked gene predictions** |
| `04_compare_metrics.py` | Score AIVA vs baselines — Top-k / MRR / nDCG |
| `05_plot_comparison.py` | Comparison figures |
| `benchmark_gene_ranking.py`, `compute_metrics.py`, `analyze_distractors.py` | Supporting AIVA analysis (ranking eval, metrics, distractor breakdown) |
| `01_prepare_inputs.py` | Build `inputs/manifest.jsonl` from the benchmark inputs |
| `00_download_setup.sh` | Download LIRICAL 2.2.1 + Exomiser 14.0.0 (baseline only) |
| `02`/`03_*lirical*`, `06_run_exomiser.py`, `08_run_clinprior.py` (+ `_clinprior_worker.R`) | Run the baseline tools (baseline only) |
| `07_validate_runs.py` | Sanity-check runs (duplicates, coverage, error logs) |

---

## Benchmarking metrics

- **Top-1 / 3 / 5 / 10 / 20 accuracy** — causal gene ranked within top *N*
- **MRR** — mean reciprocal rank of the causal gene
- **nDCG** — normalized discounted cumulative gain
- **Rank distribution** — per-tool spread of causal-gene ranks

---

## Data provenance

Cases derive from a simulated rare-disease cohort (Jacobson et al. simulator). Each VCF carries CADD/PolyPhen2/MAF annotations in a simplified CSQ format. Reference ontologies: [HPO](https://hpo.jax.org/), [OMIM](https://www.omim.org/) disease–gene associations.

---

## License

Proprietary — maintained by [Mamidi Health (MHSPL)](https://github.com/MHSPL). Contact the maintainers for licensing inquiries. Contributions/issues welcome on the `dev` branch.
