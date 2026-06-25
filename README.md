# AIVA-rank-bench

**AIVA** (AI Variant Analysis) is a clinical variant-interpretation tool developed by [Mamidi Health (MHSPL)](https://github.com/MHSPL). This repository benchmarks AIVA's **gene-prioritization** performance against three established phenotype/variant prioritization tools — **LIRICAL**, **Exomiser**, and **ClinPrior** — on a cohort of 1,396 simulated rare-disease cases.

---

## What's in this repo

A clone is self-contained for the baseline tools: the case manifest, the per-case VCFs, and the precomputed rankings/figures are all bundled.

```
AIVA-rank-bench/
├── inputs/
│   ├── manifest.jsonl              # 1,396 cases: HPO terms, true gene, VCF path (repo-relative)
│   ├── benchmark_inputs.jsonl      # raw benchmark input records (source for the manifest)
│   └── vcfs/                       # 1,396 per-case VCFs (simulated_<id>_hg19.vcf)
├── results/
│   ├── lirical_rankings.jsonl      # precomputed baseline rankings
│   ├── exomiser_rankings.jsonl
│   ├── clinprior_rankings.jsonl
│   ├── comparison_metrics.csv      # Top-k / MRR / nDCG summary
│   ├── validation_report.json
│   └── figures/                    # Top-k, MRR, nDCG, rank-distribution plots
├── scripts/                        # pipeline (see below)
├── requirements.txt
├── .env.example                    # AIVA API config template (bring your own key)
└── README.md
```

**Not bundled** (size / proprietary / bring-your-own):
- The full 42k simulated set (`final_extracted.jsonl`, ~180 MB) — only the 1,396-case benchmark subset ships here.
- **AIVA predictions** — you generate these yourself with your own API key (see [Running AIVA](#4-run-aiva-bring-your-own-api)).
- LIRICAL/Exomiser JARs + reference data (~6 GB) — fetched by `scripts/00_download_setup.sh` into `tools/`.

---

## Requirements

| | |
|---|---|
| **OS** | Linux / macOS |
| **Python** | ≥ 3.9 — `pip install -r requirements.txt` |
| **Java** | 17+ (LIRICAL, Exomiser) |
| **R** | ≥ 4.2 with `ClinPrior` and `jsonlite` (ClinPrior step only) |
| **Disk** | ~6 GB for LIRICAL + Exomiser data |

Python deps: `httpx`, `numpy`, `matplotlib` (core); `tqdm`, `python-dotenv` (optional). See [requirements.txt](requirements.txt).

---

## Configuration

All external paths are environment-overridable; defaults point at the bundled data, so **no configuration is needed for the baseline tools**. Override only if your data lives elsewhere:

| Env var | Default | Purpose |
|---|---|---|
| `AIVA_API_KEY` | _(unset)_ | AIVA API bearer token (required to run AIVA) |
| `AIVA_BENCHMARK_INPUTS` | `inputs/benchmark_inputs.jsonl` | benchmark input records |
| `AIVA_VCF_DIR` | `inputs/vcfs` | per-case VCFs |
| `AIVA_RESULTS` | `results/aiva_rankings.jsonl` | AIVA predictions (you generate) |
| `AIVA_HGNC_PATH` | `data/hgnc_complete_set.txt` | optional HGNC map (graceful if absent) |

Copy `.env.example` → `.env` and fill in `AIVA_API_KEY`.

---

## Usage

```bash
git clone https://github.com/MHSPL/AIVA-rank-bech.git
cd AIVA-rank-bech
pip install -r requirements.txt
```

### 1. Download tools (LIRICAL + Exomiser)
```bash
bash scripts/00_download_setup.sh --setup all      # or: lirical | exomiser
```

### 2. (Re)build the manifest — optional
The manifest is already committed. Rebuild it only if you change the input set:
```bash
python scripts/01_prepare_inputs.py
```

### 3. Run the baseline tools
```bash
python scripts/02_run_lirical.py            # → results/lirical_rankings.jsonl
python scripts/03_parse_lirical_results.py
python scripts/06_run_exomiser.py           # → results/exomiser_rankings.jsonl
python scripts/08_run_clinprior.py          # → results/clinprior_rankings.jsonl  (needs R + ClinPrior)
```
Each runner supports `--limit N` / `--test` for a quick subset and resumes from existing output.

### 4. Run AIVA (bring your own API)
```bash
export AIVA_API_KEY=...                      # or set it in .env
python scripts/predict_genes.py --ablation hpo+vcf
# → results/<ablation>_<prompt_style>_predictions.jsonl
# Point AIVA_RESULTS at this file for the comparison step.
```

### 5. Compare, plot, validate
```bash
python scripts/04_compare_metrics.py        # → results/comparison_metrics.csv
python scripts/05_plot_comparison.py        # → results/figures/
python scripts/07_validate_runs.py          # → results/validation_report.json
```

---

## Scripts

| Script | Role |
|---|---|
| `00_download_setup.sh` | Download LIRICAL 2.2.1 + Exomiser 14.0.0 and reference data |
| `01_prepare_inputs.py` | Build `inputs/manifest.jsonl` (HPO terms, true gene, VCF path) |
| `02_run_lirical.py` / `03_parse_lirical_results.py` | Run + parse LIRICAL |
| `06_run_exomiser.py` | Run Exomiser |
| `08_run_clinprior.py` + `_clinprior_worker.R` | Run ClinPrior (HPO-only) |
| `predict_genes.py` | Query the AIVA API for gene predictions |
| `04_compare_metrics.py` | Top-k / MRR / nDCG across all tools |
| `05_plot_comparison.py` | Comparison figures |
| `07_validate_runs.py` | Sanity-check runs (duplicates, coverage, error logs) |
| `benchmark_gene_ranking.py`, `compute_metrics.py`, `analyze_distractors.py` | Supporting analysis (ranking eval, metrics, distractor breakdown) |

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
