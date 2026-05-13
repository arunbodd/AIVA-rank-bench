# AIVA-rank-bench

**AIVA** (AI Variant Analysis) is a clinical variant interpretation benchmarking framework developed by [Mamidi Health (MHSPL)](https://github.com/MHSPL). This repository contains reproducible pipelines, scoring scripts, and result datasets used to evaluate AIVA's performance against established variant prioritization tools.

***

## Overview

AIVA-rank-bench benchmarks AIVA against state-of-the-art clinical variant interpretation tools including:

- **Exomiser** – phenotype-driven variant prioritization
- **LIRICAL** – likelihood ratio–based variant interpretation
- **InterVar** – ACMG/AMP evidence-based variant classification

The benchmark evaluates variant ranking accuracy, ACMG classification concordance, and rare disease diagnostic yield across curated patient cohorts.

***

## Repository Structure

```
AIVA-rank-bench/
├── scripts/            # Analysis and benchmarking scripts (Python, R)
├── results/            # Output files from benchmarking runs
├── data/               # Input VCF files, HPO phenotype terms (not tracked in git)
├── configs/            # Tool configuration files (Exomiser, LIRICAL YAML)
├── .gitignore
└── README.md
```

***

## Requirements

### System
- Linux / macOS
- Python ≥ 3.9
- R ≥ 4.2

### Python Dependencies
```bash
pip install pandas numpy scipy matplotlib seaborn
```

### R Dependencies
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr"))
```

### External Tools
| Tool | Version | Link |
|------|---------|------|
| Exomiser | 13.x | https://github.com/exomiser/Exomiser |
| LIRICAL | 2.x | https://github.com/TheJacksonLaboratory/LIRICAL |
| InterVar | 2.x | https://github.com/WGLab/InterVar |

***

## Usage

### 1. Clone the Repository
```bash
git clone https://github.com/MHSPL/AIVA-rank-bech.git
cd AIVA-rank-bech
```

### 2. Configure Tool Paths
Edit `configs/` YAML files to point to your local Exomiser/LIRICAL installations and reference data directories.

### 3. Run Benchmarking Scripts
```bash
# Run AIVA ranking pipeline
python scripts/run_aiva_benchmark.py --input data/cohort.vcf --hpo data/hpo_terms.txt

# Compare against Exomiser
python scripts/compare_tools.py --aiva results/aiva_output.tsv --exomiser results/exomiser_output.tsv
```

### 4. Analyze Results
```r
# In R
source("scripts/plot_benchmark_results.R")
```

***

## Benchmarking Metrics

The following metrics are computed for each tool:

- **Top-1 / Top-5 / Top-10 accuracy** – fraction of cases where causal gene is ranked within top N
- **MRR (Mean Reciprocal Rank)** – average reciprocal rank of causal gene
- **ACMG classification concordance** – agreement rate with curated ACMG classifications
- **Recall @ threshold** – sensitivity at various rank cutoffs

***

## Data

Input data (VCF files, reference genomes, phenotype datasets) are **not tracked in this repository** due to size and privacy constraints. Reference datasets used include:

- [OMIM](https://www.omim.org/) disease–gene associations
- [HPO](https://hpo.jax.org/) Human Phenotype Ontology terms
- Simulated/de-identified rare disease exome cohorts

***

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or submit a pull request on the `dev` branch.

```
dev    ← active development branch
main   ← stable releases
```

***

## License

This project is proprietary and maintained by [Mamidi Health (MHSPL)](https://github.com/MHSPL). Contact the repository maintainers for licensing inquiries.

***

## Contact

For questions or collaboration inquiries, open a GitHub issue or contact the MHSPL team via the organization page.
