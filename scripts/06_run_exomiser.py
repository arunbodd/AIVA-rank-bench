#!/usr/bin/env python3
"""
06_run_exomiser.py
Run Exomiser 14.x on every case in the manifest (results/manifest.jsonl).

For each case:
  1. Generate a per-case analysis YAML (HPO terms + VCF path + output prefix)
  2. Run Exomiser via Java with a shared application.properties
  3. Parse TSV_GENE output inline for gene rankings
  4. Append compact record to results/exomiser_rankings.jsonl (thread-safe)
  5. Delete per-case output files to save disk space

Features:
  - Parallel execution (--workers N, default 2; each JVM ~6 GB RAM)
  - Resume: skips cases already in exomiser_rankings.jsonl
  - Progress bar (tqdm if installed)
  - Errors logged to results/exomiser_errors.jsonl
  - Dry-run mode (--dry-run)

Usage:
    python scripts/06_run_exomiser.py [--workers 2] [--limit N] [--dry-run]
    python scripts/06_run_exomiser.py --test    # first 5 cases only
"""

import argparse
import json
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
EXOMISER_VERSION = "14.0.0"
EXOMISER_JAR = (
    REPO_ROOT / "tools" / f"exomiser-cli-{EXOMISER_VERSION}"
    / f"exomiser-cli-{EXOMISER_VERSION}.jar"
)
EXOMISER_DATA_DIR = REPO_ROOT / "tools" / "exomiser_data"
MANIFEST_PATH = REPO_ROOT / "inputs" / "manifest.jsonl"
OUTPUTS_DIR = REPO_ROOT / "outputs" / "exomiser"
RANKINGS_PATH = REPO_ROOT / "results" / "exomiser_rankings.jsonl"
ERRORS_PATH = REPO_ROOT / "results" / "exomiser_errors.jsonl"
APP_PROPS_PATH = OUTPUTS_DIR / "application.properties"

DATA_RELEASE = "2406"
TOP_N = 20  # max genes to store per case

# ── thread-safe writer ────────────────────────────────────────────────────────
_write_lock = threading.Lock()


def append_ranking(record: dict) -> None:
    """Thread-safe append of one compact ranking record."""
    with _write_lock:
        with open(RANKINGS_PATH, "a") as f:
            f.write(json.dumps(record) + "\n")


# ── YAML generation ───────────────────────────────────────────────────────────

# Use .format() placeholders; literal {} in YAML must be doubled as {{}}
_YAML_TEMPLATE = """\
analysis:
  genomeAssembly: hg19
  vcf: '{vcf}'
  proband: ''
  hpoIds:
{hpo_block}
  inheritanceModes:
    AUTOSOMAL_DOMINANT: 0.1
    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0
    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1
    X_DOMINANT: 0.1
    X_RECESSIVE_COMP_HET: 2.0
    X_RECESSIVE_HOM_ALT: 0.1
    MITOCHONDRIAL: 0.2
  steps:
    - variantEffectFilter:
        remove:
          - FIVE_PRIME_UTR_EXON_VARIANT
          - FIVE_PRIME_UTR_TRUNCATION
          - THREE_PRIME_UTR_EXON_VARIANT
          - THREE_PRIME_UTR_TRUNCATION
          - NON_CODING_TRANSCRIPT_EXON_VARIANT
          - UPSTREAM_GENE_VARIANT
          - INTERGENIC_VARIANT
          - REGULATORY_REGION_VARIANT
          - CODING_TRANSCRIPT_INTRON_VARIANT
          - NON_CODING_TRANSCRIPT_INTRON_VARIANT
          - DOWNSTREAM_GENE_VARIANT
          - INTRON_VARIANT
          - SYNONYMOUS_VARIANT
    - inheritanceFilter: {{}}
    - omimPrioritiser: {{}}
    - hiPhivePrioritiser: {{}}
outputOptions:
  outputContributingVariantsOnly: false
  numGenes: 0
  outputPrefix: '{output_prefix}'
  outputFormats:
    - TSV_GENE
"""

# Human-only variant: hiPhivePrioritiser restricted to human phenotype associations.
# Genome Medicine 2025 (DOI:10.1186/s13073-025-01546-1) shows human-only can
# improve ranking on some cohorts; used as an ablation run (--human-only flag).
_YAML_TEMPLATE_HUMAN = """\
analysis:
  genomeAssembly: hg19
  vcf: '{vcf}'
  proband: ''
  hpoIds:
{hpo_block}
  inheritanceModes:
    AUTOSOMAL_DOMINANT: 0.1
    AUTOSOMAL_RECESSIVE_COMP_HET: 2.0
    AUTOSOMAL_RECESSIVE_HOM_ALT: 0.1
    X_DOMINANT: 0.1
    X_RECESSIVE_COMP_HET: 2.0
    X_RECESSIVE_HOM_ALT: 0.1
    MITOCHONDRIAL: 0.2
  steps:
    - variantEffectFilter:
        remove:
          - FIVE_PRIME_UTR_EXON_VARIANT
          - FIVE_PRIME_UTR_TRUNCATION
          - THREE_PRIME_UTR_EXON_VARIANT
          - THREE_PRIME_UTR_TRUNCATION
          - NON_CODING_TRANSCRIPT_EXON_VARIANT
          - UPSTREAM_GENE_VARIANT
          - INTERGENIC_VARIANT
          - REGULATORY_REGION_VARIANT
          - CODING_TRANSCRIPT_INTRON_VARIANT
          - NON_CODING_TRANSCRIPT_INTRON_VARIANT
          - DOWNSTREAM_GENE_VARIANT
          - INTRON_VARIANT
          - SYNONYMOUS_VARIANT
    - inheritanceFilter: {{}}
    - omimPrioritiser: {{}}
    - hiPhivePrioritiser:
        runParams: 'human'
outputOptions:
  outputContributingVariantsOnly: false
  numGenes: 0
  outputPrefix: '{output_prefix}'
  outputFormats:
    - TSV_GENE
"""


def write_analysis_yaml(
    yaml_path: Path, vcf_path: str, hpo_ids: list, output_prefix: str,
    human_only: bool = False,
) -> None:
    """Write a per-case Exomiser analysis YAML.

    Args:
        human_only: If True, restricts hiPhivePrioritiser to human associations only.
                    Default (False) uses human+mouse+fish (standard Exomiser config).
    """
    hpo_block = "\n".join(f"    - '{hid}'" for hid in hpo_ids)
    template = _YAML_TEMPLATE_HUMAN if human_only else _YAML_TEMPLATE
    yaml_path.write_text(
        template.format(vcf=vcf_path, hpo_block=hpo_block, output_prefix=output_prefix)
    )


def write_application_properties(props_path: Path) -> None:
    """Write the shared Exomiser application.properties (one-time, shared across cases)."""
    props_path.write_text(
        f"exomiser.data-directory={EXOMISER_DATA_DIR}\n"
        f"exomiser.hg19.data-version={DATA_RELEASE}\n"
        f"exomiser.phenotype.data-version={DATA_RELEASE}\n"
        f"spring.jpa.open-in-view=false\n"
    )


# ── Exomiser TSV parsing ──────────────────────────────────────────────────────

def find_exomiser_tsv(out_dir: Path, pid: str) -> Path | None:
    """
    Locate the Exomiser gene TSV, handling version naming differences.
    Exomiser 14+: {prefix}.genes.tsv
    Exomiser 12-13: {prefix}-EXOMISER-genes.tsv
    """
    for candidate in [
        out_dir / f"{pid}.genes.tsv",
        out_dir / f"{pid}-EXOMISER-genes.tsv",
    ]:
        if candidate.exists():
            return candidate
    return None


def parse_exomiser_tsv(tsv_path: Path) -> list:
    """
    Parse Exomiser TSV_GENE output.

    Returns [{gene, score}, ...] sorted by EXOMISER_GENE_COMBINED_SCORE descending.
    Deduplicates by gene symbol (one gene may appear multiple times for different
    diseases — keeps the max combined score per gene).

    TSV format (Exomiser 14):
      #GENE_SYMBOL  ENTREZ_GENE_ID  EXOMISER_GENE_COMBINED_SCORE  ...
      (header line starts with '#')
    """
    gene_best: dict = {}

    with open(tsv_path) as f:
        header = None
        gene_col = score_col = None

        for raw_line in f:
            line = raw_line.rstrip("\n")
            if not line.strip():
                continue

            if line.startswith("#"):
                # Strip leading '#' and parse header
                cols = line.lstrip("#").split("\t")
                header = [c.strip() for c in cols]
                # Locate key columns by name
                gene_col = next(
                    (i for i, h in enumerate(header) if h == "GENE_SYMBOL"), None
                )
                score_col = next(
                    (i for i, h in enumerate(header) if h == "EXOMISER_GENE_COMBINED_SCORE"),
                    None,
                )
                # Fallback: partial match
                if gene_col is None:
                    gene_col = next(
                        (i for i, h in enumerate(header) if "GENE_SYMBOL" in h), None
                    )
                if score_col is None:
                    score_col = next(
                        (i for i, h in enumerate(header) if "COMBINED_SCORE" in h), None
                    )
                continue

            if header is None or gene_col is None or score_col is None:
                continue

            parts = line.split("\t")
            if len(parts) <= max(gene_col, score_col):
                continue

            try:
                gene = parts[gene_col].strip()
                score = float(parts[score_col].strip())
            except (ValueError, IndexError):
                continue

            if gene and (gene not in gene_best or score > gene_best[gene]):
                gene_best[gene] = score

    ranked = sorted(gene_best.items(), key=lambda x: x[1], reverse=True)
    return [{"gene": sym, "score": round(score, 10)} for sym, score in ranked]


def find_rank(ranked_genes: list, true_symbol: str) -> int | None:
    """Return 1-based rank of true_symbol, or None if absent. Case-insensitive fallback."""
    if not true_symbol:
        return None
    true_upper = true_symbol.upper()
    for i, entry in enumerate(ranked_genes):
        g = entry["gene"]
        if g == true_symbol or g.upper() == true_upper:
            return i + 1
    return None


# ── per-case runner ───────────────────────────────────────────────────────────

def run_case(entry: dict, done_pids: set, dry_run: bool = False,
             human_only: bool = False) -> dict:
    """
    Run Exomiser for one case, parse inline, write ranking, delete outputs.
    Returns a status dict.

    Args:
        human_only: If True, restricts hiPhivePrioritiser to human associations.
    """
    pid = entry["patient_id"]
    t_start = time.time()
    result = {"patient_id": pid, "status": None, "elapsed_s": None, "error": None}

    if pid in done_pids:
        result["status"] = "skipped"
        return result

    if not entry.get("vcf_exists", False):
        result["status"] = "error"
        result["error"] = f"VCF not found: {entry['vcf_path']}"
        return result

    observed = entry.get("observed_hpo", [])
    if not observed:
        result["status"] = "error"
        result["error"] = "No observed HPO terms"
        return result

    out_dir = OUTPUTS_DIR / str(pid)
    out_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = out_dir / "analysis.yml"
    output_prefix = str(out_dir / str(pid))

    write_analysis_yaml(yaml_path, entry["vcf_path"], observed, output_prefix,
                        human_only=human_only)

    cmd = [
        "java", "-Xmx6g",
        f"-Dspring.config.location=file:{APP_PROPS_PATH}",
        "-jar", str(EXOMISER_JAR),
        "--analysis", str(yaml_path),
    ]

    if dry_run:
        result["status"] = "dry_run"
        result["cmd"] = " ".join(cmd)
        return result

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 min per case
        )
        elapsed = time.time() - t_start
        result["elapsed_s"] = round(elapsed, 1)

        if proc.returncode != 0:
            result["status"] = "error"
            stderr_tail = proc.stderr[-1200:] if proc.stderr else ""
            stdout_tail = proc.stdout[-1200:] if proc.stdout else ""
            if stderr_tail and stdout_tail:
                result["error"] = f"stderr tail:\n{stderr_tail}\n\nstdout tail:\n{stdout_tail}"
            elif stderr_tail:
                result["error"] = stderr_tail
            elif stdout_tail:
                result["error"] = stdout_tail
            else:
                result["error"] = "non-zero exit"
            return result

        tsv_path = find_exomiser_tsv(out_dir, str(pid))
        if tsv_path is None:
            result["status"] = "error"
            result["error"] = (
                f"TSV output not found in {out_dir}. "
                f"stdout tail: {proc.stdout[-300:]}"
            )
            return result

        try:
            ranked_genes = parse_exomiser_tsv(tsv_path)
        except Exception as e:
            result["status"] = "error"
            result["error"] = f"TSV parse error: {e}"
            return result
        finally:
            # Clean up per-case files to save disk space
            for p in [tsv_path, yaml_path]:
                try:
                    if p and p.exists():
                        p.unlink()
                except OSError:
                    pass

        true_sym = entry.get("true_gene_symbol")
        rank = find_rank(ranked_genes, true_sym)

        append_ranking({
            "patient_id": pid,
            "disease_id": entry.get("disease_id"),
            "broad_category": entry.get("broad_category"),
            "true_gene_symbol": true_sym,
            "true_gene_ensg": entry.get("true_gene_ensg"),
            "exomiser_rank": rank,
            "exomiser_n_genes": len(ranked_genes),
            "exomiser_top_genes": ranked_genes[:TOP_N],
        })
        result["status"] = "ok"

    except subprocess.TimeoutExpired:
        result["status"] = "error"
        result["error"] = "Timeout (>600s)"
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)

    return result


# ── validation ────────────────────────────────────────────────────────────────

def validate_setup() -> None:
    errors = []
    if not EXOMISER_JAR.exists():
        errors.append(
            f"Exomiser JAR not found: {EXOMISER_JAR}\n"
            f"  → Run: bash scripts/00_download_setup.sh"
        )
    if not (EXOMISER_DATA_DIR / f"{DATA_RELEASE}_hg19").exists():
        errors.append(
            f"Exomiser hg19 data not found: {EXOMISER_DATA_DIR}/{DATA_RELEASE}_hg19\n"
            f"  → Run: bash scripts/00_download_setup.sh"
        )
    if not (EXOMISER_DATA_DIR / f"{DATA_RELEASE}_phenotype").exists():
        errors.append(
            f"Exomiser phenotype data not found: {EXOMISER_DATA_DIR}/{DATA_RELEASE}_phenotype\n"
            f"  → Run: bash scripts/00_download_setup.sh"
        )
    if not MANIFEST_PATH.exists():
        errors.append(
            f"Manifest not found: {MANIFEST_PATH}\n"
            f"  → Run: python scripts/01_prepare_inputs.py"
        )
    if errors:
        print("Setup errors:")
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(1)


# ── main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    global RANKINGS_PATH, ERRORS_PATH, OUTPUTS_DIR, APP_PROPS_PATH

    parser = argparse.ArgumentParser(description="Run Exomiser 14.x on all benchmark cases")
    parser.add_argument(
        "--workers", type=int, default=2,
        help="Parallel worker threads (default: 2; each JVM uses ~6 GB RAM)",
    )
    parser.add_argument(
        "--limit", type=int, default=None,
        help="Only process first N cases (for testing)",
    )
    parser.add_argument(
        "--test", action="store_true",
        help="Run first 5 cases only (alias for --limit 5)",
    )
    parser.add_argument(
        "--dry-run", action="store_true",
        help="Print commands without executing",
    )
    parser.add_argument(
        "--human-only", action="store_true",
        help=(
            "Restrict hiPhivePrioritiser to human phenotype associations only "
            "(ablation run; results go to exomiser_human_rankings.jsonl). "
            "Reference: Genome Medicine 2025 DOI:10.1186/s13073-025-01546-1"
        ),
    )
    args = parser.parse_args()

    if args.test:
        args.limit = 5

    # Redirect output paths for human-only ablation run
    if args.human_only:
        RANKINGS_PATH = REPO_ROOT / "results" / "exomiser_human_rankings.jsonl"
        ERRORS_PATH = REPO_ROOT / "results" / "exomiser_human_errors.jsonl"
        OUTPUTS_DIR = REPO_ROOT / "outputs" / "exomiser_human"
        APP_PROPS_PATH = OUTPUTS_DIR / "application.properties"
        print("[MODE] Human-only hiPhivePrioritiser (ablation run)")
        print(f"       Rankings → {RANKINGS_PATH}")

    validate_setup()

    # Load manifest
    entries = []
    with open(MANIFEST_PATH) as f:
        for line in f:
            entries.append(json.loads(line))

    if args.limit:
        entries = entries[: args.limit]
        print(f"[LIMIT] Processing first {args.limit} cases")

    # Build resume set
    done_pids: set = set()
    if RANKINGS_PATH.exists():
        with open(RANKINGS_PATH) as f:
            for line in f:
                try:
                    done_pids.add(json.loads(line)["patient_id"])
                except (json.JSONDecodeError, KeyError):
                    pass
        if done_pids:
            print(f"[RESUME] {len(done_pids)} cases already in {RANKINGS_PATH.name}")

    total = len(entries)
    to_run = sum(1 for e in entries if e["patient_id"] not in done_pids)
    print(f"Total cases   : {total}")
    print(f"Already done  : {total - to_run}")
    print(f"To run        : {to_run}")
    print(f"Memory        : {args.workers} workers × ~6 GB JVM = ~{args.workers * 6} GB RAM")

    if to_run == 0:
        print("All cases already complete. Run 04_compare_metrics.py next.")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    RANKINGS_PATH.parent.mkdir(parents=True, exist_ok=True)

    # Write shared application.properties once
    write_application_properties(APP_PROPS_PATH)
    print(f"Application properties: {APP_PROPS_PATH}")

    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
        print("[INFO] Install tqdm for a progress bar: pip install tqdm")

    stats: dict = {}
    error_log = []

    print(f"\nRunning Exomiser with {args.workers} workers...")
    print(f"Rankings appended live to: {RANKINGS_PATH}")
    print("(Per-case TSV files deleted immediately after parsing)")

    if args.dry_run:
        for entry in entries[:3]:
            r = run_case(entry, done_pids, dry_run=True, human_only=args.human_only)
            print("\nCommand:")
            print(" ", r.get("cmd", ""))
        print(f"\n[DRY RUN] Would process {to_run} cases.")
        return

    futures_map = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for entry in entries:
            future = executor.submit(run_case, entry, done_pids, False, args.human_only)
            futures_map[future] = entry["patient_id"]

        iterator = as_completed(futures_map)
        if use_tqdm:
            iterator = tqdm(iterator, total=total, desc="Exomiser", unit="case")

        for future in iterator:
            result = future.result()
            status = result["status"]
            stats[status] = stats.get(status, 0) + 1

            if status == "error":
                error_log.append(result)
                print(
                    f"\n  [ERROR] patient {result['patient_id']}: "
                    f"{result['error'][:120]}",
                    file=sys.stderr,
                )

    if error_log:
        with open(ERRORS_PATH, "w") as ef:
            for e in error_log:
                ef.write(json.dumps(e) + "\n")
        print(f"\n{len(error_log)} errors logged to: {ERRORS_PATH}")

    print(f"\n=== Run complete ===")
    print(f"  OK      : {stats.get('ok', 0)}")
    print(f"  Skipped : {stats.get('skipped', 0)}")
    print(f"  Errors  : {stats.get('error', 0)}")
    print(f"\nRankings: {RANKINGS_PATH}")
    print("Next step: python scripts/04_compare_metrics.py")


if __name__ == "__main__":
    main()
