#!/usr/bin/env python3
"""
02_run_lirical.py
Run LIRICAL v2.2.1 on every case in the manifest (results/manifest.jsonl).

For each case:
  1. Run LIRICAL with --global --output-format json
  2. Parse the output JSON in-memory immediately after completion
  3. Extract compact gene rankings (top-N by posttestProbability, deduped by gene)
  4. Append a compact record to results/lirical_rankings.jsonl (thread-safe)
  5. Delete the large JSON file  →  avoids ~32 GB / 820M lines disk bloat

Features:
  - Parallel execution (--workers N, default 4)
  - Resume: skips cases already written to lirical_rankings.jsonl
  - Progress bar (tqdm)
  - Errors logged to results/lirical_errors.jsonl
  - Dry-run mode (--dry-run)

Usage:
    python scripts/02_run_lirical.py [--workers 4] [--limit N] [--dry-run]
    python scripts/02_run_lirical.py --test    # run first 5 cases only
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
LIRICAL_VERSION = "2.2.1"
LIRICAL_JAR = REPO_ROOT / "tools" / f"lirical-cli-{LIRICAL_VERSION}" / f"lirical-cli-{LIRICAL_VERSION}.jar"
LIRICAL_DATA_DIR = REPO_ROOT / "tools" / f"lirical-cli-{LIRICAL_VERSION}" / "data"
EXOMISER_HG19_DIR = REPO_ROOT / "tools" / "exomiser_data" / "2406_hg19"
MANIFEST_PATH = REPO_ROOT / "inputs" / "manifest.jsonl"
OUTPUTS_DIR = REPO_ROOT / "outputs"
RANKINGS_PATH = REPO_ROOT / "results" / "lirical_rankings.jsonl"
ERRORS_PATH = REPO_ROOT / "results" / "lirical_errors.jsonl"

TOP_N = 20  # max genes to store per case

# ── thread-safe ranking writer ────────────────────────────────────────────────
_write_lock = threading.Lock()


def append_ranking(record: dict):
    """Thread-safe append of one compact ranking record to lirical_rankings.jsonl."""
    with _write_lock:
        with open(RANKINGS_PATH, "a") as f:
            f.write(json.dumps(record) + "\n")


# ── LIRICAL JSON parsing ──────────────────────────────────────────────────────

def parse_lirical_json(json_path: Path) -> list:
    """
    Parse a LIRICAL --global output JSON.
    Returns [{gene, prob, lr}, ...] sorted by posttestProbability DESC then
    genotypeLR.lr DESC (tiebreaker), deduplicated by gene symbol.

    Storage precision: prob stored to 10 decimal places (values can be < 5e-7
    in --global mode; round(prob, 6) would truncate ~79% of cases to 0.0).
    lr stored to 6 decimal places as secondary sort key for auditing.

    LIRICAL JSON structure:
      { "analysisResults": [
          { "posttestProbability": 0.85,
            "genotypeLR": { "geneId": {"symbol": "GENE"}, "lr": 10.0 },
            ... },
          ...
        ]
      }
    analysisResults is NOT pre-sorted — we sort here.
    """
    with open(json_path) as f:
        data = json.load(f)

    raw_results = data.get("analysisResults", [])
    gene_best: dict = {}  # gene -> (posttestProbability, genotypeLR.lr)

    for res in raw_results:
        prob = res.get("posttestProbability", 0.0)
        genotype_lr = res.get("genotypeLR")
        if genotype_lr is None:
            continue  # phenotype-only entry, no gene — skip
        lr = genotype_lr.get("lr") or 0.0
        gene_id_obj = genotype_lr.get("geneId", {})
        symbol = (gene_id_obj.get("symbol") or
                  gene_id_obj.get("name") or
                  gene_id_obj.get("value") or "")
        if not symbol:
            continue
        if symbol not in gene_best or (prob, lr) > gene_best[symbol]:
            gene_best[symbol] = (prob, lr)

    # Sort by posttestProbability DESC, then genotypeLR.lr DESC as tiebreaker
    # (when all probs are near-zero, lr provides a meaningful secondary ranking)
    ranked = sorted(gene_best.items(), key=lambda x: (-x[1][0], -x[1][1]))
    return [
        {"gene": sym, "prob": round(prob, 10), "lr": round(lr, 6)}
        for sym, (prob, lr) in ranked
    ]


def find_rank(ranked_genes: list, true_symbol: str):
    """Return 1-based rank of true_symbol in ranked_genes, or None if absent."""
    if not true_symbol:
        return None
    true_upper = true_symbol.upper()
    for i, entry in enumerate(ranked_genes):
        if entry["gene"] == true_symbol or entry["gene"].upper() == true_upper:
            return i + 1
    return None


# ── case runner ───────────────────────────────────────────────────────────────

def run_case(entry: dict, done_pids: set, dry_run: bool = False) -> dict:
    """
    Run LIRICAL for one case, parse inline, write compact ranking, delete large JSON.
    Returns a status dict for progress tracking.
    """
    pid = entry["patient_id"]
    t_start = time.time()
    result = {"patient_id": pid, "status": None, "elapsed_s": None, "error": None}

    # Resume: already written to rankings
    if pid in done_pids:
        result["status"] = "skipped"
        return result

    if not entry.get("vcf_exists", False):
        result["status"] = "error"
        result["error"] = f"VCF not found: {entry['vcf_path']}"
        return result

    observed = entry.get("observed_hpo", [])
    negated = entry.get("negated_hpo", [])

    if not observed:
        result["status"] = "error"
        result["error"] = "No observed HPO terms"
        return result

    out_dir = OUTPUTS_DIR / str(pid)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_json = out_dir / f"{pid}.json"

    cmd = [
        "java", "-Xmx4g", "-jar", str(LIRICAL_JAR),
        "prioritize",
        "--vcf", entry["vcf_path"],
        "--assembly", "hg19",
        "-ed19", str(EXOMISER_HG19_DIR),
        "-d", str(LIRICAL_DATA_DIR),
        "--observed-phenotypes", ",".join(observed),
        "--transcript-db", "REFSEQ",   # avoids broken UCSC .ser placeholder
        "--global",                    # output ALL diseases (required for correct ranking)
        "--output-format", "json",
        "--prefix", str(pid),
        "--output-directory", str(out_dir),
    ]
    if negated:
        cmd += ["--negated-phenotypes", ",".join(negated)]

    if dry_run:
        result["status"] = "dry_run"
        result["cmd"] = " ".join(cmd)
        return result

    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 min per case
        )
        elapsed = time.time() - t_start
        result["elapsed_s"] = round(elapsed, 1)

        if proc.returncode != 0:
            result["status"] = "error"
            result["error"] = proc.stderr[-2000:] if proc.stderr else "non-zero exit"
            return result

        if not out_json.exists():
            result["status"] = "error"
            result["error"] = f"Output JSON not created. stdout: {proc.stdout[-500:]}"
            return result

        # ── Parse inline, then delete large JSON ──────────────────────────
        try:
            ranked_genes = parse_lirical_json(out_json)
        except Exception as e:
            result["status"] = "error"
            result["error"] = f"JSON parse error: {e}"
            return result
        finally:
            # Always delete the large output JSON to reclaim disk space.
            # With --global, each file is ~23 MB; 1,396 cases = ~32 GB.
            try:
                out_json.unlink()
            except OSError:
                pass

        true_sym = entry.get("true_gene_symbol")
        rank = find_rank(ranked_genes, true_sym)

        ranking_record = {
            "patient_id": pid,
            "disease_id": entry.get("disease_id"),
            "broad_category": entry.get("broad_category"),
            "true_gene_symbol": true_sym,
            "true_gene_ensg": entry.get("true_gene_ensg"),
            "lirical_rank": rank,
            "lirical_n_genes": len(ranked_genes),
            "lirical_top_genes": ranked_genes[:TOP_N],
        }
        append_ranking(ranking_record)
        result["status"] = "ok"

    except subprocess.TimeoutExpired:
        result["status"] = "error"
        result["error"] = "Timeout (>300s)"
    except Exception as e:
        result["status"] = "error"
        result["error"] = str(e)

    return result


# ── validation ────────────────────────────────────────────────────────────────

def validate_setup():
    errors = []
    if not LIRICAL_JAR.exists():
        errors.append(f"LIRICAL JAR not found: {LIRICAL_JAR}\n  → Run: bash scripts/00_download_setup.sh")
    if not LIRICAL_DATA_DIR.exists():
        errors.append(f"LIRICAL data dir not found: {LIRICAL_DATA_DIR}\n  → Run: bash scripts/00_download_setup.sh")
    if not EXOMISER_HG19_DIR.exists():
        errors.append(f"Exomiser hg19 dir not found: {EXOMISER_HG19_DIR}\n  → Run: bash scripts/00_download_setup.sh")
    if not MANIFEST_PATH.exists():
        errors.append(f"Manifest not found: {MANIFEST_PATH}\n  → Run: python scripts/01_prepare_inputs.py")
    if errors:
        print("Setup errors:")
        for e in errors:
            print(f"  ERROR: {e}")
        sys.exit(1)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Run LIRICAL on all benchmark cases")
    parser.add_argument("--workers", type=int, default=4,
                        help="Parallel worker threads (default: 4)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Only process first N cases (for testing)")
    parser.add_argument("--test", action="store_true",
                        help="Run first 5 cases only (alias for --limit 5)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    args = parser.parse_args()

    if args.test:
        args.limit = 5

    validate_setup()

    # Load manifest
    entries = []
    with open(MANIFEST_PATH) as f:
        for line in f:
            entries.append(json.loads(line))

    if args.limit:
        entries = entries[: args.limit]
        print(f"[LIMIT] Processing first {args.limit} cases")

    # Build resume set from existing rankings file
    done_pids: set = set()
    if RANKINGS_PATH.exists():
        with open(RANKINGS_PATH) as f:
            for line in f:
                try:
                    rec = json.loads(line)
                    done_pids.add(rec["patient_id"])
                except (json.JSONDecodeError, KeyError):
                    pass
        if done_pids:
            print(f"[RESUME] {len(done_pids)} cases already in {RANKINGS_PATH.name}")

    total = len(entries)
    to_run = sum(1 for e in entries if e["patient_id"] not in done_pids)
    print(f"Total cases   : {total}")
    print(f"Already done  : {total - to_run}")
    print(f"To run        : {to_run}")

    if to_run == 0:
        print("All cases already complete. Run 04_compare_metrics.py next.")
        return

    OUTPUTS_DIR.mkdir(parents=True, exist_ok=True)
    RANKINGS_PATH.parent.mkdir(parents=True, exist_ok=True)
    ERRORS_PATH.parent.mkdir(parents=True, exist_ok=True)

    try:
        from tqdm import tqdm
        use_tqdm = True
    except ImportError:
        use_tqdm = False
        print("[INFO] Install tqdm for progress bar: pip install tqdm")

    stats: dict = {"ok": 0, "skipped": 0, "error": 0}
    error_log = []

    print(f"\nRunning LIRICAL with {args.workers} workers...")
    print(f"Rankings appended live to: {RANKINGS_PATH}")
    print("(Large per-case JSON files are deleted immediately after parsing)")

    if args.dry_run:
        for entry in entries[:3]:
            r = run_case(entry, done_pids, dry_run=True)
            print("\nCommand:")
            print(" ", r.get("cmd", ""))
        print(f"\n[DRY RUN] Would process {to_run} cases.")
        return

    futures_map = {}
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        for entry in entries:
            future = executor.submit(run_case, entry, done_pids, False)
            futures_map[future] = entry["patient_id"]

        iterator = as_completed(futures_map)
        if use_tqdm:
            iterator = tqdm(iterator, total=total, desc="LIRICAL", unit="case")

        for future in iterator:
            result = future.result()
            status = result["status"]
            stats[status] = stats.get(status, 0) + 1

            if status == "error":
                error_log.append(result)
                pid = result["patient_id"]
                print(f"\n  [ERROR] patient {pid}: {result['error'][:120]}", file=sys.stderr)

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
