#!/usr/bin/env python3
"""
08_run_clinprior.py
Run ClinPrior v4 in phenotype-only mode on every case in the manifest.

ClinPrior is run without variant prioritization (proteinScore + MatrixPropagation
only) because the simulated VCFs use a simplified CSQ annotation format that is
not compatible with ClinPrior's readVCF() parser, which expects full VEP output.
The resulting gene ranking reflects pure phenotype-driven network propagation.

For each case:
  1. Extract HPO terms from manifest
  2. Run ClinPrior::proteinScore() → gene-level phenotypic scores
  3. Run ClinPrior::MatrixPropagation() → propagated functional + physical scores
  4. Rank genes by mean(score_funct, score_phys) and record true gene rank
  5. Append compact record to results/clinprior_rankings.jsonl

Features:
  - Parallel execution via N independent R worker subprocesses (--workers N)
  - Resume: skips cases already in clinprior_rankings.jsonl
  - Progress bar (tqdm if installed)
  - Errors logged to results/clinprior_errors.jsonl
  - Dry-run mode (--dry-run)

Usage:
    python scripts/08_run_clinprior.py [--workers 4] [--limit N] [--dry-run]
    python scripts/08_run_clinprior.py --test    # first 5 cases only
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT   = Path(__file__).resolve().parent.parent
R_WORKER    = REPO_ROOT / "scripts" / "_clinprior_worker.R"
MANIFEST    = REPO_ROOT / "inputs" / "manifest.jsonl"
RANKINGS    = REPO_ROOT / "results" / "clinprior_rankings.jsonl"
ERRORS      = REPO_ROOT / "results" / "clinprior_errors.jsonl"

TOP_N       = 20   # max genes stored per case
DEFAULT_WORKERS = 4


# ── helpers ───────────────────────────────────────────────────────────────────

def load_done_ids(path: Path) -> set:
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            line = line.strip()
            if line:
                try:
                    done.add(json.loads(line)["patient_id"])
                except Exception:
                    pass
    return done


def split_into_batches(items, n):
    """Split items into n roughly equal batches."""
    k, r = divmod(len(items), n)
    batches, start = [], 0
    for i in range(n):
        end = start + k + (1 if i < r else 0)
        if start < end:
            batches.append(items[start:end])
        start = end
    return batches


def run_worker(batch, worker_id, results_lock, results_list, errors_list, args):
    """Write batch JSONL, call Rscript worker, collect output."""
    with tempfile.TemporaryDirectory(prefix="clinprior_") as tmpdir:
        in_path  = os.path.join(tmpdir, f"batch_{worker_id}_in.jsonl")
        out_path = os.path.join(tmpdir, f"batch_{worker_id}_out.jsonl")

        with open(in_path, "w") as fh:
            for case in batch:
                fh.write(json.dumps(case) + "\n")

        cmd = ["Rscript", "--vanilla", str(R_WORKER),
               in_path, out_path, str(TOP_N)]

        if args.dry_run:
            print(f"[worker {worker_id}] DRY-RUN: would run {len(batch)} cases")
            return

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=args.timeout_per_worker,
            )
        except subprocess.TimeoutExpired:
            msg = f"[worker {worker_id}] TIMEOUT after {args.timeout_per_worker}s"
            print(msg, file=sys.stderr)
            with results_lock:
                errors_list.append({"worker": worker_id, "error": "timeout",
                                    "pids": [c["patient_id"] for c in batch]})
            return

        # Forward stderr progress lines
        for line in proc.stderr.splitlines():
            print(f"[W{worker_id}] {line}", file=sys.stderr)

        if proc.returncode != 0:
            msg = f"[worker {worker_id}] R exited with code {proc.returncode}"
            print(msg, file=sys.stderr)
            with results_lock:
                errors_list.append({"worker": worker_id,
                                    "error": f"exit {proc.returncode}",
                                    "pids": [c["patient_id"] for c in batch]})
            return

        # Read results
        if os.path.exists(out_path):
            with open(out_path) as fh:
                lines = fh.read().splitlines()
            parsed = []
            for line in lines:
                line = line.strip()
                if not line:
                    continue
                try:
                    parsed.append(json.loads(line))
                except json.JSONDecodeError as e:
                    errors_list.append({"worker": worker_id,
                                        "error": f"json decode: {e}",
                                        "line": line[:200]})
            with results_lock:
                results_list.extend(parsed)
        else:
            print(f"[worker {worker_id}] No output file found!", file=sys.stderr)


def write_results(results, rankings_path, errors_path, errors_list):
    """Append results to rankings JSONL and errors JSONL."""
    if results:
        with open(rankings_path, "a") as fh:
            for rec in results:
                fh.write(json.dumps(rec, default=str) + "\n")
    if errors_list:
        with open(errors_path, "a") as fh:
            for err in errors_list:
                fh.write(json.dumps(err) + "\n")


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workers",  type=int, default=DEFAULT_WORKERS,
                        help=f"parallel R subprocesses (default {DEFAULT_WORKERS})")
    parser.add_argument("--limit",    type=int, default=None,
                        help="cap number of cases (debug)")
    parser.add_argument("--test",     action="store_true",
                        help="run first 5 cases only")
    parser.add_argument("--dry-run",  action="store_true",
                        help="print what would run, skip Rscript calls")
    parser.add_argument("--timeout-per-worker", type=int, default=7200,
                        help="seconds before killing a stalled worker (default 7200)")
    args = parser.parse_args()

    if args.test:
        args.limit = 5
        args.workers = 1

    # ── validate R worker script ───────────────────────────────────────────────
    if not R_WORKER.exists():
        sys.exit(f"ERROR: R worker not found at {R_WORKER}")

    # ── load manifest ──────────────────────────────────────────────────────────
    if not MANIFEST.exists():
        sys.exit(f"ERROR: manifest not found at {MANIFEST}")

    manifest = []
    for line in MANIFEST.read_text().splitlines():
        line = line.strip()
        if line:
            manifest.append(json.loads(line))
    print(f"[ClinPrior] Manifest: {len(manifest)} cases")

    # ── filter to-do cases ─────────────────────────────────────────────────────
    done_ids = load_done_ids(RANKINGS)
    if done_ids:
        print(f"[ClinPrior] Resuming: {len(done_ids)} already done, skipping")

    todo = [c for c in manifest if c["patient_id"] not in done_ids]
    if args.limit:
        todo = todo[:args.limit]

    if not todo:
        print("[ClinPrior] All cases already complete.")
        return

    print(f"[ClinPrior] Cases to run: {len(todo)}  workers: {args.workers}")

    # ── split into batches ─────────────────────────────────────────────────────
    n_workers = min(args.workers, len(todo))
    batches   = split_into_batches(todo, n_workers)
    print(f"[ClinPrior] Batch sizes: {[len(b) for b in batches]}")

    if args.dry_run:
        print("[ClinPrior] DRY-RUN mode: no Rscript calls.")
        return

    # ── run workers in parallel threads ────────────────────────────────────────
    results_lock = threading.Lock()
    results_list = []
    errors_list  = []

    threads = []
    t_start = time.time()
    for wid, batch in enumerate(batches):
        t = threading.Thread(
            target=run_worker,
            args=(batch, wid, results_lock, results_list, errors_list, args),
            daemon=True,
        )
        threads.append(t)
        t.start()

    try:
        for t in threads:
            t.join()
    except KeyboardInterrupt:
        print("\n[ClinPrior] Interrupted — saving partial results...",
              file=sys.stderr)

    elapsed = time.time() - t_start
    print(f"\n[ClinPrior] Workers finished in {elapsed/60:.1f} min")
    print(f"[ClinPrior] Results collected: {len(results_list)}")

    # ── write results ──────────────────────────────────────────────────────────
    write_results(results_list, RANKINGS, ERRORS, errors_list)
    print(f"[ClinPrior] Rankings written to {RANKINGS}")
    if errors_list:
        print(f"[ClinPrior] {len(errors_list)} errors written to {ERRORS}",
              file=sys.stderr)

    # ── summary ────────────────────────────────────────────────────────────────
    total_done = sum(1 for r in results_list if r.get("clinprior_rank") is not None)
    total_miss = sum(1 for r in results_list if r.get("clinprior_rank") is None)
    print(f"[ClinPrior] True gene ranked: {total_done}  absent: {total_miss}")


if __name__ == "__main__":
    main()
