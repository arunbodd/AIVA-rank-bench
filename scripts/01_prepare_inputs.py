#!/usr/bin/env python3
"""
01_prepare_inputs.py
Build a manifest (results/manifest.jsonl) from the AIVA benchmark input JSONL.
Each manifest entry contains all data needed to run one LIRICAL case.

Key actions:
  - Extract HPO IDs (observed + negated) from input records
  - Resolve true_genes ENSG IDs → gene symbols using:
      1. Variants field in JSONL
      2. VCF INFO fields (handles mitochondrial gene ID swap)
      3. mygene.info API (fallback for the remaining ~60 unresolved)
  - Verify VCF file exists per case
  - Write results/manifest.jsonl

Usage:
    python scripts/01_prepare_inputs.py [--dry-run]
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

# ── paths ─────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent
AIVA_ROOT = Path("/data/AIVA_Benchmarking/AIVA-var-bench/Simulated_Datasets")
INPUT_JSONL = AIVA_ROOT / "data" / "final_extracted_broad_cat_v2_inputs.jsonl"
VCF_DIR = AIVA_ROOT / "data" / "generated_vcfs" / "simulated"
RESULTS_DIR = REPO_ROOT / "results"
INPUTS_DIR = REPO_ROOT / "inputs"
MANIFEST_PATH = INPUTS_DIR / "manifest.jsonl"


# ── ENSG → symbol resolution ──────────────────────────────────────────────────

def build_symbol_map_from_jsonl(input_jsonl: Path) -> dict:
    """Pass 1: collect ENSG→symbol from variants fields in the JSONL."""
    mapping = {}
    with open(input_jsonl) as f:
        for line in f:
            rec = json.loads(line)
            for v in rec.get("variants", []):
                gid = v.get("gene_id", "")
                sym = v.get("gene_symbol", "")
                if gid and sym:
                    # Normal: gene_id is ENSG..., gene_symbol is symbol
                    if gid.startswith("ENSG") or gid.startswith("ENS"):
                        mapping[gid] = sym
                    # Mitochondrial swap: gene_id is symbol, gene_symbol is ENSG
                    elif sym.startswith("ENSG") or sym.startswith("ENS"):
                        mapping[sym] = gid
    return mapping


def parse_vcf_info(info_str: str) -> dict:
    """Parse VCF INFO field into key→value dict."""
    result = {}
    for part in info_str.split(";"):
        if "=" in part:
            k, _, v = part.partition("=")
            result[k.strip()] = v.strip()
    return result


def build_symbol_map_from_vcf(vcf_path: Path) -> dict:
    """
    Parse a VCF file and extract ENSG→symbol pairs from INFO fields.
    Handles the mitochondrial swap where ensembl_id=SYMBOL and name=ENSG.
    """
    mapping = {}
    try:
        with open(vcf_path) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 8:
                    continue
                info = parse_vcf_info(parts[7])
                eid = info.get("ensembl_id", "")
                name = info.get("name", "")
                if not eid or not name:
                    continue
                # Normal case: ensembl_id=ENSG..., name=SYMBOL
                if eid.startswith("ENSG") or eid.startswith("ENS"):
                    mapping[eid] = name
                # Mitochondrial swap: ensembl_id=SYMBOL (e.g. MT-TE), name=ENSG...
                elif name.startswith("ENSG") or name.startswith("ENS"):
                    mapping[name] = eid
    except (OSError, IOError):
        pass
    return mapping


def resolve_via_mygene(ensg_ids: list) -> dict:
    """
    Use mygene.info REST API to resolve remaining unresolved ENSG/gene IDs.
    Returns a dict of id → symbol.
    """
    if not ensg_ids:
        return {}
    try:
        import urllib.request
        import urllib.parse
        mapping = {}
        # Batch in chunks of 100
        chunk_size = 100
        for i in range(0, len(ensg_ids), chunk_size):
            chunk = ensg_ids[i : i + chunk_size]
            query = ",".join(chunk)
            url = "https://mygene.info/v3/gene"
            data = urllib.parse.urlencode({
                "ids": query,
                "fields": "symbol,ensembl.gene",
                "species": "human",
            }).encode()
            req = urllib.request.Request(url, data=data, method="POST")
            req.add_header("Content-Type", "application/x-www-form-urlencoded")
            with urllib.request.urlopen(req, timeout=30) as resp:
                results = json.loads(resp.read())
            for entry in results:
                if "symbol" in entry:
                    # Map by query id (could be ENSG or numeric)
                    qid = entry.get("query", "")
                    mapping[qid] = entry["symbol"]
        return mapping
    except Exception as e:
        print(f"  [WARN] mygene.info lookup failed: {e}", file=sys.stderr)
        return {}


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Prepare LIRICAL input manifest")
    parser.add_argument("--dry-run", action="store_true",
                        help="Parse and report stats without writing output")
    args = parser.parse_args()

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    print("[1/4] Building ENSG→symbol map from input JSONL variants...")
    symbol_map = build_symbol_map_from_jsonl(INPUT_JSONL)
    print(f"      {len(symbol_map)} mappings collected from variants field")

    # Load all records first
    records = []
    with open(INPUT_JSONL) as f:
        for line in f:
            records.append(json.loads(line))
    print(f"      Loaded {len(records)} records")

    # Find unresolved true genes
    unresolved = set()
    for rec in records:
        for tg in rec.get("true_genes", []):
            if tg not in symbol_map:
                unresolved.add(tg)

    print(f"\n[2/4] Resolving {len(unresolved)} unresolved true genes via VCF files...")
    vcf_resolved = 0
    for rec in records:
        pid = rec["id"]
        vcf_path = VCF_DIR / f"simulated_{pid}_hg19.vcf"
        if not vcf_path.exists():
            continue
        # Only parse VCF if it might help resolve an unresolved gene
        needs_vcf = any(tg in unresolved for tg in rec.get("true_genes", []))
        if needs_vcf:
            vcf_map = build_symbol_map_from_vcf(vcf_path)
            for ensg, sym in vcf_map.items():
                if ensg in unresolved:
                    symbol_map[ensg] = sym
                    unresolved.discard(ensg)
                    vcf_resolved += 1
                else:
                    symbol_map[ensg] = sym  # add anyway for completeness

    print(f"      {vcf_resolved} additional mappings from VCF files")
    print(f"      Still unresolved: {len(unresolved)}")

    if unresolved:
        print(f"\n[3/4] Fetching {len(unresolved)} remaining symbols from mygene.info...")
        mg_map = resolve_via_mygene(sorted(unresolved))
        symbol_map.update(mg_map)
        still_missing = unresolved - set(mg_map.keys())
        print(f"      Resolved via API: {len(mg_map)}")
        if still_missing:
            print(f"      [WARN] Could not resolve {len(still_missing)} IDs: {sorted(still_missing)[:10]}")
    else:
        print("\n[3/4] No mygene.info lookups needed.")

    # ── build manifest ────────────────────────────────────────────────────────
    print("\n[4/4] Building manifest...")
    manifest = []
    stats = {"total": 0, "vcf_missing": 0, "gene_unresolved": 0, "ok": 0}

    for rec in records:
        pid = rec["id"]
        stats["total"] += 1

        # HPO terms
        observed = list(rec.get("positive_phenotypes", {}).keys())
        negated = list(rec.get("negative_phenotypes", {}).keys())

        # True gene symbol
        true_ensg = rec.get("true_genes", [None])[0]
        true_symbol = symbol_map.get(true_ensg, None) if true_ensg else None
        if true_symbol is None:
            stats["gene_unresolved"] += 1

        # VCF path
        vcf_path = VCF_DIR / f"simulated_{pid}_hg19.vcf"
        vcf_exists = vcf_path.exists()
        if not vcf_exists:
            stats["vcf_missing"] += 1

        entry = {
            "patient_id": pid,
            "disease_id": rec.get("disease_id"),
            "broad_category": rec.get("broad_category"),
            "category": rec.get("category"),
            "observed_hpo": observed,
            "negated_hpo": negated,
            "vcf_path": str(vcf_path),
            "vcf_exists": vcf_exists,
            "true_gene_ensg": true_ensg,
            "true_gene_symbol": true_symbol,
        }
        manifest.append(entry)
        if vcf_exists:
            stats["ok"] += 1

    # Stats report
    print(f"\n  Total cases     : {stats['total']}")
    print(f"  VCF missing     : {stats['vcf_missing']}")
    print(f"  Gene unresolved : {stats['gene_unresolved']}")
    print(f"  Ready to run    : {stats['ok']}")

    if args.dry_run:
        print("\n[DRY RUN] Not writing manifest.")
        return

    with open(MANIFEST_PATH, "w") as out:
        for entry in manifest:
            out.write(json.dumps(entry) + "\n")

    print(f"\nManifest written: {MANIFEST_PATH} ({len(manifest)} records)")
    print("Next step: python scripts/02_run_lirical.py")


if __name__ == "__main__":
    main()
