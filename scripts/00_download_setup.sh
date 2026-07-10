#!/usr/bin/env bash
# 00_download_setup.sh
# Downloads and sets up LIRICAL v2.2.1 and/or Exomiser 14.0.0 with required data.
#
# Usage:
#   bash scripts/00_download_setup.sh                    # set up everything (default: all)
#   bash scripts/00_download_setup.sh --setup lirical    # LIRICAL only
#   bash scripts/00_download_setup.sh --setup exomiser   # Exomiser only
#   bash scripts/00_download_setup.sh --setup all        # both tools + shared data
#
# KNOWN ISSUE (LIRICAL): `lirical download` fails with HTTP 403 on UCSC Jannovar files
# (storage.googleapis.com/ielis bucket is no longer publicly accessible).
# Workaround:
#   1. Run `lirical download` with || true to capture non-UCSC files (hp.json etc.)
#   2. Download Jannovar .ser files from Zenodo (non-UCSC) manually
#   3. Create dummy UCSC .ser placeholders (never read since --transcript-db REFSEQ)
#
# Requirements: Java 17+, wget or curl, unzip
# Disk:  LIRICAL only  ~1 GB  (+4 GB hg19 data shared with Exomiser)
#        Exomiser only ~500 MB CLI + ~4 GB hg19 + ~200 MB phenotype
#        all           ~6 GB total

#set -euo pipefail

# ── argument parsing ──────────────────────────────────────────────────────────
SETUP_TARGET="all"

usage() {
    echo "Usage: bash scripts/00_download_setup.sh [--setup lirical|exomiser|all]"
    echo "  --setup lirical    Download LIRICAL JAR + data + shared hg19 data"
    echo "  --setup exomiser   Download Exomiser CLI + shared hg19 data + phenotype data"
    echo "  --setup all        Download everything (default)"
    exit 0
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --setup)
            SETUP_TARGET="$(echo "$2" | tr '[:upper:]' '[:lower:]')"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "ERROR: Unknown argument: $1"
            usage
            ;;
    esac
done

case "$SETUP_TARGET" in
    lirical|exomiser|all) ;;
    *)
        echo "ERROR: --setup must be 'lirical', 'exomiser', or 'all'  (got: '$SETUP_TARGET')"
        exit 1
        ;;
esac

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
TOOLS_DIR="$ROOT_DIR/tools"
EXOMISER_DATA_DIR="$TOOLS_DIR/exomiser_data"   # shared hg19 + phenotype data

# ── version / URL constants ───────────────────────────────────────────────────
# LIRICAL
LIRICAL_VERSION="2.2.1"
LIRICAL_DIST_ZIP="lirical-cli-${LIRICAL_VERSION}-distribution.zip"
LIRICAL_DIST_URL="https://github.com/TheJacksonLaboratory/LIRICAL/releases/download/v${LIRICAL_VERSION}/${LIRICAL_DIST_ZIP}"
LIRICAL_JAR_PATH="$TOOLS_DIR/lirical-cli-${LIRICAL_VERSION}/lirical-cli-${LIRICAL_VERSION}.jar"
LIRICAL_DATA_DIR="$TOOLS_DIR/lirical-cli-${LIRICAL_VERSION}/data"

# Zenodo record 5410367 — Jannovar transcript .ser files (non-UCSC)
ZENODO_BASE="https://zenodo.org/record/5410367/files"

# Exomiser  (docs: https://exomiser.readthedocs.io/en/stable/installation.html)
EXOMISER_VERSION="14.0.0"
# Distribution ZIP is hosted at data.monarchinitiative.org, NOT GitHub.
# The filename includes the '-distribution' suffix.
EXOMISER_DIST_ZIP="exomiser-cli-${EXOMISER_VERSION}-distribution.zip"
EXOMISER_DIST_URL="https://data.monarchinitiative.org/exomiser/latest/${EXOMISER_DIST_ZIP}"
EXOMISER_EXTRACT_DIR="$TOOLS_DIR/exomiser-cli-${EXOMISER_VERSION}"
EXOMISER_JAR="$EXOMISER_EXTRACT_DIR/exomiser-cli-${EXOMISER_VERSION}.jar"

# Shared data bundles (used by both LIRICAL -ed19 and Exomiser)
DATA_RELEASE="2406"
EXOMISER_DATA_BASE="https://data.monarchinitiative.org/exomiser/latest"

# ── helpers ───────────────────────────────────────────────────────────────────
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
warn() { echo "[$(date '+%H:%M:%S')] WARN: $*"; }
section() { echo ""; log "══ $* ══"; }

check_java() {
    if ! command -v java &>/dev/null; then
        echo "ERROR: Java not found. Install Java 17+ and retry."
        exit 1
    fi
    java_version=$(java -version 2>&1 | awk -F '"' '/version/ {print $2}' | cut -d. -f1)
    if [[ "$java_version" -lt 17 ]]; then
        echo "ERROR: Java 17+ required (found Java $java_version)."
        exit 1
    fi
    log "Java $java_version found."
}

download_file() {
    local url="$1"
    local dest="$2"
    if [[ -f "$dest" ]] && [[ -s "$dest" ]]; then
        log "Already exists, skipping: $(basename "$dest")"
        return
    fi
    log "Downloading $(basename "$dest") ..."
    if command -v wget &>/dev/null; then
        wget -q --show-progress -O "$dest" "$url"
    elif command -v curl &>/dev/null; then
        curl -L --progress-bar -o "$dest" "$url"
    else
        echo "ERROR: Neither wget nor curl found. Install one and retry."
        exit 1
    fi
    log "Done: $(basename "$dest")"
}

unzip_if_needed() {
    local zip_file="$1"
    local dest_dir="$2"
    if [[ -d "$dest_dir" ]]; then
        log "Already extracted, skipping: $(basename "$zip_file")"
        return
    fi
    log "Extracting $(basename "$zip_file") ..."
    unzip -q "$zip_file" -d "$(dirname "$dest_dir")"
    log "Extracted to $dest_dir"
}

# ── pre-flight ────────────────────────────────────────────────────────────────
log "Setup target : $SETUP_TARGET"
check_java
mkdir -p "$TOOLS_DIR" "$EXOMISER_DATA_DIR"

# ─────────────────────────────────────────────────────────────────────────────
# LIRICAL SETUP (steps 1-3)
# ─────────────────────────────────────────────────────────────────────────────
if [[ "$SETUP_TARGET" == "lirical" || "$SETUP_TARGET" == "all" ]]; then

    # ── 1. LIRICAL distribution ZIP ──────────────────────────────────────────
    section "1. LIRICAL CLI"
    LIRICAL_ZIP_DEST="$TOOLS_DIR/$LIRICAL_DIST_ZIP"
    download_file "$LIRICAL_DIST_URL" "$LIRICAL_ZIP_DEST"
    unzip_if_needed "$LIRICAL_ZIP_DEST" "$TOOLS_DIR/lirical-cli-${LIRICAL_VERSION}"

    if [[ ! -f "$LIRICAL_JAR_PATH" ]]; then
        echo "ERROR: LIRICAL JAR not found at $LIRICAL_JAR_PATH"
        echo "  Contents of $TOOLS_DIR/lirical-cli-${LIRICAL_VERSION}/:"
        ls "$TOOLS_DIR/lirical-cli-${LIRICAL_VERSION}/" 2>/dev/null || true
        exit 1
    fi
    log "LIRICAL JAR confirmed: $LIRICAL_JAR_PATH"

    # ── 2. LIRICAL phenotype data (hp.json, phenotype.hpoa, hgnc, …) ─────────
    section "2. LIRICAL phenotype data"
    # `lirical download` FAILS on UCSC .ser files (GCS bucket 403) — expected.
    # We use || true so the script continues; the important hp.json files come first.
    mkdir -p "$LIRICAL_DATA_DIR"
    if [[ -f "$LIRICAL_DATA_DIR/hp.json" ]]; then
        log "LIRICAL phenotype data (hp.json) already present, skipping."
    else
        log "Running 'lirical download' (UCSC transcript 403 is expected)..."
        java -jar "$LIRICAL_JAR_PATH" download -d "$LIRICAL_DATA_DIR" 2>&1 || true
        log "lirical download finished (UCSC 403 error handled below)."
    fi

    # ── 3. Jannovar transcript .ser files from Zenodo ─────────────────────────
    section "3. Jannovar transcript .ser files"
    # UCSC files are inaccessible (GCS 403); download all others from Zenodo.
    # Dummy placeholders created for UCSC (never read: --transcript-db REFSEQ).
    log "Downloading from Zenodo record 5410367..."
    download_file "${ZENODO_BASE}/refseq_105_hg19.ser"          "$LIRICAL_DATA_DIR/hg19_refseq.ser"
    download_file "${ZENODO_BASE}/ensembl_87_hg19.ser"          "$LIRICAL_DATA_DIR/hg19_ensembl.ser"
    download_file "${ZENODO_BASE}/refseq_curated_105_hg19.ser"  "$LIRICAL_DATA_DIR/hg19_refseq_curated.ser"
    download_file "${ZENODO_BASE}/refseq_109_hg38.ser"          "$LIRICAL_DATA_DIR/hg38_refseq.ser"
    download_file "${ZENODO_BASE}/ensembl_91_hg38.ser"          "$LIRICAL_DATA_DIR/hg38_ensembl.ser"
    download_file "${ZENODO_BASE}/refseq_curated_109_hg38.ser"  "$LIRICAL_DATA_DIR/hg38_refseq_curated.ser"

    for ucsc_file in hg19_ucsc.ser hg38_ucsc.ser; do
        ucsc_path="$LIRICAL_DATA_DIR/$ucsc_file"
        if [[ ! -f "$ucsc_path" ]]; then
            warn "Creating placeholder for $ucsc_file (GCS bucket 403-inaccessible)"
            warn "NOT used — LIRICAL runs with --transcript-db REFSEQ"
            touch "$ucsc_path"
        fi
    done

    if [[ ! -f "$LIRICAL_DATA_DIR/hp.json" ]]; then
        echo ""
        echo "ERROR: hp.json is missing from $LIRICAL_DATA_DIR"
        echo "  Try running manually:"
        echo "    java -jar $LIRICAL_JAR_PATH download -d $LIRICAL_DATA_DIR"
        exit 1
    fi

    log "LIRICAL data directory:"
    ls -lh "$LIRICAL_DATA_DIR" | awk '{print "  "$0}'

fi  # end LIRICAL setup

# ─────────────────────────────────────────────────────────────────────────────
# SHARED DATA: hg19 variant data (needed by LIRICAL -ed19 AND Exomiser)
# ─────────────────────────────────────────────────────────────────────────────
# Always download for both tools; Exomiser-only setups also need it.
section "4. Shared hg19 variant data bundle (~4 GB)"
HG19_ZIP="$EXOMISER_DATA_DIR/${DATA_RELEASE}_hg19.zip"
download_file "${EXOMISER_DATA_BASE}/${DATA_RELEASE}_hg19.zip" "$HG19_ZIP"
unzip_if_needed "$HG19_ZIP" "$EXOMISER_DATA_DIR/${DATA_RELEASE}_hg19"

# ─────────────────────────────────────────────────────────────────────────────
# EXOMISER SETUP (steps 5-6)
# ─────────────────────────────────────────────────────────────────────────────
if [[ "$SETUP_TARGET" == "exomiser" || "$SETUP_TARGET" == "all" ]]; then

    # ── 5. Exomiser phenotype data bundle ─────────────────────────────────────
    section "5. Exomiser phenotype data bundle (~200 MB)"
    PHENO_ZIP="$EXOMISER_DATA_DIR/${DATA_RELEASE}_phenotype.zip"
    download_file "${EXOMISER_DATA_BASE}/${DATA_RELEASE}_phenotype.zip" "$PHENO_ZIP"
    unzip_if_needed "$PHENO_ZIP" "$EXOMISER_DATA_DIR/${DATA_RELEASE}_phenotype"

    # ── 6. Exomiser CLI distribution ──────────────────────────────────────────
    # Source: https://exomiser.readthedocs.io/en/stable/installation.html
    # ZIP URL  : https://data.monarchinitiative.org/exomiser/latest/exomiser-cli-14.0.0-distribution.zip
    # After unzip: exomiser-cli-14.0.0/exomiser-cli-14.0.0.jar
    #              exomiser-cli-14.0.0/application.properties   (template)
    #              exomiser-cli-14.0.0/examples/
    section "6. Exomiser CLI $EXOMISER_VERSION"
    EXOMISER_ZIP_DEST="$TOOLS_DIR/$EXOMISER_DIST_ZIP"
    download_file "$EXOMISER_DIST_URL" "$EXOMISER_ZIP_DEST"
    unzip_if_needed "$EXOMISER_ZIP_DEST" "$EXOMISER_EXTRACT_DIR"

    if [[ ! -f "$EXOMISER_JAR" ]]; then
        echo "ERROR: Exomiser JAR not found at $EXOMISER_JAR"
        echo "  Contents of $EXOMISER_EXTRACT_DIR/:"
        ls "$EXOMISER_EXTRACT_DIR/" 2>/dev/null || true
        echo ""
        echo "  The ZIP may extract to a differently-named directory."
        echo "  Check the above listing and update EXOMISER_VERSION if needed."
        exit 1
    fi
    log "Exomiser JAR confirmed: $EXOMISER_JAR"

    # ── configure application.properties (in distribution directory) ──────────
    # This lets users run Exomiser manually from the tools dir.
    # The benchmark runner (06_run_exomiser.py) writes its own copy at runtime.
    APP_PROPS="$EXOMISER_EXTRACT_DIR/application.properties"
    if [[ -f "$APP_PROPS" ]]; then
        # Patch the data-directory into the bundled application.properties
        if grep -q "exomiser.data-directory" "$APP_PROPS"; then
            sed -i.bak "s|.*exomiser\.data-directory.*|exomiser.data-directory=${EXOMISER_DATA_DIR}|" "$APP_PROPS"
        else
            echo "exomiser.data-directory=${EXOMISER_DATA_DIR}" >> "$APP_PROPS"
        fi
        # Set data version (add if missing)
        grep -q "exomiser.hg19.data-version" "$APP_PROPS" || \
            echo "exomiser.hg19.data-version=${DATA_RELEASE}" >> "$APP_PROPS"
        grep -q "exomiser.phenotype.data-version" "$APP_PROPS" || \
            echo "exomiser.phenotype.data-version=${DATA_RELEASE}" >> "$APP_PROPS"
        log "Patched application.properties: $APP_PROPS"
    else
        # Create minimal application.properties if none bundled
        cat > "$APP_PROPS" <<EOF
exomiser.data-directory=${EXOMISER_DATA_DIR}
exomiser.hg19.data-version=${DATA_RELEASE}
exomiser.phenotype.data-version=${DATA_RELEASE}
spring.jpa.open-in-view=false
EOF
        log "Created application.properties: $APP_PROPS"
    fi

    log "Exomiser data directory:"
    ls -lh "$EXOMISER_DATA_DIR" | awk '{print "  "$0}'

fi  # end Exomiser setup

# ── summary ───────────────────────────────────────────────────────────────────
log ""
log "=== Setup complete ($SETUP_TARGET) ==="

if [[ "$SETUP_TARGET" == "lirical" || "$SETUP_TARGET" == "all" ]]; then
    log "  LIRICAL JAR    : $LIRICAL_JAR_PATH"
    log "  LIRICAL data   : $LIRICAL_DATA_DIR"
fi

log "  hg19 data      : $EXOMISER_DATA_DIR/${DATA_RELEASE}_hg19"

if [[ "$SETUP_TARGET" == "exomiser" || "$SETUP_TARGET" == "all" ]]; then
    log "  Phenotype data : $EXOMISER_DATA_DIR/${DATA_RELEASE}_phenotype"
    log "  Exomiser JAR   : $EXOMISER_JAR"
    log "  app.properties : $EXOMISER_EXTRACT_DIR/application.properties"
fi

log ""
if [[ "$SETUP_TARGET" == "lirical" || "$SETUP_TARGET" == "all" ]]; then
    log "NOTE: UCSC .ser files are placeholders — LIRICAL runs with --transcript-db REFSEQ."
fi
log ""
log "Next steps:"
[[ "$SETUP_TARGET" == "lirical" || "$SETUP_TARGET" == "all" ]] && \
    log "  python scripts/01_prepare_inputs.py"
[[ "$SETUP_TARGET" == "lirical" || "$SETUP_TARGET" == "all" ]] && \
    log "  python scripts/02_run_lirical.py"
[[ "$SETUP_TARGET" == "exomiser" || "$SETUP_TARGET" == "all" ]] && \
    log "  python scripts/06_run_exomiser.py"
