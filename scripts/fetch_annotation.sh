#!/usr/bin/env bash
# Download a GENCODE human (GRCh38/hg38) annotation and install it where
# Arachnia looks for one.
#
# Environment overrides:
#   ARACHNIA_ANNOTATION_DIR     destination directory (default: XDG data dir/arachnia)
#   ARACHNIA_ANNOTATION_FORMAT  gtf (default) or gff3
#   ARACHNIA_ANNOTATION_KIND    all (default), basic, or primary_assembly
#   ARACHNIA_GENCODE_RELEASE    pin a GENCODE release, for example 49
#   ARACHNIA_GENCODE_URL        release directory that contains MD5SUMS
#   ARACHNIA_ANNOTATION_URL     direct gzip-compressed annotation URL
#   ARACHNIA_ANNOTATION_MD5_URL direct MD5SUMS URL for a custom annotation
set -euo pipefail

FORMAT="$(printf '%s' "${ARACHNIA_ANNOTATION_FORMAT:-gtf}" | tr '[:upper:]' '[:lower:]')"
KIND="$(printf '%s' "${ARACHNIA_ANNOTATION_KIND:-all}" | tr '[:upper:]' '[:lower:]')"

case "$FORMAT" in
    gtf | gff3) ;;
    *)
        echo "Error: ARACHNIA_ANNOTATION_FORMAT must be gtf or gff3 (got '$FORMAT')." >&2
        exit 1
        ;;
esac

case "$KIND" in
    all) VARIANT="annotation" ;;
    basic) VARIANT="basic.annotation" ;;
    primary_assembly) VARIANT="primary_assembly.annotation" ;;
    *)
        echo "Error: ARACHNIA_ANNOTATION_KIND must be all, basic, or primary_assembly (got '$KIND')." >&2
        exit 1
        ;;
esac

DEST_DIR="${ARACHNIA_ANNOTATION_DIR:-${XDG_DATA_HOME:-$HOME/.local/share}/arachnia}"
DEST="$DEST_DIR/gencode.annotation.$FORMAT"
GZ="${DEST}.gz.part"
PART="${DEST}.part"
SUMS="${DEST}.md5.part"
DEFAULT_BASE="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human"

if [[ -s "$DEST" ]]; then
    echo "Annotation is already installed: $DEST"
    echo "Remove it first to download a different release or format."
    exit 0
fi

mkdir -p "$DEST_DIR"

cleanup() {
    rm -f "$GZ" "$PART" "$SUMS"
}
trap cleanup EXIT

fetch_text() {
    local url="$1"
    local output="$2"
    if command -v curl >/dev/null 2>&1; then
        curl --fail --silent --show-error --location --retry 3 --retry-delay 2 \
            "$url" --output "$output"
    elif command -v wget >/dev/null 2>&1; then
        wget --quiet --tries=3 --output-document="$output" "$url"
    else
        echo "Error: curl or wget is required to download the annotation." >&2
        return 1
    fi
}

download() {
    local url="$1"
    local output="$2"
    echo "Downloading $url"
    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --retry 3 --retry-delay 2 --progress-bar \
            "$url" --output "$output"
    elif command -v wget >/dev/null 2>&1; then
        wget --tries=3 --output-document="$output" "$url"
    else
        echo "Error: curl or wget is required to download the annotation." >&2
        return 1
    fi
}

URL="${ARACHNIA_ANNOTATION_URL:-}"
MD5_URL="${ARACHNIA_ANNOTATION_MD5_URL:-}"
EXPECTED=""
FILENAME=""

if [[ -z "$URL" ]]; then
    if [[ -n "${ARACHNIA_GENCODE_RELEASE:-}" ]]; then
        RELEASE_DIR="${ARACHNIA_GENCODE_URL:-$DEFAULT_BASE/release_${ARACHNIA_GENCODE_RELEASE}}"
        echo "Looking up GENCODE release ${ARACHNIA_GENCODE_RELEASE} …"
    else
        RELEASE_DIR="${ARACHNIA_GENCODE_URL:-$DEFAULT_BASE/latest_release}"
        echo "Looking up the latest GENCODE human annotation …"
    fi

    SUMS_URL="${MD5_URL:-$RELEASE_DIR/MD5SUMS}"
    fetch_text "$SUMS_URL" "$SUMS"

    REGEX="^gencode[.]v[0-9]+[.]${VARIANT}[.]${FORMAT}[.]gz$"
    MATCH="$(awk -v re="$REGEX" '$2 ~ re { print $1 "\t" $2; exit }' "$SUMS")"
    if [[ -z "$MATCH" ]]; then
        echo "Error: $SUMS_URL does not list a GENCODE $KIND $FORMAT annotation." >&2
        exit 1
    fi
    EXPECTED="${MATCH%%$'\t'*}"
    FILENAME="${MATCH#*$'\t'}"
    URL="$RELEASE_DIR/$FILENAME"
    echo "Selected $FILENAME"
    if [[ "$KIND" == "all" ]]; then
        echo "Tip: ARACHNIA_ANNOTATION_KIND=basic installs a smaller annotation."
    fi
else
    FILENAME="$(basename "$URL")"
    if [[ -n "$MD5_URL" ]]; then
        fetch_text "$MD5_URL" "$SUMS"
        EXPECTED="$(awk -v name="$FILENAME" '$2 == name { print $1; exit }' "$SUMS")"
        if [[ -z "$EXPECTED" ]]; then
            echo "Warning: $MD5_URL does not list $FILENAME; skipping checksum." >&2
        fi
    else
        echo "Warning: no ARACHNIA_ANNOTATION_MD5_URL set; skipping checksum." >&2
    fi
fi

download "$URL" "$GZ"

if [[ -n "$EXPECTED" ]]; then
    if command -v md5sum >/dev/null 2>&1; then
        ACTUAL="$(md5sum "$GZ" | awk '{ print $1 }')"
    elif command -v md5 >/dev/null 2>&1; then
        ACTUAL="$(md5 -q "$GZ")"
    else
        echo "Warning: neither md5sum nor md5 is installed; skipping checksum." >&2
        ACTUAL=""
    fi
    if [[ -n "$ACTUAL" && "$ACTUAL" != "$EXPECTED" ]]; then
        echo "Checksum mismatch for $GZ" >&2
        echo "  expected $EXPECTED" >&2
        echo "  actual   $ACTUAL" >&2
        exit 1
    fi
    [[ -n "$ACTUAL" ]] && echo "Checksum verified."
fi

if ! command -v gzip >/dev/null 2>&1; then
    echo "Error: gzip is required to decompress the annotation." >&2
    exit 1
fi

echo "Decompressing to $DEST"
gzip -dc "$GZ" >"$PART"
if [[ ! -s "$PART" ]]; then
    echo "Error: downloaded annotation is empty." >&2
    exit 1
fi
mv "$PART" "$DEST"
cleanup
trap - EXIT

echo "Annotation installed: $DEST"
echo "Arachnia searches this path automatically. To use it elsewhere:"
echo "  export ARACHNIA_ANNOTATION=\"$DEST\""
