#!/usr/bin/env bash
# Download the UCSC hg38 FASTA and install it where Arachnia looks for it.
#
# Environment overrides:
#   ARACHNIA_HG38_DIR   destination directory (default: XDG data dir/arachnia)
#   ARACHNIA_HG38_URL   alternate compressed FASTA URL
#   ARACHNIA_HG38_MD5_URL  alternate md5sum.txt URL
set -euo pipefail

URL="${ARACHNIA_HG38_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz}"
MD5_URL="${ARACHNIA_HG38_MD5_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/md5sum.txt}"
DEST_DIR="${ARACHNIA_HG38_DIR:-${XDG_DATA_HOME:-$HOME/.local/share}/arachnia}"
DEST="${DEST_DIR}/hg38.fa"
GZ="${DEST}.gz.part"
PART="${DEST}.part"

if [[ -s "$DEST" ]]; then
    echo "hg38 is already installed: $DEST"
    exit 0
fi

mkdir -p "$DEST_DIR"

cleanup() {
    rm -f "$GZ" "$PART" "$GZ.md5"
}
trap cleanup EXIT

echo "Downloading hg38 from $URL"
echo "About 940 MiB compressed; the unpacked FASTA is about 3.1 GiB."
if command -v curl >/dev/null 2>&1; then
    curl --fail --location --retry 3 --retry-delay 2 --progress-bar "$URL" --output "$GZ"
else
    wget --tries=3 --output-document="$GZ" "$URL"
fi

if command -v md5sum >/dev/null 2>&1 && command -v curl >/dev/null 2>&1; then
    if curl --fail --silent --show-error --location "$MD5_URL" --output "$GZ.md5"; then
        expected="$(awk '$2 ~ /hg38\.fa\.gz$/ { print $1; exit }' "$GZ.md5")"
        if [[ -n "$expected" ]]; then
            actual="$(md5sum "$GZ" | awk '{ print $1 }')"
            if [[ "$actual" != "$expected" ]]; then
                echo "Checksum mismatch for $GZ" >&2
                echo "  expected $expected" >&2
                echo "  actual   $actual" >&2
                exit 1
            fi
            echo "Checksum verified."
        fi
    else
        echo "Warning: could not download md5sum.txt; skipping checksum." >&2
    fi
fi

echo "Decompressing to $DEST"
gzip -dc "$GZ" >"$PART"
mv "$PART" "$DEST"
cleanup
trap - EXIT

echo "hg38 installed: $DEST"
echo "Arachnia detects this path automatically. To use it elsewhere:"
echo "  export ARACHNIA_HG38=\"$DEST\""
