#!/usr/bin/env bash
# Download every raw file listed in metadata/sources.tsv, check it, and make it read-only.
#
# - Run from the project folder:  bash scripts/download.sh
# - Files that are already present are not downloaded again.
# - If the provider publishes an md5 (column `md5`), the download is checked against it.
# - curl, md5sum and sha256sum come from the Pixi environment (curl and coreutils in
#   pixi.toml), so the script works the same on macOS, Linux and in the container.
# - If metadata/checksums.sha256 exists, every file is checked against it at the end,
#   so you know you have exactly the same data as the author.
set -euo pipefail

SOURCES=metadata/sources.tsv
CHECKSUMS=metadata/checksums.sha256

# skip the header line, then read the five tab-separated columns
tail -n +2 "$SOURCES" | while IFS=$'\t' read -r name url md5 unpack dest; do
    mkdir -p "$dest"
    file="$dest/$(basename "$url")"
    if [[ -f "$file" ]]; then
        echo "already have  $file"
        continue
    fi
    echo "downloading   $name -> $file"
    curl --fail --location --silent --show-error --retry 3 -o "$file.part" "$url"
    if [[ "$md5" != "NA" ]]; then
        got=$(md5sum "$file.part" | cut -d' ' -f1)
        if [[ "$got" != "$md5" ]]; then
            echo "md5 mismatch for $file: expected $md5, got $got" >&2
            exit 1
        fi
    fi
    mv "$file.part" "$file"
    chmod a-w "$file"               # raw data is read-only
done

if [[ -f "$CHECKSUMS" ]]; then
    echo "checking files against $CHECKSUMS"
    sha256sum -c "$CHECKSUMS"
else
    echo "no $CHECKSUMS yet: create it once with"
    echo "  sha256sum data/raw/geo/* data/raw/go/* > $CHECKSUMS"
fi
