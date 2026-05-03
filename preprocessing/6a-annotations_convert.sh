#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

annotations_dir="${ANNOTATIONS_DIR:-$repo_root/annotations}"
gtf_file="${GTF_FILE:-$annotations_dir/gencode.v46.basic.annotation.gtf}"
exons_bed="$annotations_dir/gencode.v46.exons.bed"

awk 'BEGIN {OFS="\t"} $3 == "exon" { start = $4 - 1; end = $5; print $1, start, end, ".", ".", $7 }' "$gtf_file" > "$exons_bed"

echo "Created exon BED: $exons_bed"
