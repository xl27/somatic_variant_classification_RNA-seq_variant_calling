#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
annotations_dir="${ANNOTATIONS_DIR:-$repo_root/annotations}"
samples_file="${SAMPLES_FILE:-$repo_root/features/samples_id.txt}"
variants_dir="${VARIANTS_DIR:-$repo_root/features_new}"
output_dir="$annotations_dir/exons"
mkdir -p "$output_dir"

while IFS= read -r sample_id; do
  input_csv="$variants_dir/${sample_id}_features_filtered_GT_RD_MQ_AF_vcf+gnomAD.csv"
  if [[ ! -f "$input_csv" ]]; then
    echo "Missing $input_csv" >&2
    continue
  fi

  bed_file="$output_dir/${sample_id}.variants.bed"
  awk -F',' 'NR > 1 { split($1, a, "_"); print a[1] "\t" (a[2]-1) "\t" a[2] "\t" a[3] "\t" a[4] }' "$input_csv" > "$bed_file"

  bedtools intersect -a "$bed_file" -b "$annotations_dir/gencode.v46.exons.bed" -wa | awk '{print $1 "_" $3 "_" $4 "_" $5}' | sort -u > "$output_dir/${sample_id}.variants.exons.txt"
  bedtools intersect -a "$bed_file" -b "$annotations_dir/exons_edges_merged.bed" -wa | awk '{print $1 "_" $3 "_" $4 "_" $5}' | sort -u > "$output_dir/${sample_id}.variants.exon_edges.txt"
  echo "$sample_id"
done < "$samples_file"
