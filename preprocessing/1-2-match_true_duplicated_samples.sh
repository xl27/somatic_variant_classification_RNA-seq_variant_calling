#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
data_dir="${DATA_DIR:-/gpfs/commons/groups/gursoy_lab/cwalker/projects/somatic_variant_prediction/data}"
save_dir="${SAVE_DIR:-$repo_root/match_ture_tcga}"

if [[ ! -f "$save_dir/tcga_id.txt" ]]; then
  echo "Missing $save_dir/tcga_id.txt" >&2
  exit 1
fi

awk '{ print $2 }' "$save_dir/tcga_id.txt" | sort | uniq -d > "$save_dir/tcga_id_duplicated.txt"

while IFS= read -r tcga_id; do
  sample_dir="$save_dir/$tcga_id"
  rm -rf "$sample_dir"
  mkdir -p "$sample_dir"

  zcat "${data_dir}/seven_bridges/processed_vcfs/raw/${tcga_id}.vcf.gz" | grep -v '^##' | sed '1d' | awk '{ print $1 "_" $2 }' > "$sample_dir/all_variants.txt"

  rm -f "$sample_dir/tcga_true_variants_dup.txt"
  while IFS=' ' read -r file_id sample_tcga; do
    if [[ "$sample_tcga" == "$tcga_id" ]]; then
      zcat "${data_dir}/tcga/true_somatic_mutations/${file_id}"/*.maf.gz | grep -v '^#' | sed '1d' | awk '{ print $5 "_" $6 }' >> "$sample_dir/tcga_true_variants_dup.txt"
    fi
done < "$save_dir/tcga_id.txt"

  sort -u "$sample_dir/tcga_true_variants_dup.txt" > "$sample_dir/tcga_true_variants.txt"
  rm -f "$sample_dir/tcga_true_variants_dup.txt"

  rm -f "$sample_dir/all_variants_if_match.txt"
  while IFS= read -r var; do
    if grep -qx "$var" "$sample_dir/tcga_true_variants.txt"; then
      echo 1 >> "$sample_dir/all_variants_if_match.txt"
    else
      echo 0 >> "$sample_dir/all_variants_if_match.txt"
    fi
done < "$sample_dir/all_variants.txt"

  paste "$sample_dir/all_variants.txt" "$sample_dir/all_variants_if_match.txt" > "$sample_dir/all_variants_matched.txt"
  echo "Processed duplicated sample $tcga_id"
done < "$save_dir/tcga_id_duplicated.txt"
