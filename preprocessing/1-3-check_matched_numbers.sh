#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
save_dir="${SAVE_DIR:-$repo_root/match_ture_tcga}"

if [[ ! -f "$save_dir/tcga_id.txt" ]]; then
  echo "Missing $save_dir/tcga_id.txt" >&2
  exit 1
fi

rm -f "$save_dir/matched_summary.txt"
while IFS=' ' read -r _ tcga_id; do
  called=$(wc -l < "$save_dir/$tcga_id/all_variants.txt")
  true=$(wc -l < "$save_dir/$tcga_id/tcga_true_variants.txt")
  matched=$(grep -c '^1$' "$save_dir/$tcga_id/all_variants_if_match.txt")
  echo "$tcga_id $called $true $matched" >> "$save_dir/matched_summary.txt"
done < <(awk '{ print $2 }' "$save_dir/tcga_id.txt" | sort | uniq)

echo "Wrote $save_dir/matched_summary.txt"
