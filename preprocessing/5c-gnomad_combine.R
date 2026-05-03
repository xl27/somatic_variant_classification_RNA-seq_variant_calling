#!/usr/bin/env Rscript
repo_root <- normalizePath(getwd())
gnomad_dir <- Sys.getenv('GNOMAD_DIR', unset = file.path(repo_root, 'gnomAD'))
sample_files <- read.table(file.path(repo_root, 'features', 'samples_id.txt'), stringsAsFactors = FALSE)$V1
chr_ls <- paste0('chr', c(1:22, 'X', 'Y'))

for (sample_id in sample_files) {
  features_chr_ls <- lapply(chr_ls, function(chr) {
    path <- file.path(gnomad_dir, 'temp', chr, paste0(sample_id, '.joint_AF_pop.csv'))
    if (!file.exists(path)) {
      stop('Missing file: ', path)
    }
    read.csv(path, row.names = 1, stringsAsFactors = FALSE)
  })
  features_df <- do.call(rbind, features_chr_ls)
  output_dir <- file.path(gnomad_dir, 'features_jointAFpop')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(features_df, file.path(output_dir, paste0(sample_id, '.joint_AF_pop.csv')), row.names = TRUE, quote = FALSE)
  message('Combined ', sample_id)
}
