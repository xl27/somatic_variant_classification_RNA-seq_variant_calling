#!/usr/bin/env Rscript
repo_root <- normalizePath(getwd())
annotations_dir <- Sys.getenv('ANNOTATIONS_DIR', unset = file.path(repo_root, 'annotations'))
sample_files <- read.table(file.path(repo_root, 'features', 'samples_id.txt'), stringsAsFactors = FALSE)$V1

for (sample_id in sample_files) {
  input_file <- file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_max1e-4.csv'))
  exons_file <- file.path(annotations_dir, 'exons', paste0(sample_id, '.variants.exons.txt'))
  if (!file.exists(input_file) || !file.exists(exons_file)) {
    message('Skipping ', sample_id, ' because required file is missing.')
    next
  }

  features_df <- read.csv(input_file, row.names = 1)
  exons_var <- read.table(exons_file, stringsAsFactors = FALSE)$V1
  features_fil <- features_df[rownames(features_df) %in% exons_var, ]
  output_file <- file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_exons_jointAF1e-4.csv'))
  write.csv(features_fil, output_file, row.names = TRUE, quote = FALSE)
  message('Filtered ', sample_id)
}
