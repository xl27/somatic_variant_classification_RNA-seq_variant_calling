#!/usr/bin/env Rscript
repo_root <- normalizePath(getwd())
gnomad_dir <- Sys.getenv('GNOMAD_DIR', unset = file.path(repo_root, 'gnomAD'))
chr <- commandArgs(trailingOnly = TRUE)
if (length(chr) != 1) {
  stop('Usage: Rscript 5b-gnomad_match_AF_chr.R chr')
}

out_dir <- file.path(gnomad_dir, 'temp', chr)
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

input_file <- file.path(gnomad_dir, 'joint', paste0('gnomad.joint.v4.1.sites.', chr, '.AF_joint_pop.tsv'))
if (!file.exists(input_file)) stop('Missing input: ', input_file)

gnomad_chr <- read.delim(input_file, header = FALSE,
                         col.names = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'AF_joint', 'AF_joint_afr', 'AF_joint_ami', 'AF_joint_amr', 'AF_joint_asj', 'AF_joint_eas', 'AF_joint_fin', 'AF_joint_mid', 'AF_joint_nfe', 'AF_joint_remaining', 'AF_joint_sas'))
rownames(gnomad_chr) <- paste(gnomad_chr$CHROM, gnomad_chr$POS, gnomad_chr$REF, gnomad_chr$ALT, sep = '_')

gnomad_chr[, 6:ncol(gnomad_chr)] <- lapply(gnomad_chr[, 6:ncol(gnomad_chr)], function(x) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  x
})

match_chr_AF <- function(features_df) {
  CHROM <- sapply(strsplit(rownames(features_df), '_'), '[[', 1)
  features_chr <- features_df[CHROM == chr, ]
  features_chr_AF <- gnomad_chr[match(rownames(features_chr), rownames(gnomad_chr)), -c(1:5)]
  cbind(if_true_somatic = features_chr$if_true_somatic, features_chr_AF)
}

sample_files <- read.table(file.path(repo_root, 'features', 'samples_id.txt'), stringsAsFactors = FALSE)$V1
for (sample_id in sample_files) {
  input_file <- file.path(repo_root, 'features_new', paste0(sample_id, '_features_filtered_GT_RD_MQ_AF_vcf_exons+jointAF.csv'))
  if (!file.exists(input_file)) {
    message('Skipping ', sample_id, ': missing input file.')
    next
  }
  features_fil <- read.csv(input_file, row.names = 1)
  AF_chr <- match_chr_AF(features_df = features_fil)
  write.csv(AF_chr, file.path(out_dir, paste0(sample_id, '.joint_AF_pop.csv')), row.names = TRUE, quote = FALSE)
  message('Processed ', sample_id)
}
