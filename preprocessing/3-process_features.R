#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop('Usage: Rscript 3-process_features.R TCGA_ID')
}

repo_root <- normalizePath(getwd())
data_dir <- Sys.getenv('DATA_DIR', file.path(repo_root, 'data'))
tcga_id <- args[1]

library(vcfR)
library(BSgenome.Hsapiens.UCSC.hg38)

raw_vcf_path <- file.path(data_dir, 'seven_bridges', 'processed_vcfs', 'raw', paste0(tcga_id, '.vcf.gz'))
match_file <- file.path(repo_root, 'match_ture_tcga', tcga_id, 'all_variants_matched.txt')

raw_variants <- read.vcfR(raw_vcf_path)
raw_variants_fix <- as.data.frame(raw_variants@fix)
rownames(raw_variants_fix) <- paste(raw_variants_fix$CHROM, raw_variants_fix$POS, sep = ">")
raw_variants_fix$POS <- as.numeric(raw_variants_fix$POS)

match_true <- read.table(match_file, col.names = c('locus', 'if_true_somatic'))
gt <- sapply(strsplit(raw_variants@gt[,2], ':'), '[[', 1)
raw_variants_df <- cbind(if_true_somatic = match_true$if_true_somatic, raw_variants_fix[, c(1:2, 4:6)], gt = gt)
raw_variants_df <- raw_variants_df[(nchar(raw_variants_df$REF) == 1) & (nchar(raw_variants_df$ALT) == 1), ]
raw_variants_df$trinucleotide_alt <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, raw_variants_df$CHROM, start = raw_variants_df$POS - 1, end = raw_variants_df$POS + 1))
substr(raw_variants_df$trinucleotide_alt, 2, 2) <- raw_variants_df$ALT
raw_variants_df <- raw_variants_df[, -c(2, 3)]

position_features <- read.csv(file.path(data_dir, 'seven_bridges', 'position_level', paste0(tcga_id, '_position_level.csv')))
rownames(position_features) <- position_features$locus
position_features <- position_features[rownames(raw_variants_df), ]

read_features <- read.csv(file.path(data_dir, 'seven_bridges', 'read_level', paste0(tcga_id, '_read_level.csv')))
read_features <- read_features[read_features$locus %in% rownames(raw_variants_df), ]

position_features_depth <- t(sapply(position_features$X21_bp_depth, function(x) as.numeric(strsplit(x, ',')[[1]])))
rownames(position_features_depth) <- position_features$locus
colnames(position_features_depth) <- paste('read_depth', 1:21, sep = '-')
position_features$ref_freq <- sapply(seq_len(nrow(position_features)), function(i) {
  position_features[i, raw_variants_df[i, 'REF']] / (position_features[i, raw_variants_df[i, 'REF']] + position_features[i, raw_variants_df[i, 'ALT']])
})
position_features_summary <- as.data.frame(cbind(ref_freq = position_features$ref_freq, position_features_depth))

read_features_21bp_qual <- t(sapply(read_features$X21_bp_quality, function(x) strsplit(x, ',')[[1]]))
read_features_21bp_qual <- apply(read_features_21bp_qual, 2, function(x) as.numeric(gsub('n', 0, x)))
rownames(read_features_21bp_qual) <- read_features$locus

read_features_summary <- as.data.frame(cbind(
  apply(read_features_21bp_qual, 2, function(x) tapply(x, read_features$locus, mean)),
  apply(read_features_21bp_qual, 2, function(x) tapply(x, read_features$locus, sd))
))
colnames(read_features_summary) <- paste(rep(c('qual_mean', 'qual_sd'), each = 21), rep(1:21, 2), sep = '-')

read_features_mq_mean <- tapply(read_features$mapping_quality, list(read_features$locus, read_features$base), mean)
for (i in seq_len(ncol(read_features_mq_mean))) {
  read_features_mq_mean[is.na(read_features_mq_mean[, i]), i] <- 0
}
read_features_summary$ref_mapping_qual_mean <- sapply(seq_len(nrow(read_features_mq_mean)), function(i) read_features_mq_mean[i, raw_variants_df[rownames(read_features_mq_mean)[i], 'REF']])
read_features_summary$alt_mapping_qual_mean <- sapply(seq_len(nrow(read_features_mq_mean)), function(i) read_features_mq_mean[i, raw_variants_df[rownames(read_features_mq_mean)[i], 'ALT']])

read_features_mq_sd <- tapply(read_features$mapping_quality, list(read_features$locus, read_features$base), sd)
for (i in seq_len(ncol(read_features_mq_sd))) {
  read_features_mq_sd[is.na(read_features_mq_sd[, i]), i] <- 0
}
read_features_summary$ref_mapping_qual_sd <- sapply(seq_len(nrow(read_features_mq_sd)), function(i) read_features_mq_sd[i, raw_variants_df[rownames(read_features_mq_sd)[i], 'REF']])
read_features_summary$alt_mapping_qual_sd <- sapply(seq_len(nrow(read_features_mq_sd)), function(i) read_features_mq_sd[i, raw_variants_df[rownames(read_features_mq_sd)[i], 'ALT']])

read_features_rev <- tapply(read_features$reverse, read_features$locus, function(x) sum(x == 'True') / length(x))
read_features_summary$reverse_ratio <- read_features_rev
read_features_summary <- read_features_summary[rownames(raw_variants_df), ]
features_df <- cbind(raw_variants_df, position_features_summary, read_features_summary)

pos_df <- as.data.frame(do.call(rbind, strsplit(rownames(features_df), split = '>')))
pos_df$V2 <- as.numeric(pos_df$V2)
ref_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, pos_df$V1, start = pos_df$V2 - 10, end = pos_df$V2 + 10)
freqs <- alphabetFrequency(ref_seq)
features_df$GC <- (freqs[, 'C'] + freqs[, 'G']) / rowSums(freqs)

features_dir <- file.path(repo_root, 'features')
dir.create(features_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(features_df, file = file.path(features_dir, paste0(tcga_id, '_features_processed.rds')))

sample_ids <- list.files(features_dir, pattern = 'features_processed.rds')
sample_ids <- sub('_features_processed.rds$', '', sample_ids)
write.table(sample_ids, file.path(features_dir, 'samples_id.txt'), row.names = FALSE, col.names = FALSE, quote = FALSE)
