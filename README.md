# Somatic Variant Classification from RNA-seq Variant Calling using Tumor-only Samples

## Structure

- `preprocessing/`: preprocessing pipeline scripts
- `modeling/`: machine learning model training and evaluation scripts for somatic variant prediction
- `workflows/`: Nextflow workflow definitions for feature generation, filtering and modeling

## Preprocessing overview

The preprocessing stage prepares the somatic variant feature set for modeling:

- Identify TCGA sample IDs and match somatic variants to true TCGA calls.
- Remove duplicate sample mappings and verify matched sample counts.
- Process 1000 Genomes variant data for germline filtering.
- Generate per-sample feature tables from filtered variant calls.
- Apply genotype, depth, mapping quality, AF, and VCF filters.
- Annotate variants with gnomAD joint allele frequencies.
- Convert gene/exon annotations into exon interval filters and apply exon-based variant filtering.
- Add VAF/CN-derived neighbor features and combine all per-sample results into a final deduplicated dataset.

## Example run

```bash
cd Github
bash preprocessing/0-find_tcga_sample_id.sh
bash preprocessing/1-1-match_true_somatic_mutations.sh
bash preprocessing/1-2-match_true_duplicated_samples.sh
bash preprocessing/1-3-check_matched_numbers.sh
bash preprocessing/2-process_oneKG.sh
bash preprocessing/3-process_features.sh
bash preprocessing/4-filter_new.sh
bash preprocessing/5b-gnomad_match_AF_chr.sh
bash preprocessing/5c-gnomad_combine.sh
bash preprocessing/6a-annotations_convert.sh
bash preprocessing/6b-annotations_extract_exons.sh
bash preprocessing/6c-annotations_filter_exons.sh
bash preprocessing/7-add_vaf_cn.sh
python3 preprocessing/8-combine_nodup.py
```

## Run workflows

```bash
cd Github
nextflow run workflows/generate_features_pipeline.nf -resume
nextflow run workflows/post_filtering_pipeline.nf -resume
nextflow run workflows/modeling_pipeline.nf -resume
```

## Model Training

The `modeling/` folder contains XGBoost model training scripts with SMOTE oversampling and class weighting strategies. Scripts are organized by filtering approach:

- `*_noSampleFilter.py`: No additional sample filtering
- `*_Filter-num20.py`: Samples with >20 somatic mutations
- `*_Filter-tcga.py`: TCGA quantile-based filtering (0.1-0.9)
- `*_perDisease.py`: Disease-specific training

Run scripts from the repository root with appropriate environment variables set for data paths.

## Notes

- Data directories are intentionally excluded from this repo.
- Use environment variables like `DATA_DIR`, `SAVE_DIR`, `ONEKG_DIR`, `ANNOTATIONS_DIR`, `GNOMAD_DIR`, `TCGA_INFO_DIR`, `FEATURES_NEW_DIR` to configure paths.
- `preprocessing/5a-gnomad_download.sh` is optional and can be used to download raw gnomAD VCFs.
