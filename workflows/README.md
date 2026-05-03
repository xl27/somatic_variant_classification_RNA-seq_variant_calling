# Pipeline workflows

This folder contains Nextflow workflow definitions for the somatic discovery pipeline.

- `generate_features_pipeline.nf`: builds TCGA matched variant sets, deduplicates, checks summary counts, and processes 1000 Genomes.
- `post_filtering_pipeline.nf`: runs feature filtering, population allele-frequency annotation, exon filtering, VAF augmentation, and final combination.
- `modeling_pipeline.nf`: runs downstream model training using the combined feature set.

Run from the repository root:

```bash
cd Github
nextflow run workflows/generate_features_pipeline.nf -resume
nextflow run workflows/post_filtering_pipeline.nf -resume
nextflow run workflows/modeling_pipeline.nf -resume
```
