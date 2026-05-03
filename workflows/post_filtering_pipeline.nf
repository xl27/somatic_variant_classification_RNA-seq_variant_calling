#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.repo_root = params.repo_root ?: "${baseDir}"
params.features_dir = params.features_dir ?: "${params.repo_root}/features"

process FILTER_FEATURES {
    executor 'slurm'
    memory '50 GB'
    time '12h'
    module 'R/4.3.1'

    output:
    path "${params.repo_root}/features_new/*"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/4-filter_new.sh
    """
}

process GNOMAD_FILTER {
    executor 'slurm'
    memory '50 GB'
    time '12h'
    module 'R/4.3.1'

    output:
    path "${params.repo_root}/gnomAD/features_jointAFpop/*"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/5b-gnomad_match_AF_chr.sh
    bash preprocessing/5c-gnomad_combine.sh
    """
}

process ANNOTATION_FILTER {
    executor 'slurm'
    memory '10 GB'
    time '12h'

    output:
    path "${params.repo_root}/features_new/*_features_filtered_GT_RD_MQ_AF_vcf_exons_jointAF1e-4.csv"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/6a-annotations_convert.sh
    bash preprocessing/6b-annotations_extract_exons.sh
    bash preprocessing/6c-annotations_filter_exons.sh
    """
}

process ADD_VAF {
    executor 'slurm'
    module 'R/4.3.1'

    output:
    path "${params.repo_root}/features_new/*_vafs.csv"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/7-add_vaf_cn.sh
    """
}

process COMBINE_NODUP {
    executor 'slurm'
    module 'python/3.10'

    output:
    path "${params.repo_root}/features_new/combined_nodup.csv"

    script:
    """
    cd ${params.repo_root}
    python3 preprocessing/8-combine_nodup.py
    """
}

process MODELING {
    executor 'slurm'
    cpus 4
    memory '100 GB'
    time '1d'

    output:
    path "${params.repo_root}/modeling/*"

    script:
    """
    cd ${params.repo_root}
    bash modeling/jobs.sh
    """
}

workflow {
    FILTER_FEATURES()
    GNOMAD_FILTER()
    ANNOTATION_FILTER()
    ADD_VAF()
    COMBINE_NODUP()
    MODELING()
}
