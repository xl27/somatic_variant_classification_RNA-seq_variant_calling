#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.repo_root = params.repo_root ?: "${baseDir}"
params.data_dir = params.data_dir ?: "${params.repo_root}/data"
params.save_dir = params.save_dir ?: "${params.repo_root}/match_ture_tcga"
params.onekg_dir = params.onekg_dir ?: "${params.repo_root}/oneKG_germline"

process FIND_TCGA_SAMPLE_ID {
    executor 'slurm'
    cpus 1
    memory '2 GB'
    time '10m'

    output:
    path "${params.save_dir}/tcga_id.txt" into tcga_ids_file

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/0-find_tcga_sample_id.sh
    cp ${params.save_dir}/tcga_id.txt .
    """
}

process MATCH_TRUE_SOMATIC {
    executor 'slurm'
    cpus 1
    memory '10 GB'
    time '36h'

    input:
    path tcga_ids_file

    output:
    path "${params.save_dir}/**"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/1-1-match_true_somatic_mutations.sh
    """
}

process MATCH_DUPLICATES {
    executor 'slurm'
    cpus 1
    memory '1 GB'
    time '4h'

    output:
    path "${params.save_dir}/**"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/1-2-match_true_duplicated_samples.sh
    """
}

process CHECK_MATCHED_NUMBERS {
    executor 'slurm'
    cpus 1
    memory '1 GB'
    time '10m'

    output:
    path "${params.save_dir}/matched_summary.txt" into matched_summary_file

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/1-3-check_matched_numbers.sh
    cp ${params.save_dir}/matched_summary.txt .
    """
}

process PROCESS_ONEKG {
    executor 'slurm'
    cpus 1
    memory '8 GB'
    time '2h'

    output:
    path "${params.onekg_dir}/All_SNPs_filtered.txt"
    path "${params.onekg_dir}/All_SNPs_filtered.rds"

    script:
    """
    cd ${params.repo_root}
    bash preprocessing/2-process_oneKG.sh
    cp ${params.onekg_dir}/All_SNPs_filtered.txt .
    cp ${params.onekg_dir}/All_SNPs_filtered.rds .
    """
}

workflow {
    FIND_TCGA_SAMPLE_ID()
    MATCH_TRUE_SOMATIC()
    MATCH_DUPLICATES()
    CHECK_MATCHED_NUMBERS()
    PROCESS_ONEKG()
}
