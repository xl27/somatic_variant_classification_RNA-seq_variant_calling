#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.repo_root = params.repo_root ?: "${baseDir}"

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
    MODELING()
}
