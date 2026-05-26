process SAMTOOLS_SORT_INDEX {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(bam)
    val(interleave)

    output:
    tuple val(meta), path("*sorted.bam")                              , optional:true, emit: bam
    tuple val(meta), path("*sorted.bam.bai")                          , optional:true, emit: bai
    tuple val(meta), path("*sorted.bam"), path("*sorted.bam.bai")     , optional:true, emit: bam_bai
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "> ${prefix}_quality_filtered.bam"

    """
    samtools sort $bam > ${prefix}_quality_filtered_sorted.bam

    samtools index ${prefix}_quality_filtered_sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}