process SAMTOOLS_FASTQ {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.19.2--h50ea8bc_0' :
        'biocontainers/samtools:1.19.2--h50ea8bc_0' }"

    input:
    tuple val(meta), path(input)
    val(interleave)

    output:
    tuple val(meta), path("*.fastq")      , optional:true, emit: fastq
    path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output = "> ${prefix}_filtered.fastq"

    """
    samtools \\
        fastq \\
        $args \\
        --threads 4 \\
        $input \\
        $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}