process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3a70f8bc7e17b723591f6132418640cfdbc88246-0' :
        'biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3a70f8bc7e17b723591f6132418640cfdbc88246-0' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.csi"), optional: true, emit: csi
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_output = " | samtools sort -@ ${task.cpus-1} -o ${prefix}.bam ${args2}"
    """
    if ${params.memory_saver}; then
    
        start="false"

        while [[ \$start == "false" ]]; do

            if [[ \$(ls ${projectDir}/queue/minimap2 | wc -l) -gt 1 ]]; then

                sleep 5
            
            else

                start="true"
                touch ${projectDir}/queue/minimap2/$prefix-minimap2

            fi

        done
    fi

    if [[ ${reference} == *".mmi" ]]; then

        minimap2 \\
            $args -a ${reference} ${reads} \\
            $bam_output

    else

        minimap2 \\
            $args -ax map-ont ${reference} ${reads} \\
            $bam_output

    fi    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS

    if ${params.memory_saver}; then

        rm -f ${projectDir}/queue/minimap2/$prefix-minimap2
    
    fi

    """
}