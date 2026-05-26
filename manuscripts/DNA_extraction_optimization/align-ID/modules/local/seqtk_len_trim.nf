process SEQTK_LEN_TRIM {
    tag "$meta.id"
    label 'process_low'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/scicomp/groups-pure/WDPB/EMEL/singularity/seqtk/seqtk%3A1.4--he4a0461_2' :
        'https://depot.galaxyproject.org/singularity/seqtk%3A1.4--he4a0461_2' }"

    input:
    tuple val(meta), path(fastq)
    val(length)

    output:
    tuple val(meta), path('*fastq.gz'), emit: fastq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (meta.single_end) {

    """
    
    seqtk seq -L $length $fastq > ${prefix}_length_filtered.fastq
    gzip ${prefix}_length_filtered.fastq

    """

    } else {

    """
    
    seqtk seq -L $length ${fastq[0]} > ${prefix}_R1_length_filtered.fastq
    gzip ${prefix}_R1_length_filtered.fastq

    seqtk seq -L $length ${fastq[1]} > ${prefix}_R2_length_filtered.fastq
    gzip ${prefix}_R2_length_filtered.fastq

    """        

    }
}
