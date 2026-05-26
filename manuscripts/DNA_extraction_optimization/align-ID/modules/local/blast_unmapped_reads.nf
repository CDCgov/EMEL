process BLAST_UNMAPPED_READS {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(bam)
    path  db

    output:
    tuple val(meta), path('*.txt'), emit: txt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """

    bash "${projectDir}/bin/blast_unmapped_reads.sh" ${prefix} $bam ${params.blast_db} ${params.blast_evalue} ${params.blast_perc_identity} ${params.blast_target_seqs}

    """

}
