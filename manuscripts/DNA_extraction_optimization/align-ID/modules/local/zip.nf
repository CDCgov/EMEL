process ZIP {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq.gz')    , optional:true, emit: zipped_reads

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_to_zip = "$reads"
    """
    gzip --force $reads_to_zip
    """
}
