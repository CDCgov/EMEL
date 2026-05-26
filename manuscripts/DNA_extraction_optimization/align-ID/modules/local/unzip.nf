process UNZIP {
    tag "$meta.id"
    label 'process_low'

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path('*')    , optional:true, emit: unzip_contigs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def my_contigs = "$contigs"
    """
    gzip -d --force $my_contigs
    """
}
