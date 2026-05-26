process NANOPLOT {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.41.6--pyhdfd78af_0' :
        'biocontainers/nanoplot:1.41.6--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    //tuple val(meta), path("*report.html")                , emit: html
    //tuple val(meta), path("*dot.png") , optional: true, emit: dot
    //tuple val(meta), path("*kde.png") , optional: true, emit: kde
    tuple val(meta), path("*.txt")                 , emit: txt
    //tuple val(meta), path("*.log")                 , emit: log
    //path  "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith(".fastq.gz")) ? "--fastq ${ontfile}" : ("$ontfile".endsWith(".fastq")) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith(".txt")) ? "--summary ${ontfile}" : ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    NanoPlot \\
        $args \\
        -t 4 \\
        --prefix "${prefix}_" \\
        --fastq ${ontfile}

    if [[ "$input_file" == *"chopper"* ]]; then

        mv ${prefix}_NanoStats.txt ${prefix}_NanoStats_chopper.txt

    elif [[ "$input_file" == *".trim."* ]]; then

        mv ${prefix}_NanoStats.txt ${prefix}_NanoStats_chopper-cutadapt.txt
    
    elif [[ "$input_file" == *"filtered"* ]]; then

        mv ${prefix}_NanoStats.txt ${prefix}_NanoStats_mapq-filtered.txt

    fi

    if [[ -f ${prefix}_NanoStats_post_filtering.txt ]]; then

        rm ${prefix}_NanoStats_post_filtering.txt
    
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """

}
