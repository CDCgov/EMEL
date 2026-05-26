process ALIGNMENT_CLASSIFY {
    tag "$meta.id"
    label 'process_high'
    errorStrategy 'ignore'

    input:
    tuple val(meta), path(bam)
    path(seq2tax_map)
    val(filter_alignment_by_id)
    path(my_tax_ids)
    val(include_children)


    output:
    tuple val(meta), path('*unmapped.bam') , optional:true, emit: unmapped_bam
    tuple val(meta), path('*-classified-plus-filtered-sorted.bam') , optional:true, emit: classified_plus_filtered_bam
    tuple val(meta), path('*summary.tsv') , optional:true, emit: summary_tsv
    tuple val(meta), path('*primary-all-sorted.bam') , optional:true, emit: primary_all
    tuple val(meta), path('*primary_unambiguous-sorted.bam') , optional:true, emit: primary_unambiguous
    tuple val(meta), path('*primary_ambiguous_single_genome-sorted.bam') , optional:true, emit: primary_ambiguous_single_genome
    tuple val(meta), path('*primary_ambiguous_multi_genome-sorted.bam') , optional:true, emit: primary_ambiguous_multi_genome    

    script:
    def prefix = "${meta.id}"

    if (params.workflow == "long-read") {
    
        """

        bash "${projectDir}/bin/alignment_classify_SE.sh" $prefix $bam $seq2tax_map $filter_alignment_by_id $my_tax_ids $include_children ${params.ncbi_taxonomy_nodes} ${params.ncbi_taxonomy_names} ${projectDir} ${params.local_nodes_db} ${params.non_standard_reference} ${params.mapping_quality} ${meta.single_end}

        """

    } else if (params.workflow == "short-read") {

        """
        if [ ${meta.single_end} == 'true' ]; then

            bash "${projectDir}/bin/alignment_classify_SE.sh" $prefix $bam $seq2tax_map $filter_alignment_by_id $my_tax_ids $include_children ${params.ncbi_taxonomy_nodes} ${params.ncbi_taxonomy_names} ${projectDir} ${params.local_nodes_db} ${params.non_standard_reference} ${params.mapping_quality} ${meta.single_end}

        else

            bash "${projectDir}/bin/alignment_classify_PE.sh" $prefix $bam $seq2tax_map $filter_alignment_by_id $my_tax_ids $include_children ${params.ncbi_taxonomy_nodes} ${params.ncbi_taxonomy_names} ${projectDir} ${params.local_nodes_db} ${params.non_standard_reference} ${params.mapping_quality} ${meta.single_end}

        fi

        """

    }
}