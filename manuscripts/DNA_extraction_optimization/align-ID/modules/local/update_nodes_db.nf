process UPDATE_NODES_DB {
    label 'process_medium'

    input:
    path(local_nodes_db)
    path(ncbi_nodes)
    path(my_tax_ids)

    output:
    path('complete.txt') , emit: complete

    script:

    """

    bash "${projectDir}/bin/update_nodes_db.sh" $local_nodes_db $ncbi_nodes $my_tax_ids

    """

}