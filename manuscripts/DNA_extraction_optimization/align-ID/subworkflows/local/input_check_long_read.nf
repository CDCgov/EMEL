//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check_long_read'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    complete // placeholder file to prevent INPUT_CHECK from running until UPDATE_NODES_DB has finished

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
    versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

// Function to get list of [ meta, [ fastq_long, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = true
    meta.short_1 = row.short_1
    meta.short_2 = row.short_2

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    if (!file(row.fastq_long).exists()) {
        exit 1, "ERROR: Please check input samplesheet -> FastQ file does not exist!\n${row.fastq_long}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ file(row.fastq_long) ] ]
    }
    return fastq_meta
}
