/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

//skip the line below to prevent "--fasta not specified" errors
//WorkflowMetashort.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/short-read/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/short-read/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check_short_read'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC as FASTQC_RAW                           } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED                       } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_ALIGNED_FILTERED              } from '../modules/nf-core/fastqc/main'
include { SAMTOOLS_DEPTH                                 } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE                              } from '../modules/nf-core/samtools/coverage/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// MODULE: custom, local modules
//
include { UPDATE_NODES_DB                                 } from '../modules/local/update_nodes_db'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_SUBSAMPLE      } from '../modules/local/bbmap_reformat_subsample'
include { TRIMMOMATIC                                     } from '../modules/local/trimmomatic'
include { BBMAP_BBMERGE                                   } from '../modules/local/bbmap_merge'
include { FASTP                                           } from '../modules/local/fastp'
include { MINIMAP2_ALIGN as ALIGN_READS_MINIMAP2          } from '../modules/local/minimap2_short'
include { ALIGNMENT_CLASSIFY                              } from '../modules/local/alignment_classify'
include { BLAST_UNMAPPED_READS                            } from '../modules/local/blast_unmapped_reads'
include { SAMTOOLS_STATS                                  } from '../modules/local/samtools_stats'
include { SAMTOOLS_SORT_INDEX                             } from '../modules/local/samtools_sort_index'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../modules/local/samtools_fastq'
include { ALIGNMENT_CLASSIFICATION_GRAPH as ALIGNMENT_CLASSIFICATION_GRAPH_READS } from '../modules/local/alignment_classification_graph'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//clearing out kraken and minimap2 queues if memory_saver mode is enabled (only required for local compute; memory allocation should generally be handled by the job scheduler when using the cluster)
if (params.memory_saver) {

    def kraken_queue = new File("${projectDir}/queue/kraken")

    if (kraken_queue.exists() && kraken_queue.isDirectory()) {
        kraken_queue.eachFile { file ->
            file.delete()
        }
    }

    def minimap2_queue = new File("${projectDir}/queue/minimap2")

    if (minimap2_queue.exists() && minimap2_queue.isDirectory()) {
        minimap2_queue.eachFile { file ->
            file.delete()
        }
    }

}


// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow SHORT_READ_ID {

    if (!params.skip_alignment_based_filtering && params.filter_alignment_by_id) {

        UPDATE_NODES_DB (

            params.local_nodes_db,
            params.ncbi_taxonomy_nodes,
            params.my_tax_ids

        )

        placeholder = UPDATE_NODES_DB.out.complete

    } else {

        placeholder = []

    }

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input),
        placeholder
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // MODULE: Run FastQC
    //

    if (!params.skip_subsample) {

        BBMAP_REFORMAT_SUBSAMPLE (
            INPUT_CHECK.out.reads,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT_SUBSAMPLE.out.fastq

    } else {

        raw_reads = INPUT_CHECK.out.reads

    }

    FASTQC_RAW (
        raw_reads
    )
    ch_versions = ch_versions.mix(FASTQC_RAW.out.versions.first())


    if (!params.skip_trimming) {

        if (params.trim_tool == "fastp") {

            if (!params.adapter_auto_detect) {

                FASTP (
                    raw_reads,
                    ["${params.adapt_ref}"],
                    [],
                    [],
                    []
                )

            } else {

                FASTP (
                    raw_reads,
                    [],
                    [],
                    [],
                    []
                )

            }

            trimmed_reads = FASTP.out.reads

        } else if (params.trim_tool == "trimmomatic") {

            TRIMMOMATIC (

                raw_reads

            )

            trimmed_reads = TRIMMOMATIC.out.trimmed_reads

        }

        if (!params.skip_merging) {

            BBMAP_BBMERGE (
                trimmed_reads,
                false
            )

            trimmed_reads = BBMAP_BBMERGE.out.merged
        }

        FASTQC_TRIMMED (

            trimmed_reads

        )

        FASTQC_CH = FASTQC_TRIMMED
        ch_versions = ch_versions.mix(FASTQC_TRIMMED.out.versions.first())

    } else {

        FASTQC_CH = FASTQC_RAW
        trimmed_reads = raw_reads

    }


    ALIGN_READS_MINIMAP2 (

        trimmed_reads,
        [[params.minimap2_meta],[params.minimap2_genome_reference]]

    )

    aligner_bam_ch = ALIGN_READS_MINIMAP2.out.bam

    ALIGNMENT_CLASSIFY (

        aligner_bam_ch,
        params.seqid2taxid_map,
        params.filter_alignment_by_id,
        params.my_tax_ids,
        params.include_children

    )

    if (!params.skip_blast_unmapped) {

        BLAST_UNMAPPED_READS (
            ALIGNMENT_CLASSIFY.out.unmapped_bam,
            params.blast_db
        )

    }

    if (params.filter_alignment_by_id) {

        alignment_classified_bam = ALIGNMENT_CLASSIFY.out.classified_plus_filtered_bam

    } else {

        alignment_classified_bam = ALIGNMENT_CLASSIFY.out.primary_all

    }

    SAMTOOLS_STATS (

        alignment_classified_bam

    )

    SAMTOOLS_SORT_INDEX (

        alignment_classified_bam,
        false

    )

    SAMTOOLS_DEPTH (

        alignment_classified_bam

    )

    SAMTOOLS_COVERAGE (

        alignment_classified_bam.join(SAMTOOLS_SORT_INDEX.out.bai)

    )

    // capturing aligned reads and converting to fastq
    SAMTOOLS_FASTQ_MAPPED (

        alignment_classified_bam,
        false

    )

    filtered_reads = SAMTOOLS_FASTQ_MAPPED.out.fastq

    FASTQC_ALIGNED_FILTERED (

        filtered_reads

    )


    ALIGNMENT_CLASSIFICATION_GRAPH_READS (

        ALIGNMENT_CLASSIFY.out.summary_tsv,
        "Read"

    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowAlignID.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowAlignID.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect { it[1] }.ifEmpty([]))

    if (!params.skip_trimming) {
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect { it[1] }.ifEmpty([]))

        if (params.trim_tool == "trimmomatic") {

            ch_multiqc_files = ch_multiqc_files.mix(TRIMMOMATIC.out.trim_log.collect{it[1]}.ifEmpty([]))

        } else if (params.trim_tool == "fastp") {

            ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect { it[1] }.ifEmpty([]))

        }
    }

    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_CLASSIFICATION_GRAPH_READS.out.plot.collect{it[1]}.ifEmpty([]))

    MULTIQC (

        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()

    )

    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
