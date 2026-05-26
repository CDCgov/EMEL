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


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/long-read/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/long-read/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check_long_read'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { PORECHOP_PORECHOP                              } from '../modules/nf-core/porechop/porechop/main'
include { CHOPPER                                        } from '../modules/nf-core/chopper/main'
include { CUTADAPT                                      } from '../modules/nf-core/cutadapt/main'
include { SAMTOOLS_DEPTH                                 } from '../modules/nf-core/samtools/depth/main'
include { SAMTOOLS_COVERAGE                              } from '../modules/nf-core/samtools/coverage/main'
include { MULTIQC                                        } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                    } from '../modules/nf-core/custom/dumpsoftwareversions/main'


//
// MODULE: custom, local modules
//
include { UPDATE_NODES_DB                                 } from '../modules/local/update_nodes_db'
include { BBMAP_REFORMAT as BBMAP_REFORMAT_SUBSAMPLE      } from '../modules/local/bbmap_reformat_subsample'
include { NANOPLOT as NANOPLOT_RAW                        } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_TRIMMED                    } from '../modules/local/nanoplot'
include { NANOPLOT as NANOPLOT_ALIGNMENT_TAXON_FILTERED   } from '../modules/local/nanoplot'
include { MINIMAP2_ALIGN as ALIGN_READS                   } from '../modules/local/minimap2_long'
include { SAMTOOLS_STATS                                  } from '../modules/local/samtools_stats'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED       } from '../modules/local/samtools_fastq'
include { ALIGNMENT_CLASSIFY                              } from '../modules/local/alignment_classify'
include { BLAST_UNMAPPED_READS                            } from '../modules/local/blast_unmapped_reads'
include { SAMTOOLS_SORT_INDEX                             } from '../modules/local/samtools_sort_index'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../modules/local/samtools_fastq'
include { ALIGNMENT_CLASSIFICATION_GRAPH as ALIGNMENT_CLASSIFICATION_GRAPH_READS } from '../modules/local/alignment_classification_graph'

//clearing out minimap2 queues if memory_saver mode is enabled (only required for local compute; memory allocation should generally be handled by the job scheduler when submitting to the cluster)
if (params.memory_saver) {

    def minimap2_queue = new File("${projectDir}/queue/minimap2")

    if (minimap2_queue.exists() && minimap2_queue.isDirectory()) {
        minimap2_queue.eachFile { file ->
            file.delete()
        }
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

ch_versions = Channel.empty()
ch_multiqc_files = Channel.empty()

workflow LONG_READ_ID {

    if (params.filter_alignment_by_id) {

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

    //System.exit(1)
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    if (!params.skip_subsample) {

        BBMAP_REFORMAT_SUBSAMPLE (
            INPUT_CHECK.out.reads,
            params.num_subsamples
        )

        raw_reads = BBMAP_REFORMAT_SUBSAMPLE.out.fastq

    } else {

        raw_reads = INPUT_CHECK.out.reads

    }

    NANOPLOT_RAW (

        raw_reads

    )

    if (params.skip_porechop) {

        trimmed_reads = raw_reads

    } else {

        PORECHOP_PORECHOP (

            raw_reads

        )

        trimmed_reads = PORECHOP_PORECHOP.out.reads

    }

    CHOPPER (

        trimmed_reads

    )

    if (!params.skip_cutadapt) {

        CUTADAPT (

            CHOPPER.out.fastq

        )

        trimmed_reads = CUTADAPT.out.reads

    } else {

        trimmed_reads = CHOPPER.out.fastq

    }

    NANOPLOT_TRIMMED (

        trimmed_reads

    )

    ALIGN_READS (

        trimmed_reads,
        [[params.minimap2_meta],[params.minimap2_genome_reference]]

    )

    // custom module that parses a bam file and pulls out each alignment that corresponds to a reference sequence under a specified tax id
    ALIGNMENT_CLASSIFY (

        ALIGN_READS.out.bam,
        params.seqid2taxid_map,
        params.filter_alignment_by_id,
        params.my_tax_ids,
        params.include_children

    )

    BLAST_UNMAPPED_READS (

        ALIGNMENT_CLASSIFY.out.unmapped_bam,
        params.blast_db

    )

    if (params.filter_alignment_by_id) {

        alignment_classified_bam = ALIGNMENT_CLASSIFY.out.classified_plus_filtered_bam

    } else {

        alignment_classified_bam = ALIGNMENT_CLASSIFY.out.primary_all

    }

    SAMTOOLS_STATS (

        alignment_classified_bam

    )

    SAMTOOLS_DEPTH (

        alignment_classified_bam

    )

    SAMTOOLS_SORT_INDEX (

        alignment_classified_bam,
        false

    )

    SAMTOOLS_COVERAGE (

        alignment_classified_bam.join(SAMTOOLS_SORT_INDEX.out.bai)

    )

    // capturing aligned reads and converting to fastq
    SAMTOOLS_FASTQ_MAPPED (

        SAMTOOLS_SORT_INDEX.out.bam,
        false

    )

    filtered_reads = SAMTOOLS_FASTQ_MAPPED.out.fastq

    // running nanoplot again to compare read stat pre and post filter
    NANOPLOT_ALIGNMENT_TAXON_FILTERED (

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

    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_RAW.out.txt.collect{it[1]}.ifEmpty([]))
    if (!params.skip_porechop) {
        ch_multiqc_files = ch_multiqc_files.mix(PORECHOP_PORECHOP.out.log.collect().ifEmpty([]))
    }

    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.txt.map{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_DEPTH.out.tsv.map{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_COVERAGE.out.coverage.map{it[1]})
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_ALIGNMENT_TAXON_FILTERED.out.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ALIGNMENT_CLASSIFICATION_GRAPH_READS.out.plot.collect{it[1]}.ifEmpty([]))

    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT_TRIMMED.out.txt.collect{it[1]}.ifEmpty([]))

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
