/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_novast_pipeline'
include { MINIMAP2_ALIGN                            } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                            } from '../modules/nf-core/minimap2/index/main'
include { SEQKIT_SPLIT2                             } from '../modules/nf-core/seqkit/split2/main'
include { CAT_FASTQ                                 } from '../modules/nf-core/cat/fastq/main'
include { UMITOOLS_DEDUP                            } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRE      } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_POST     } from '../modules/nf-core/samtools/index/main'
include { FLEXIFORMATTER                            } from '../modules/local/flexiformatter/main'

/*
 * Import subworkflows
 */
include { QCFASTQ_NANOPLOT_FASTQC                   } from '../subworkflows/nf-core/toulligqc_nanoplot_fastqc'
include { PREPARE_REFERENCE_FILES                   } from '../subworkflows/local/prepare_reference_files'
include { RUN_FLEXIPLEX                             } from '../subworkflows/local/run_flexiplex'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOVAST {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // SUBWORKFLOW: FastQC, Nanoplot, ToulligQC
    //
    //TODO: Will there be only pretrim QC or posttrim QC as well?
    // If only pretrim, change names for clarity
    ch_fastqc_multiqc_pretrim = Channel.empty()
    if (!params.skip_qc){

        QCFASTQ_NANOPLOT_FASTQC ( ch_samplesheet, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = QCFASTQ_NANOPLOT_FASTQC.out.fastqc_multiqc.ifEmpty([])
    }
    
    //
    // SUBWORKFLOW: Prepare reference files
    //

    PREPARE_REFERENCE_FILES ( "",
                              "",
                              params.fasta,
                              params.gtf )

    fasta = PREPARE_REFERENCE_FILES.out.prepped_fasta
    fai = PREPARE_REFERENCE_FILES.out.prepped_fai
    gtf = PREPARE_REFERENCE_FILES.out.prepped_gtf
    
    
    //
    // MODULE: Run SEQKIT_SPLIT2
    //
    SEQKIT_SPLIT2 (
        ch_samplesheet
    )
    
    // Transpose channel and add part to metadata
    SEQKIT_SPLIT2.out.reads
        | transpose
        | map { meta, reads ->
              part = (reads =~ /.*part_(\d+)\.fastq(?:\.gz)?$/)[0][1]
              newmap = [part: part]
              [meta + newmap, reads] }
        | set { ch_split_fastq }
            
    ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)
    
    //
    // SUBWORKFLOW: RUN_FLEXIPLEX
    //
    RUN_FLEXIPLEX (
        ch_split_fastq,
        params.adapter_5prime,
        params.adapter_3prime,
        params.spacer_3prime,
        params.ed_adapter_5prime,
        params.ed_adapter_3prime,
        params.barcode_length,
        params.umi_length
    )
    
    // Group by ID for CATFASTQ
    RUN_FLEXIPLEX.out.flexiplex_fastq
        | map { meta, reads ->
               [meta.subMap('id', 'single_end'), meta.part, reads] }
        | groupTuple
        | map { meta, part, reads -> [meta + [partcount: part.size()], reads] }
        | set { ch_grouped_flexiplex_fastq }

        //.set { ch_split_flexiplex_fastq }

    ch_versions = ch_versions.mix(RUN_FLEXIPLEX.out.versions)
    
    //
    // MODULE: CATFASTQ
    //  
    CAT_FASTQ (
        ch_grouped_flexiplex_fastq
    )
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions)
    
    // Set channel
    CAT_FASTQ.out.reads
        | set { ch_merged_flexiplex_fastq}
    
    //
    // MODULE: Run MINIMAP2_INDEX
    //
    ch_minimap_index = fasta
    if (!params.skip_save_minimap2_index) {
        
        MINIMAP2_INDEX ( fasta )
        ch_minimap_index = MINIMAP2_INDEX.out.index
        ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions)
    }

    //
    // MODULE: Run MINIMAP2_ALIGN
    //
    MINIMAP2_ALIGN (
        ch_merged_flexiplex_fastq,
        ch_minimap_index,
        true,
        'bai',
        "",
        ""
    )

    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    MINIMAP2_ALIGN.out.bam 
        | set { ch_minimap_bam }

    //
    // MODULE: Run FLEXI_FORMATTER
    //
    FLEXIFORMATTER (
        ch_minimap_bam
    )
    ch_versions = ch_versions.mix(FLEXIFORMATTER.out.versions)
    FLEXIFORMATTER.out.bam
        | set { ch_tagged_bam }
    
    //
    // MODULE: Run SAMTOOLS_INDEX
    //
    SAMTOOLS_INDEX_PRE (
        ch_tagged_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_PRE.out.versions)
    SAMTOOLS_INDEX_PRE.out.bai
        | set { ch_bai }
    
    //
    // MODULE: Run UMITOOLS
    //  
    UMITOOLS (
        ch_tagged_bam.join(ch_bai, by: [0]),
        true
    )
    ch_versions = ch_versions.mix(UMITOOLS.out.versions)
    UMITOOLS.out.bam
        | set { ch_deduped_bam }

    
    //
    // MODULE: Run SAMTOOLS_INDEX_POSTDEDUP
    //
    SAMTOOLS_INDEX_POST (
        ch_deduped_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX_POST.out.versions)
    SAMTOOLS_INDEX_POST.out.bai
        | set { ch_dedup_bai }
    
    ch_deduped_bam.join(ch_dedup_bai, by: [0])
        | set { ch_deduped_bam_bai }
    
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_replace_names = params.multiqc_replace_names ?
        Channel.fromPath(params.multiqc_replace_names, checkIfExists: true) :
        Channel.empty()
    ch_multi_qc_sample_names = params.multiqc_sample_names ?
        Channel.fromPath(params.multiqc_sample_names, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )
    
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        ch_multiqc_replace_names.toList(),
        ch_multi_qc_sample_names.toList()
    )
    
    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
