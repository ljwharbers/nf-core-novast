/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config                       = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config                = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo                         = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description   = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC as MULTIQC_RAWQC                  } from "../modules/nf-core/multiqc/main"
include { MULTIQC as MULTIQC_FINALQC                } from "../modules/nf-core/multiqc/main"
include { paramsSummaryMap                          } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_novast_pipeline'

include { MINIMAP2_ALIGN                            } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX                            } from '../modules/nf-core/minimap2/index/main'
include { SEQKIT_SPLIT2                             } from '../modules/nf-core/seqkit/split2/main'
include { CAT_FASTQ as CAT_FASTQ_SAMPLE             } from '../modules/nf-core/cat/fastq/main'
include { CAT_FASTQ as CAT_FASTQ_SEQSPLIT           } from '../modules/nf-core/cat/fastq/main'
include { UMITOOLS_DEDUP                            } from '../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                            } from '../modules/nf-core/samtools/index/main'
include { FLEXIFORMATTER                            } from '../modules/local/flexiformatter/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS               } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { NANOCOMP as NANOCOMP_FASTQ                } from '../modules/nf-core/nanocomp/main' 
include { NANOCOMP as NANOCOMP_BAM                  } from '../modules/nf-core/nanocomp/main' 

/*
 * Import subworkflows
 */
include { QCFASTQ_NANOPLOT_FASTQC                   } from '../subworkflows/nf-core/toulligqc_nanoplot_fastqc'
include { PREPARE_REFERENCE_FILES                   } from '../subworkflows/local/prepare_reference_files'
include { RUN_FLEXIPLEX                             } from '../subworkflows/local/run_flexiplex'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../subworkflows/nf-core/bam_sort_stats_samtools/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow NOVAST {
    //TODO: Check the different workflows/processes and make sure they're written in the same way.
    // Currently some are copied from scnanoseq/nanoseq and are coded slightly different

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
     
    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    //TODO: Samplesheet from scnanoseq has cell count but mine doesnt. 
    ch_samplesheet
        .branch{
            meta, fastq ->
                single: fastq.size() == 1
                    return [ meta, fastq.flatten() ]
                multiple: fastq.size() > 1
                    return [ meta, fastq.flatten() ]
        }
        .set { ch_fastqs }
    
    //
    // MODULE: Combine fastqs from the same sample
    //
    CAT_FASTQ_SAMPLE ( ch_fastqs.multiple )
        .reads
        .mix ( ch_fastqs.single )
        .set { ch_cat_fastq }

    ch_versions = ch_versions.mix (CAT_FASTQ_SAMPLE.out.versions.first().ifEmpty(null))
    
    //
    // SUBWORKFLOW: FastQC, Nanoplot, ToulligQC
    //
    //TODO: Will there be only pretrim QC or posttrim QC as well?
    // If only pretrim, change names for clarity
    ch_fastqc_multiqc_pretrim = Channel.empty()
    if (!params.skip_qc){

        QCFASTQ_NANOPLOT_FASTQC ( ch_cat_fastq, params.skip_nanoplot, params.skip_toulligqc, params.skip_fastqc )

        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.nanoplot_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.toulligqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(QCFASTQ_NANOPLOT_FASTQC.out.fastqc_version.first().ifEmpty(null))

        ch_fastqc_multiqc_pretrim = QCFASTQ_NANOPLOT_FASTQC.out.fastqc_multiqc.ifEmpty([])
    }
    
    //
    // MODULE: NanoComp for FastQ files
    //

    ch_nanocomp_fastq_html = Channel.empty()
    ch_nanocomp_fastq_txt = Channel.empty()
    if (!params.skip_qc && !params.skip_fastq_nanocomp) {

        NANOCOMP_FASTQ (
            ch_cat_fastq
                .collect{it[1]}
                .map{
                    [ [ 'id': 'nanocomp_fastq.' ] , it ]
                }
        )

        ch_nanocomp_fastq_html = NANOCOMP_FASTQ.out.report_html
        ch_nanocomp_fastq_txt = NANOCOMP_FASTQ.out.stats_txt

        ch_versions = ch_versions.mix( NANOCOMP_FASTQ.out.versions )

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
        ch_cat_fastq
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
    CAT_FASTQ_SEQSPLIT (
        ch_grouped_flexiplex_fastq
    )
    ch_versions = ch_versions.mix(CAT_FASTQ_SEQSPLIT.out.versions)
    
    // Set channel
    CAT_FASTQ_SEQSPLIT.out.reads
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
    // SUBWORKFLOW: BAM_SORT_STATS_SAMTOOLS
    // 
    BAM_SORT_STATS_SAMTOOLS_MINIMAP (
        ch_tagged_bam,
        fasta )
    
    ch_minimap_sorted_bam = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bam
    ch_minimap_sorted_bai = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.bai

    // these stats go for multiqc
    ch_minimap_sorted_stats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.stats
    ch_minimap_sorted_flagstat = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.flagstat
    ch_minimap_sorted_idxstats = BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.idxstats
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS_MINIMAP.out.versions)
    

    //
    // MODULE: NanoComp for BAM files (unfiltered for QC purposes)
    //
    ch_nanocomp_bam_html = Channel.empty()
    ch_nanocomp_bam_txt = Channel.empty()

    if (!params.skip_qc && !params.skip_bam_nanocomp) {

        NANOCOMP_BAM (
            ch_tagged_bam
                .collect{it[1]}
                .map{
                    [ [ 'id': 'nanocomp_bam.' ] , it ]
                }

        )

        ch_nanocomp_bam_html = NANOCOMP_BAM.out.report_html
        ch_nanocomp_bam_txt = NANOCOMP_BAM.out.stats_txt
        ch_versions = ch_versions.mix( NANOCOMP_BAM.out.versions )
    }
    
    
    //
    // MODULE: Run UMITOOLS_DEDUP
    //  
    UMITOOLS_DEDUP (
        ch_tagged_bam.join(ch_bai, by: [0]),
        true
    )
    ch_versions = ch_versions.mix(UMITOOLS_DEDUP.out.versions)
    UMITOOLS_DEDUP.out.bam
        | set { ch_deduped_bam }

    
    //
    // MODULE: Run SAMTOOLS_INDEX_POSTDEDUP
    //
    SAMTOOLS_INDEX (
        ch_deduped_bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)
    SAMTOOLS_INDEX.out.bai
        | set { ch_dedup_bai }
    
    ch_deduped_bam.join(ch_dedup_bai, by: [0])
        | set { ch_deduped_bam_bai }
    
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
            ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

    if (!params.skip_qc && !params.skip_multiqc){

        //
        // MODULE: MultiQC for raw data
        //

        ch_multiqc_rawqc_files = Channel.empty()
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_config)
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_fastqc_multiqc_pretrim.collect().ifEmpty([]))
        ch_multiqc_rawqc_files = ch_multiqc_rawqc_files.mix(ch_nanocomp_fastq_txt.collect{it[1]}.ifEmpty([]))

        MULTIQC_RAWQC (
            ch_multiqc_rawqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([]),
            [],
            []
        )

        //
        // MODULE: MultiQC for final pipeline outputs
        //
        summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
        ch_workflow_summary    = Channel.value(paramsSummaryMultiqc(summary_params))

        ch_multiqc_finalqc_files = Channel.empty()
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_config)
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

        //ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postrim.collect().ifEmpty([]))
        //ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_fastqc_multiqc_postextract.collect().ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_stats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_minimap_sorted_idxstats.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_rseqc_read_dist.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_nanocomp_bam_txt.collect{it[1]}.ifEmpty([]))

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_tagged_sorted_flagstat.collect{it[1]}.ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_tagged_sorted_idxstats.collect{it[1]}.ifEmpty([]))

        if (!params.skip_dedup) {
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_sorted_flagstat.collect{it[1]}.ifEmpty([]))
            ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_dedup_sorted_idxstats.collect{it[1]}.ifEmpty([]))
        }

        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_read_counts.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_gene_stats_combined.collect().ifEmpty([]))
        ch_multiqc_finalqc_files = ch_multiqc_finalqc_files.mix(ch_transcript_stats_combined.collect().ifEmpty([]))

        MULTIQC_FINALQC (
            ch_multiqc_finalqc_files.collect(),
            ch_multiqc_config,
            ch_multiqc_custom_config.collect().ifEmpty([]),
            ch_multiqc_logo.collect().ifEmpty([]),
            [],
            []
        )
        ch_multiqc_report = MULTIQC_FINALQC.out.report
        ch_versions    = ch_versions.mix(MULTIQC_FINALQC.out.versions)
    }

    emit:
    multiqc_report = ch_multiqc_report.toList()
    versions       = ch_versions
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
