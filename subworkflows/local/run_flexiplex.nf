//
// Creates gtfs to that add introns as features
//

include { PIGZ_UNCOMPRESS } from '../../modules/nf-core/pigz/uncompress/main'
include { PIGZ_COMPRESS } from '../../modules/nf-core/pigz/compress/main'
include { FLEXIPLEX       } from '../../modules/local/flexiplex/main'

workflow RUN_FLEXIPLEX {
    take:
        reads
        adapter_5prime
        adapter_3prime
        spacer_3prime
        ed_adapter_5prime
        ed_adapter_3prime
        barcode_length
        umi_length

    main:
        ch_versions = Channel.empty()

		//
        // Check if reads are zipped
        //
        gzipped = reads.map { meta, fastq_list -> 
            def all_gzipped = fastq_list.every { it.endsWith('.gz') }
            def none_gzipped = fastq_list.every { !it.endsWith('.gz') }

            if (!all_gzipped && !none_gzipped) {
                throw new IllegalArgumentException("Error: Mixed gzipped and non-gzipped files in ${fastq_list}")
            }
            
            return all_gzipped
        }

    
        ch_prepared_fasta = Channel.empty()
        if (gzipped){
            PIGZ_UNCOMPRESS( reads )

            ch_reads = PIGZ_UNCOMPRESS.out.file
            ch_versions = ch_versions.mix(PIGZ_UNCOMPRESS.out.versions)
        } else {
            ch_reads = reads
        }


        //
        // MODULE: Run flexiplex
        //
        FLEXIPLEX (
            ch_reads,
            params.adapter_5prime,
            params.adapter_3prime,
            params.spacer_3prime,
            params.ed_adapter_5prime,
            params.ed_adapter_3prime,
            params.barcode_length,
            params.umi_length
    	)
        
        //
        // MODULE: Compress Fastq
        //
        PIGZ_COMPRESS ( FLEXIPLEX.out.reads )
        
        
    emit:
        flexiplex_fastq = PIGZ_COMPRESS.out.archive
        versions = ch_versions
}