//
// run pseudo alignment pipelines
//


include { KALLISTO_QUANT } from '../../modules/nf-core/kallisto/quant' 
include { KALLISTO_INDEX } from '../../modules/nf-core/kallisto/index' 
//include { SALMON_INDEX } from '../../modules/nf-core/salmon/index/main'

include { SALMON_INDEX } from '../../modules/local/salmon_index_custom'
include { SALMON_QUANT } from '../../modules/nf-core/salmon/quant'
include { CUSTOM_TX2GENE   } from '../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../modules/nf-core/tximeta/tximport'




workflow RUN_PSEUDO_ALIGNMENT {
    
    take:
    reads      	// channel: [ val(meta), [ reads ] ]
    aligner		 //  string: aligner method
    fasta              
    gtf              
    

    main:

    ch_versions = Channel.empty()


	//
    // create index
    //
	if (aligner =='kallisto') {
		KALLISTO_INDEX(fasta)
	} else if (aligner == 'salmon') {
		SALMON_INDEX(fasta)
	}


	//
    // run quantification
    //
    if (aligner=='kallisto') {
    	KALLISTO_QUANT(reads, KALLISTO_INDEX.out.index, gtf, [], [], [])
        ch_pseudo_results = KALLISTO_QUANT.out.results
        ch_pseudo_multiqc = KALLISTO_QUANT.out.log
        ch_versions = ch_versions.mix(KALLISTO_QUANT.out.versions.first())
    } else if (aligner == 'salmon') {
        def lib_type = 'A'
    	def alignment_mode = false
        SALMON_QUANT( reads, SALMON_INDEX.out.index,  gtf,  fasta, alignment_mode, lib_type )
        ch_pseudo_results = SALMON_QUANT.out.results
        ch_pseudo_multiqc = ch_pseudo_results
        ch_versions = ch_versions.mix(SALMON_QUANT.out.versions.first())
	}


	//
    // sum up transcripts to genes
    //
    CUSTOM_TX2GENE (
        gtf.map { [ [:], it ] },
        ch_pseudo_results.collect{ it[1] }.map { [ [:], it ] },
        aligner,
        "gene_id","gene_name"
    )
    ch_versions = ch_versions.mix(CUSTOM_TX2GENE.out.versions)

    TXIMETA_TXIMPORT (
        ch_pseudo_results.collect{ it[1] }.map { [ ['id': 'all_samples'], it ] },
        CUSTOM_TX2GENE.out.tx2gene,
        aligner
    )
    ch_versions = ch_versions.mix(TXIMETA_TXIMPORT.out.versions)


    emit:
    results                       = ch_pseudo_results                              // channel: [ val(meta), results_dir ]
    multiqc                       = ch_pseudo_multiqc                              // channel: [ val(meta), files_for_multiqc ]

    tpm_gene                      = TXIMETA_TXIMPORT.out.tpm_gene                  //    path: *gene_tpm.tsv
    counts_gene                   = TXIMETA_TXIMPORT.out.counts_gene               //    path: *gene_counts.tsv
    lengths_gene                  = TXIMETA_TXIMPORT.out.lengths_gene              //    path: *gene_lengths.tsv
    counts_gene_length_scaled     = TXIMETA_TXIMPORT.out.counts_gene_length_scaled //    path: *gene_counts_length_scaled.tsv
    counts_gene_scaled            = TXIMETA_TXIMPORT.out.counts_gene_scaled        //    path: *gene_counts_scaled.tsv
    tpm_transcript                = TXIMETA_TXIMPORT.out.tpm_transcript            //    path: *gene_tpm.tsv
    counts_transcript             = TXIMETA_TXIMPORT.out.counts_transcript         //    path: *transcript_counts.tsv
    lengths_transcript            = TXIMETA_TXIMPORT.out.lengths_transcript        //    path: *transcript_lengths.tsv

    versions                      = ch_versions                                    // channel: [ versions.yml ]
}






