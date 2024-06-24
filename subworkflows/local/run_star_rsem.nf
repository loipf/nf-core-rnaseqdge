//
// run STAR alignment pipelines
//


include { FASTQ_ALIGN_STAR } from '../../subworkflows/nf-core/fastq_align_star/main'
include { STAR_GENOMEGENERATE } from '../../modules/nf-core/star/genomegenerate'


include { RSEM_PREPAREREFERENCE } from '../../modules/nf-core/rsem/preparereference/main'
include { RSEM_CALCULATEEXPRESSION } from '../../modules/nf-core/rsem/calculateexpression/main'


include { QUANTIFY_RSEM } from '../../subworkflows/local/quantify_rsem'



include { SALMON_INDEX } from '../../modules/local/salmon_index_custom'
include { SALMON_QUANT } from '../../modules/nf-core/salmon/quant'
include { CUSTOM_TX2GENE   } from '../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../modules/nf-core/tximeta/tximport'




workflow RUN_STAR_RSEM {
    
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
	RSEM_PREPAREREFERENCE(fasta.map{ fasta -> [fasta[1]] }.collect(), gtf)
	
	
	//
    // align reads and run quantification
    //
	QUANTIFY_RSEM(reads, RSEM_PREPAREREFERENCE.out.index, fasta)
	ch_versions = ch_versions.mix(QUANTIFY_RSEM.out.versions)
    //ch_pseudo_multiqc = QUANTIFY_RSEM.out.stat


    emit:
    counts_gene						= QUANTIFY_RSEM.out.merged_counts_gene.collect{ it }.map { [ ['id': 'all_samples'], it ] }
    tpm_gene						= QUANTIFY_RSEM.out.merged_tpm_gene.collect{ it }.map { [ ['id': 'all_samples'], it ] }
    counts_transcript				= QUANTIFY_RSEM.out.merged_counts_transcript.collect{ it }.map { [ ['id': 'all_samples'], it ] }
    tpm_transcript					= QUANTIFY_RSEM.out.merged_tpm_transcript.collect{ it }.map { [ ['id': 'all_samples'], it ] }
    
    versions                      = ch_versions                                    // channel: [ versions.yml ]

}





