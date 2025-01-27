/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { INPUT_CHECK                   } from '../subworkflows/local/input_check' // validate the input samplesheet.csv and prepare input channels
include { OBTAIN_TRANSCRIPTOME          } from '../subworkflows/local/obtain_transcriptome' // prepare transcriptome
include { RUN_PSEUDO_ALIGNMENT          } from '../subworkflows/local/run_pseudo_alignment.nf'
include { RUN_STAR_SALMON	            } from '../subworkflows/local/run_star_salmon.nf'
include { RUN_STAR_RSEM            } from '../subworkflows/local/run_star_rsem.nf'
include { DGE_ANALYSIS_DESEQ2           } from '../modules/local/dge_analysis_deseq2.nf'
include { DGE_ANALYSIS_EDGER            } from '../modules/local/dge_analysis_edger.nf'


include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ              } from '../modules/nf-core/cat/fastq/main'


include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_rnaseqdge_pipeline'





/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQDGE {

    take:
    input_samplesheet	//  file: path to input samplesheet
    aligner				//  string: aligner method
    genome_fasta
    genome_gtf
    ensembl_release
    

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()



    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        input_samplesheet
    )
    .reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)


    //
    // MODULE: concatenate fastq files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))


    //
    // MODULE: run FastQC
    //
    FASTQC (
        ch_cat_fastq
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    //
    // SUBWORKFLOW: obtain genome/transcriptome
    //
    OBTAIN_TRANSCRIPTOME (
    	aligner,
		genome_fasta,
		genome_gtf,
		ensembl_release
    )
    ch_versions = ch_versions.mix(OBTAIN_TRANSCRIPTOME.out.versions.first())


    //
    // SUBWORKFLOW: mapping step
    //
	if (aligner == "star_rsem") {
		ch_gene_counts = RUN_STAR_RSEM(ch_cat_fastq, aligner, OBTAIN_TRANSCRIPTOME.out.fasta, OBTAIN_TRANSCRIPTOME.out.gtf).counts_gene
	} else if (aligner == "star_salmon") {
		ch_gene_counts = RUN_STAR_SALMON(ch_cat_fastq, aligner, OBTAIN_TRANSCRIPTOME.out.fasta, OBTAIN_TRANSCRIPTOME.out.gtf).counts_gene
	} else if (aligner in ["salmon","kallisto"]) {
		ch_gene_counts = RUN_PSEUDO_ALIGNMENT(ch_cat_fastq, aligner, OBTAIN_TRANSCRIPTOME.out.fasta, OBTAIN_TRANSCRIPTOME.out.gtf).counts_gene
	} else {
		println "no matching aligner found"
	}


    //
    // MODULE: differential gene expression analysis
    //
	DGE_ANALYSIS_DESEQ2(ch_gene_counts, input_samplesheet)
	DGE_ANALYSIS_EDGER(ch_gene_counts, input_samplesheet)


/*
    //
    // TODO: Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }
*/
    
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
    //ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)  ## TODO uncomment and fix versions
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
        ch_multiqc_logo.toList()
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
