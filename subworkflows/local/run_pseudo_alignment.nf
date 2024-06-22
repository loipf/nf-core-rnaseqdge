//
// run pseudo alignment pipelines
//


include { QUANTIFY_PSEUDO_ALIGNMENT } from '../../subworkflows/nf-core/quantify_pseudo_alignment/main' 
include { QUANTIFY_PSEUDO_ALIGNMENT } from '../../subworkflows/nf-core/quantify_pseudo_alignment/main' 





workflow RUN_PSEUDO_ALIGNMENT {
    
    take:
    ch_fastq  // channel: 
    aligner		 //  string: aligner method
    fasta                //      file: /path/to/genome.fasta
    gtf                  //      file: /path/to/genome.gtf
    

    main:

    ch_versions = Channel.empty()







    //
    // download genome/transcriptome from ensembl if not given
    //
    if(fasta==null) {
		gtf = DOWNLOAD_ENSEMBL_GTF(ensembl_release).gtf
		
		if (aligner in ["star_rsem", "star_salmon"]) {
			fasta = DOWNLOAD_ENSEMBL_GENOME(ensembl_release).fasta
		} else if (aligner in ["salmon","kallisto"]) {
			fasta = DOWNLOAD_ENSEMBL_TRANSCRIPTOME(ensembl_release).fasta
		} else {
			println "no matching aligner found"
		}
    }
    

    //
    // Uncompress genome fasta file if required
    //
    if (fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            fasta.map{ it -> [[id:it[0].baseName], it] }
        )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map{ meta, fasta -> [fasta] }.collect()
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = fasta.collect()
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gffread_version = Channel.empty()
    if (gtf.endsWith('.gz')) {
        GUNZIP_GTF (
            gtf.map{ it -> [[id:it[0].baseName], it] }
        )
        ch_gtf = GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect()
        ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
    } else {
        ch_gtf = gtf.collect()
    }



    emit:
    fasta            = ch_fasta            // path: genome.fasta
    gtf              = ch_gtf              // path: genome.gtf
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}






def getIndexVersion( index_path ) {
    genomeParameters = new File("$index_path/genomeParameters.txt")
    if ( genomeParameters.exists() ) {
        for(line: genomeParameters.readLines()){
            if(line.startsWith("versionGenome")){
                return line.split("\t")[1].trim()
            }
        }
    }
}







/*

include { GATK4_CREATESEQUENCEDICTIONARY }    from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main' //addParams(options: params.genome_options)
include { GFFREAD }                           from '../../modules/nf-core/modules/gffread/main'                        //addParams(options: params.gffread_options)
include { GTF2BED }                           from '../../modules/local/gtf2bed'                                       //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_FASTA }            from '../../modules/nf-core/modules/gunzip/main'                         //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GENE_BED }         from '../../modules/nf-core/modules/gunzip/main'                         //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GFF }              from '../../modules/nf-core/modules/gunzip/main'                         //addParams(options: params.genome_options)
include { GUNZIP as GUNZIP_GTF }              from '../../modules/nf-core/modules/gunzip/main'                         //addParams(options: params.genome_options)
include { SAMTOOLS_FAIDX }                    from '../../modules/nf-core/modules/samtools/faidx/main'                 //addParams(options: params.genome_options)
include { STAR_GENOMEGENERATE }               from '../../modules/nf-core/modules/star/genomegenerate/main'            //addParams(options: params.star_index_options)
include { UNTAR as UNTAR_STAR_INDEX }         from '../../modules/nf-core/modules/untar/main'                          //addParams(options: params.star_untar_options)


workflow PREPARE_GENOME {
    take:
    prepare_tool_indices
    
    additional_fasta     //      file: /path/to/additional.fasta
    transcript_fasta     //      file: /path/to/transcript.fasta
    gene_bed             //      file: /path/to/gene.bed
    splicesites          //      file: /path/to/splicesites.txt
    bbsplit_fasta_list   //      file: /path/to/bbsplit_fasta_list.txt
    star_index           // directory: /path/to/star/index/
    rsem_index           // directory: /path/to/rsem/index/
    salmon_index         // directory: /path/to/salmon/index/
    kallisto_index       // directory: /path/to/kallisto/index/
    hisat2_index         // directory: /path/to/hisat2/index/
    bbsplit_index        // directory: /path/to/rsem/index/
    gencode              //   boolean: whether the genome is from GENCODE
    is_aws_igenome       //   boolean: whether the genome files are from AWS iGenomes
    biotype              //    string: if additional fasta file is provided biotype value to use when appending entries to GTF file
    prepare_tool_indices //      list: tools to prepare indices for
    filter_gtf           //   boolean: whether to filter GTF file

    
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        GUNZIP_FASTA (
            Channel.fromPath(params.fasta).map{ it -> [[id:it[0].baseName], it] }
        )
        ch_fasta = GUNZIP_FASTA.out.gunzip.map{ meta, fasta -> [fasta] }.collect()
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.fromPath(params.fasta).collect()
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    ch_gffread_version = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            GUNZIP_GTF (
                Channel.fromPath(params.gtf).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gtf = GUNZIP_GTF.out.gunzip.map{ meta, gtf -> [gtf] }.collect()
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.fromPath(params.gtf).collect()
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            GUNZIP_GFF (
                Channel.fromPath(params.gff).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gff = GUNZIP_GFF.out.gunzip.map{ meta, gff -> [gff] }.collect()
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.fromPath(params.gff).collect()
        }

        GFFREAD (
            ch_gff
        )
        .gtf
        .set { ch_gtf }

        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress exon BED annotation file or create from GTF if required
    //
    if (params.exon_bed) {
        if (params.exon_bed.endsWith('.gz')) {
            GUNZIP_GENE_BED (
                Channel.fromPath(params.exon_bed).map{ it -> [[id:it[0].baseName], it] }
            )
            ch_gene_bed = GUNZIP_GENE_BED.out.gunzip.map{ meta, bed -> [bed] }.collect()
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.fromPath(params.exon_bed).collect()
        }
    } else {
        ch_exon_bed = GTF2BED ( ch_gtf ).bed.collect()
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    // Index the genome fasta
    ch_fasta_fai = Channel.empty()
    if (params.fasta_fai) ch_fasta_fai = Channel.fromPath(params.fasta_fai).collect()
    if (!params.fasta_fai) {
        SAMTOOLS_FAIDX(
            ch_fasta.map{ it -> [[id:it[0].getName()], it]}
        )
        ch_fasta_fai    = SAMTOOLS_FAIDX.out.fai.map{ meta, fai -> [fai] }.collect()
        ch_versions     = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    }

    // Create dictionary file for the genome fasta
    ch_fasta_dict = Channel.empty()
    if (params.dict) ch_fasta_dict = Channel.fromPath(params.dict).collect()
    else ch_fasta_dict = GATK4_CREATESEQUENCEDICTIONARY(ch_fasta).dict

    //
    // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                UNTAR_STAR_INDEX (
                    Channel.fromPath(params.star_index).map{ it -> [[id:it[0].baseName], it] }
                )
                ch_star_index = UNTAR_STAR_INDEX.out.untar.map{ meta, star_index -> [star_index] }.collect()
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.fromPath(params.star_index).collect()
            }
        }
        else {
            STAR_GENOMEGENERATE (
                ch_fasta,ch_gtf
            )
            .index
            .set { ch_star_index }
            ch_versions     = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }

        //if((!ch_star_index) || getIndexVersion(ch_star_index) != '2.7.4a'){
        //    ch_star_index   = STAR_GENOMEGENERATE(ch_fasta,ch_gtf).index
        //    ch_versions     = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        //}
    }


    emit:
    fasta            = ch_fasta            // path: genome.fasta
    fai              = ch_fasta_fai        // path: genome.fasta.fai
    dict             = ch_fasta_dict       // path: genome.fasta.dict
    gtf              = ch_gtf              // path: genome.gtf
    exon_bed         = ch_exon_bed         // path: exon.bed
    star_index       = ch_star_index       // path: star/index/
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}

def getIndexVersion( index_path ) {
    genomeParameters = new File("$index_path/genomeParameters.txt")
    if ( genomeParameters.exists() ) {
        for(line: genomeParameters.readLines()){
            if(line.startsWith("versionGenome")){
                return line.split("\t")[1].trim()
            }
        }
    }
}


*/
