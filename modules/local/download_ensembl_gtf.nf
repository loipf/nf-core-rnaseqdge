process DOWNLOAD_ENSEMBL_GTF {
    tag "$ensembl_release"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val ensembl_release

    output:
    path '*.gtf.gz'       , emit: gtf

    when:
    task.ext.when == null || task.ext.when

    script: 
    """
    
    wget -O Homo_sapiens.GRCh38.gtf.gz ftp://ftp.ensembl.org/pub/release-$ensembl_release/gtf/homo_sapiens/Homo_sapiens.GRCh38.112.chr_patch_hapl_scaff.gtf.gz

    """
}
