process DOWNLOAD_ENSEMBL_TRANSCRIPTOME {
    tag "$ensembl_release"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"

    input:
    val ensembl_release

    output:
    path '*.fa.gz'       , emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script: 
    """
    
    wget -O Homo_sapiens.GRCh38.cdna.fa.gz ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz

    """
}
