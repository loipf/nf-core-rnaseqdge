process DOWNLOAD_ENSEMBL_GENOME {
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
    
    wget -O Homo_sapiens.GRCh38.primary_assembly.fa.gz  ftp://ftp.ensembl.org/pub/release-$ensembl_release/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    
    """
}
