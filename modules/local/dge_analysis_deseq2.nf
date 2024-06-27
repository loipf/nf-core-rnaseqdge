process DGE_ANALYSIS_DESEQ2 {
    tag "$gene_counts"
    label "process_medium"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-sartools:1.8.1--r43hdfd78af_2' :
        'quay.io/biocontainers/r-sartools:1.8.1--r43hdfd78af_2' }"

    input:
    tuple val(meta), path(gene_counts)
    path(sample_sheet)

    output:
    path "*"
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when


    beforeScript 'apt-get update && apt-get install -y libgmp10 && apt-get clean'


	script: 
    """
    
    ### edit for sartools
    mv $gene_counts "all_samples.gene_counts.tsv"
    cut -d',' -f1,2,4- "$sample_sheet" | tr ',' '\t' > sample_sheet.tsv
    
    rnaseqdge_sartools_script_DESeq2.R
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(echo \$(R --version | head -n 1) | sed -e "s/R version //g")
    END_VERSIONS
    """
}
