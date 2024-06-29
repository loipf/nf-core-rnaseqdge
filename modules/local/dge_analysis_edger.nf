process DGE_ANALYSIS_EDGER {
    tag "$gene_counts"
    label "process_medium"

    container 'docker.io/loipf/nf_rnaseqdge_dge_analysis'

    input:
    tuple val(meta), path(gene_counts)
    path(sample_sheet)

    output:
    path "dge_edger_analysis.html"
    path "dge_edger_results.tsv"
    path "gene_counts_tmm_cpm.tsv"

    when:
    task.ext.when == null || task.ext.when


    script: 
    """
    
    Rscript -e "
    curr_dir = getwd();
    rmarkdown_file_path = system('which run_dge_edger.R', intern = TRUE);
	params_list = list(curr_dir = curr_dir, gene_count_matrix_path='$gene_counts', sample_anno_path='$sample_sheet');
	print(params_list);
	rmarkdown::render(rmarkdown_file_path, output_file=file.path(curr_dir,'dge_edger_analysis.html'), params = params_list );
	"

    """
    
}
