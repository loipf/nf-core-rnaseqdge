#!/usr/bin/env Rscript


############################################
### common parameter and functions

#install.packages("pacman")
pacman::p_load(knitr, data.table, plotly, ggplot2, pheatmap, vsn, hexbin, LSD, cowplot, DESeq2, edgeR)



############################################
### general functions


load_dataframe = function(file_path, ... ) {
  df = data.frame(fread(file_path, header=T, stringsAsFactors = F), row.names = 1, check.names = F, ...)
  return(df)
}

save_dataframe = function(df, file_path, ...) {
  if(!is.null(file_path)) {
    file_path_ext = tools::file_ext(file_path)
    default_sep = switch(file_path_ext,
                         "csv"=",",
                         "tsv"="\t",
                         "txt"="\t",
                         ",")
    if(file_path_ext=="") {  file_path = paste0(file_path,".tsv")  }  ### default
    
    
    ### numeric rownames sometimes still cause problems in fwrite
    ### see: https://github.com/Rdatatable/data.table/issues/4957
    suppressWarnings({ numeric_rownames = as.numeric(rownames(df)) })
    if(is.numeric(numeric_rownames) & !all(is.na(numeric_rownames)) ) {
      warning('numeric rownames: rownames in output file may be incorrect')
      
      # ### alternative
      # write.table(as.data.frame(df), file_path, row.names = T, col.names = NA, sep = default_sep, quote = F, ...)
    }
    
    fwrite(as.data.frame(df), file_path, row.names=T, sep=default_sep, ...)
  } else {
    warning("file_path not given, file will not be saved")
  }
}

### print with timestamp
print_time = function(txt){
  print(paste0("### ",Sys.time(),":   ",txt))
}


############################################
### read in parser functions


sa_find_group_column = function(sa, group_name=NA){
  ### rownames: sample, first 2 colnames: fastq_1, fastq_2
  available_colnames = grep("fastq", colnames(sa), invert=T, value = T, ignore.case = T)
  if(length(available_colnames) == 0) {
    stop("no group column found")
  } 
  if(!is.na(group_name) & group_name %in% available_colnames){
      return(group_name)
  } else {
      group_column = grep("(group|treatment|condition)", colnames(sa), value = T, ignore.case = T)
      print_time(paste0("automatic group column detection: ",group_column))
      return(group_column)
  }
}

sa_group_set_reference = function(sa_group, reference_group=NA){
  if(length(unique(sa_group))!=2) {
    stop("more than 2 groups found")
  }
  group_levels = unique(sa_group)
  
  if(is.na(reference_group)) {
    reference_group = grep("(baseline|control|normal|untreated)", group_levels, value = T, ignore.case = T)
    if(length(reference_group) == 0){  ### no reference_group label found
      reference_group = group_levels[0]
    }
  } 
  print_time(paste0("automatic reference group detection: ",reference_group))
  group_levels = c(reference_group, setdiff(group_levels, reference_group))
  sa_group = factor(sa_group, levels = group_levels)
  return(sa_group)
}

### preprocess sample_anno_df
read_in_sample_anno = function(sample_anno_path, group_name=NA, reference_group=NA) {
  sample_anno_df = load_dataframe(sample_anno_path)
  group_column = sa_find_group_column(sample_anno_df, group_name)
  sample_anno_df = sample_anno_df[!is.na(sample_anno_df[[group_column]]) | sample_anno_df[[group_column]] %in% c("","not_available"),]
  sample_anno_df[["group"]] = sa_group_set_reference(sample_anno_df[[group_column]], reference_group)
  return(sample_anno_df)
}


### preprocess counts data
read_in_gene_df = function(gene_count_matrix_path) {
  gene_counts_df = load_dataframe(gene_count_matrix_path)
  gene_counts_df = gene_counts_df[,colnames(gene_counts_df) != "gene_name"]
  print_time(paste0("raw gene counts (samples x genes): ", ncol(gene_counts_df)," x ", nrow(gene_counts_df)))
  
  ### filter low expressed genes: 10 counts in 1% of the samples
  gene_counts_df = gene_counts_df[rowSums(gene_counts_df>10) > ncol(gene_counts_df)*0.01,]
  print_time(paste0("filtered gene counts (samples x genes): ", ncol(gene_counts_df)," x ", nrow(gene_counts_df)))
  return(gene_counts_df)
}




############################################
### plot functions

center_rowwise <- function(x, na.rm=T) {
  x_center = x - rowMeans(x, na.rm=na.rm)
  return(x_center)
}

plot_heatscatter <- function(x, y, ...) {
  plot_df = data.frame(x=as.vector(x), y=as.vector(y))
  plot_df = na.omit(plot_df)
  LSD::heatscatter(plot_df$x, plot_df$y, cor=TRUE, colpal="standard", ...)
  abline(0,1, col="orange")
}

plot_density_genes <- function(gene_matrix, main="", xlab="values") {
  plot_df = stack(gene_matrix[seq(1, nrow(gene_matrix), 10),])   ### subset to avoid overload
  colnames(plot_df)[2] = "sample_id"
  p <- ggplot(plot_df, aes(x=values)) + labs(x=xlab, title=main) +
    stat_density(aes(color=sample_id),position="identity",geom="line") +
    theme_bw() + theme(legend.position='none')
  p
}

plot_density_genes_cumulative <- function(gene_matrix, main="", xlab="values") {
  gene_matrix = gene_matrix[seq(1, nrow(gene_matrix), 10),]   ### subset to avoid overload
  plot_df_list = sapply(colnames(gene_matrix), function(sample_id) {
    x_value = sort(gene_matrix[,sample_id])
    cumul_value = cumsum(x_value)/sum(x_value)
    data.frame(x=x_value, y=cumul_value, sample_id=sample_id)
  }, simplify = F, USE.NAMES = F)
  plot_df = do.call("rbind", plot_df_list)
  
  p = ggplot(plot_df, aes(x=x, y=y, color=sample_id)) + geom_line() +
    theme_bw() + labs(y="density", x=xlab, title=main) + theme(legend.position='none')
  p
}



