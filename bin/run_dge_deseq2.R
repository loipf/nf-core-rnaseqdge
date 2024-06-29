#' ---
#' title: "differential gene expression analysis using DESeq2"
#' author: "loipf"
#'
#' output: html_document
#' code_folding: hide
#' 
#' params:
#'    curr_dir: NULL
#'    gene_count_matrix_path: NULL
#'    sample_anno_path: NULL
#' ---

source("dge_shared_functions.R")


#' ***
#' # read in

print(params)

GENE_COUNT_MATRIX_PATH = file.path(params["curr_dir"], params["gene_count_matrix_path"])
SAMPLE_ANNO_PATH = file.path(params["curr_dir"], params["sample_anno_path"])

sample_anno_df = read_in_sample_anno(SAMPLE_ANNO_PATH)
gene_counts_df = read_in_gene_df(GENE_COUNT_MATRIX_PATH)

common_samples = intersect(rownames(sample_anno_df), colnames(gene_counts_df))
sample_anno_df = sample_anno_df[common_samples,]
gene_counts_df = gene_counts_df[,common_samples]


#' ***
#' # dge analysis

dds = DESeqDataSetFromMatrix(countData = round(gene_counts_df), colData = sample_anno_df, design = ~ group)
dds = DESeq(dds, fitType = "parametric", sfType="ratio")
res = results(dds)

print(summary(res))

print(res[order(res[["pvalue"]]),][1:10,])


### get different normalized gene expression dfs
gene_counts_vst_dds = varianceStabilizingTransformation(dds, blind=F)
gene_counts_vst_df = data.frame(assay(gene_counts_vst_dds), check.names = F)
gene_counts_log2_df = data.frame(log2(gene_counts_df+1),check.names = F)


### export
save_dataframe(res, file.path(params["curr_dir"], "dge_deseq2_results.tsv"))
save_dataframe(gene_counts_vst_df, file.path(params["curr_dir"], "gene_counts_sf_vst.tsv"))



#' ***
#' # QC plots

#' ***
sf_df = data.frame("size_factors"=sizeFactors(dds))
sf_df[["sample_id"]] = rownames(sf_df)

#' ### DESeq2 size factors (should be around 0.7-1.3)
p = ggplot(sf_df, aes(x = size_factors, label=sample_id)) +
  geom_rug(sides="b") + geom_line(stat="density") +
  theme_bw() + labs(x="size factors", title="DESeq2 size factor estimation") +
  geom_vline(aes(xintercept = 1),col='orange',linewidth=0.5)
ggplotly(p)

#' ***

#' ### sample correlation
COLOR_PALETTE <- c("#2166AC","#67A9CF","#D1E5F0","#F7F7F7","#FDDBC7","#EF8A62","#B2182B") ## rev(brewer.pal(n = 7, name = "RdBu"))
gene_counts_vst_df_center = center_rowwise(gene_counts_vst_df)
corr_matrix = cor(gene_counts_vst_df_center)
#+ fig.width=10, fig.height=9
pheatmap(corr_matrix, main="sample correlation gene_counts",
         breaks = seq(-1,1, length.out = 100), na_col = "grey", color=colorRampPalette(COLOR_PALETTE)(100) )


#' ***

#' ### gene count distribution
#+ fig.width=10, fig.height=7
mean_gene_counts_log2_df = mean(apply(gene_counts_log2_df,2,median))
boxplot(gene_counts_log2_df, main="gene count distribution", ylab="log2(expression+1)", las=2)
abline(mean_gene_counts_log2_df, 0, col="orange", lty=1)

mean_gene_counts_vst_df = mean(apply(gene_counts_vst_df,2,median))
boxplot(gene_counts_vst_df, main="gene count distribution [sf+vst]", ylab="normalized expression", las=2)
abline(mean_gene_counts_vst_df, 0, col="orange", lty=1)

#' ***

#' ### compare first 2 samples
#+ fig.width=7, fig.height=7
plot_heatscatter(gene_counts_log2_df[,1], gene_counts_log2_df[,2], main="gene count sample comparison [log2(expression+1)]",
                 xlab=colnames(gene_counts_log2_df)[1], ylab=colnames(gene_counts_log2_df)[2])
plot_heatscatter(gene_counts_vst_df[,1], gene_counts_vst_df[,2], main="gene count sample comparison [sf + vst]",
                 xlab=colnames(gene_counts_log2_df)[1], ylab=colnames(gene_counts_log2_df)[2])

#' ***

#' ### raw to vst heteroscedasticity
non_zero_count_genes = names(!rowSums(gene_counts_df)==0)

#+ fig.width=10, fig.height=5
p1 = vsn::meanSdPlot(as.matrix(gene_counts_df[non_zero_count_genes,]), plot=F)$gg + ggtitle("raw") + theme_bw()
p2 = vsn::meanSdPlot(as.matrix(gene_counts_log2_df[non_zero_count_genes,]), plot=F)$gg + ggtitle("log2 normalized") + theme_bw()
p3 = vsn::meanSdPlot(as.matrix(gene_counts_vst_df[non_zero_count_genes,]), plot=F)$gg + ggtitle("sf+vst normalized") + theme_bw()
p_grid = plot_grid(p1,p2,p3, labels="AUTO", ncol=3)
title <- ggdraw() + draw_label("gene counts mean-sd", fontface='bold')
plot_grid(title, p_grid, ncol=1, rel_heights=c(0.1, 1))

#' ***

#' ### gene count density plots
ggplotly(plot_density_genes(gene_counts_vst_df, main="gene count densities [sf+vst]", xlab="sf+vst expression"))

#' ***

#' ### gene count density plots cumulative
ggplotly(plot_density_genes_cumulative(gene_counts_vst_df, main="gene count densities [sf+vst]", xlab="sf+vst expression"))


#' ***
#' ### dge results

plotDispEsts(dds)
plotMA(res, main="MA-plot", ylim=c(-2,2))
hist(res$pvalue, breaks=50, main="pvalue distribution", xlab="pvalue")


#' ### volcano plot
res_df = as.data.frame(res)
res_df$significance = ifelse(res_df$padj < 0.05, "significant", "not significant")
ggplot(res_df, aes(x=log2FoldChange, y=-log10(pvalue), color=significance)) +
  geom_point(alpha=0.4, pch=15) +
  scale_color_manual(values = c("grey", "red")) +
  theme_minimal() +
  ggtitle("volcano plot") +
  xlab("log2 fold change") +
  ylab("-log10 p-value")

#' ### heatmap of top differentially expressed genes
top_genes = head(order(res$padj), 20)
plot_top_genes_df = gene_counts_vst_df[top_genes, ]
plot_top_genes_df = center_rowwise(plot_top_genes_df)
pheatmap(plot_top_genes_df, annotation_col=sample_anno_df[,"group",drop=F], main="top 20 DEGs")

#' ### PCA plot
pca_data = plotPCA(gene_counts_vst_dds, intgroup="group", returnData=TRUE)
percent_var = round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(PC1, PC2, color=group)) +
  geom_point(size=2) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  ggtitle("PCA of vst data") +
  theme_minimal()


#' ***
#' ### session info
print(sessionInfo())


