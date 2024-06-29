FROM pegi3s/r_data-analysis:latest

MAINTAINER loipf


RUN apt-get --allow-releaseinfo-change update && apt install -y libcurl4-openssl-dev libssl-dev gzip curl libxml2-dev pandoc

RUN R -e 'install.packages("pacman",repos = "http://cran.us.r-project.org")'
RUN R -e 'pacman::p_load(BiocManager, plotly, pheatmap, vsn, hexbin, LSD, cowplot, DESeq2, edgeR)'




# ### docker image for DGE_ANALYSIS_DESEQ2 and DGE_ANALYSIS_EDGER modules
# docker tag nfcore_dge_analysis:latest loipf/nf_rnaseqdge_dge_analysis:latest
# docker push loipf/nf_rnaseqdge_dge_analysis:latest



