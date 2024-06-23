#!/usr/bin/env Rscript

################################################################################
### R script to compare several conditions with the SARTools and edgeR packages
### Hugo Varet
### March 23rd, 2022
### designed to be executed with SARTools 1.8.1
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
rm(list=ls())                                        # remove all the objects from the R session

workDir <- "."      # working directory for the R session
# workDir <- "/home/stefanloipfinger/Desktop/startools/"      # working directory for the R session

projectName <- "rnaseqdge_startools_edger"                         # name of the project
author <- "loipf"                                # author of the statistical analysis/report

targetFile <- "sample_sheet.tsv"                     # path to the design/target file
rawDir <- "/here"                                      # path to the directory containing raw counts files
featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

counts_file_path = "all_samples.gene_counts.tsv"  # gene counts file: genes x samples


varInt <- "group"                                    # factor of interest
condRef <- "control"                                   # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

alpha <- 0.05                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

colors <- c("#f3c300", "#875692", "#f38400",         # vector of colors of each biological condition on the plots
            "#a1caf1", "#be0032", "#c2b280",
            "#848482", "#008856", "#e68fac",
            "#0067a5", "#f99379", "#604e97")

forceCairoGraph <- FALSE

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library(SARTools)
if (forceCairoGraph) options(bitmapType="cairo")

# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# # loading counts
# counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)
# counts = data.frame(data.table::fread(counts_file_path), row.names = 1)
counts = data.frame(read.table(counts_file_path, sep = "\t", quote = "\"",header=T, stringsAsFactors = FALSE), row.names=1)
counts$gene_name = NULL
counts = counts[,target$sample]
counts = round(counts)


dir.create("figures", showWarnings = F)
# description plots - fails all the time - do manual
# majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)
counts=counts; group=target[,varInt]; col=colors;ggplot_theme = theme_light()
barplotTotal(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
barplotNull(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
densityPlot(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
majSequences <- majSequences(counts = counts, group = group, col = col, ggplot_theme = ggplot_theme)
# cat("Matrix of SERE statistics:\n")
# print(tabSERE(counts))
# pairwiseScatterPlots(counts = counts, ggplot_theme = ggplot_theme)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
tryCatch( {

### docker quay.io/biocontainers/r-sartools:1.8.1--r43hdfd78af_2 missing library:
# /usr/local/bin/pandoc: error while loading shared libraries: libgmp.so.10: cannot open shared object file: No such file or directory
# Error in system(paste(shQuote(path), "--version"), intern = TRUE) : 
# 	error in running command
# Calls: writeReport.DESeq2 ... get_pandoc_version -> with_pandoc_safe_environment -> force -> system

	writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
		              majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
		              targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
		              condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, cpmCutoff=cpmCutoff,
		              colors=colors, gene.selection=gene.selection, normalizationMethod=normalizationMethod)
    },
    error = function(cond) {
        message(conditionMessage(cond))
    }
)

