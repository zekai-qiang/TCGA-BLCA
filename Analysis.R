setwd("/Volumes/MEDICINE/TCGA_BLCA")

if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TCGAbiolinks")
BiocManager::install("remotes")
BiocManager::install("Rhtslib")
BiocManager::install("DNAcopy")
BiocManager::install("maftools")
BiocManager::install("sparseMatrixStats")
BiocManager::install("DelayedMatrixStats")
BiocManager::install("edgeR")

source("https://mac.R-project.org/bin/install.R")

install.packages("ggcorrplot")
install.packages("microbenchmark")
install.packages("rbenchmark")
install.packages("coin")
install.packages("fitdistrplus")
install.packages("boot")

library(TCGAbiolinks)
library(dplyr)
library(DT)
library(remotes)
library(Rhtslib)
library(DNAcopy)
library(maftools)
library(SummarizedExperiment)
library(readr)
library(ggplot2)
library(ggpubr)
library(ggcorrplot)
library(tidyverse)
library(viridis)
library(coin)
library(pastecs)
library(rbenchmark)
library(edgeR)
library(fitdistrplus)
library(boot)
library(magrittr)
library(tidyverse)

packageVersion("TCGAbiolinks")

TCGAbiolinks:::getGDCprojects()$project_id
TCGAbiolinks:::getProjectSummary("TCGA-BLCA")

query.blca_exp <- GDCquery(
  project = "TCGA-BLCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
blca_exp <- GDCprepare(query.blca_exp)
blca_exp

assays(blca_exp)
rownames(blca_exp)
rowRanges(blca_exp)
colnames(blca_exp)
colData(blca_exp)
metadata(blca_exp)

assay <- assay(blca_exp) %>% as.data.frame()
rowranges <- rowRanges(blca_exp) %>% as.data.frame()
coldata <- colData(blca_exp) %>% as.data.frame()

null.fc_filter_default <- function(list) {
  fltrd_df.deg <- lapply(1:length(list), 
                         function(f) list[[f]][list[[f]]$fold.change > 1 | 
                                                 list[[f]]$fold.change < -1,])
  dim <- lapply(1:length(fltrd_df.deg), function(e) dim(fltrd_df.deg[[e]]))
  effect <- which(lapply(1:length(fltrd_df.deg), function(d) dim[[d]][1] != 0) == TRUE)
  fltrd_df.deg <- fltrd_df.deg[effect]
  return(fltrd_df.deg)
} #filters dim 0 and fc > 1, < -1

null.fc_filter <- function(list, ufc, lfc) {
  fltrd_df.deg <- lapply(1:length(list), 
                         function(f) list[[f]][list[[f]]$fold.change > ufc | 
                                                 list[[f]]$fold.change < lfc,])
  dim <- lapply(1:length(fltrd_df.deg), function(e) dim(fltrd_df.deg[[e]]))
  effect <- which(lapply(1:length(fltrd_df.deg), function(d) dim[[d]][1] != 0) == TRUE)
  fltrd_df.deg <- fltrd_df.deg[effect]
  return(fltrd_df.deg)
} #filters dim 0 and fc > ufc, < lfc

incl.coldata <- coldata[!is.na(coldata$paper_patient),] %>%
  filter(sample_type == "Primary Tumor")
incl.barcode <- incl.coldata$barcode #412
incl.assay <- assay[,incl.barcode]
hist(log2(as.numeric(as.matrix(incl.assay))), freq = FALSE)

genes.paper_mutation <- as.list(colnames(incl.coldata %>% 
                                   select(starts_with("paper_mutation.in."))))
str(genes.paper_mutation)
genes.paper_mutation[[1]]

paper_mutation.in.bar <- function(m) {
  mut.bar <- incl.coldata$barcode[incl.coldata[,m] == "yes"]
  nmut.bar <- incl.coldata$barcode[incl.coldata[,m] == "no"]
  group.assay <- cbind(incl.assay[,mut.bar], incl.assay[,nmut.bar])
  dgelist.group <- c(rep("mut", length(mut.bar)), rep("nmut", length(nmut.bar)))
  dgelist <- DGEList(counts = group.assay, group = dgelist.group)
  fltrd.dgelist <- dgelist[filterByExpr(dgelist), ,keep.lib.sizes = FALSE]
  norm.dgelist <- as.data.frame(cpm(calcNormFactors(fltrd.dgelist, method = "TMM")))
  fltrd.gene <- rownames(norm.dgelist)
  mut.assay <- norm.dgelist[,mut.bar]
  nmut.assay <- norm.dgelist[,nmut.bar]
  adj.p_fun <- sapply(1:length(fltrd.gene), 
                      function(g) wilcox.test(as.numeric(mut.assay[g,]),
                                              as.numeric(nmut.assay[g,]),
                                              alternative = "two.sided",
                                              exact = FALSE)$p.value) %>%
    {p.adjust(., method = "BH")}
  deg <- fltrd.gene[adj.p_fun < 0.05]
  df.deg <- data.frame(mut.mean = rowMeans(norm.dgelist[deg, mut.bar]),
                       nmut.mean = rowMeans(norm.dgelist[deg, nmut.bar]))
  df.deg$fold.change = df.deg %>% {log2(.$mut.mean/.$nmut.mean)}
  return(df.deg)
}

paper_mutation.in.bar(genes.paper_mutation[[1]]) #test

list_df.deg <- lapply(1:length(genes.paper_mutation), 
                      function(l) paper_mutation.in.bar(genes.paper_mutation[[l]]))
for(i in 1:length(list_df.deg)) {
  attr(list_df.deg[[i]], "paper_mutation") <- genes.paper_mutation[[i]] 
}
str(list_df.deg)

fltrd.list <- null.fc_filter_default(list_df.deg)
str(fltrd.list)

#quartiles apobec
DEG.APOBEC_load <- function(apobec) {
  apobec_quartiles <- quantile(apobec, prob = c(.25, .5, .75))
  apobec.high <- incl.coldata[apobec > apobec_quartiles[[3]], "barcode"]
  apobec.low <- incl.coldata[apobec < apobec_quartiles[[1]], "barcode"]
  group.apobec <- cbind(incl.assay[,apobec.high], incl.assay[,apobec.low])
  apobec.dgegroup <- c(rep("high", length(apobec.high)), rep("low", length(apobec.low)))
  dge.apobec <- DGEList(counts = group.apobec, group = apobec.dgegroup)
  fltrddge.apobec <- dge.apobec[filterByExpr(dge.apobec), ,keep.lib.sizes = FALSE]
  normdge.apobec <- as.data.frame(cpm(calcNormFactors(fltrddge.apobec, method = "TMM")))
  fltrdgene.apobec <- rownames(normdge.apobec)
  apobech.assay <- normdge.apobec[,apobec.high]
  apobecl.assay <- normdge.apobec[,apobec.low]
  adj.p.apobec <- sapply(1:length(fltrdgene.apobec),
                         function(g) wilcox.test(as.numeric(apobech.assay[g,]),
                                                 as.numeric(apobecl.assay[g,]),
                                                 alternative = "two.sided",
                                                 exact = FALSE)$p.value) %>%
    {p.adjust(., method = "BH")}
  deg.apobec <- fltrdgene.apobec[adj.p.apobec < 0.05]
  df.degapobec <- data.frame(highmean = rowMeans(normdge.apobec[deg.apobec, apobec.high]),
                             lowmean = rowMeans(normdge.apobec[deg.apobec, apobec.low]))
  df.degapobec$fold.change = df.degapobec %>% {log2(.$highmean/.$lowmean)}
  return(df.degapobec)
}

DEG.APOBEC <- DEG.APOBEC_load(incl.coldata$paper_APOBEC.induced.mutation.load..PMACD)
str(DEG.APOBEC[abs(DEG.APOBEC$fold.change) > 1,])
