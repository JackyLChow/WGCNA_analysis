## Revisit Chow et al., PNAS 2020 bulk RNA seq data
# CRAN
library(readr)
library(dplyr)
library(WGCNA)

# Bioconductor
library(DESeq2)
library(org.Hs.eg.db)
library(ReactomePA)
library(clusterProfiler)

################################################################################
#
# load data
#
################################################################################

setwd("~/Documents/bfx_proj/Revisit_RNAseq_2020_Chow_PNAS/")

### counts data ---
r_count <- read_csv("Data/Chow_PNAS_rawcounts.csv")
r_count <- r_count[!duplicated(r_count$gene), ] # remove dupilicate genes
cnt_mtx <- as.matrix(r_count[, -1]) # make numerical matrix
rownames(cnt_mtx) <- r_count$gene

### metadata ---
sp_meta <- read_csv("Data/Chow_PNAS_meta.csv")

### gene synonym reference ---
hs <- org.Hs.eg.db
hs <- AnnotationDbi::select(hs,
                            keys = rownames(cnt_mtx),
                            columns = c("ENTREZID"),
                            keytype = "SYMBOL")
hs <- hs[!duplicated(hs$SYMBOL), ]























