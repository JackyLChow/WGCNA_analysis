## Revisit Chow et al., PNAS 2020 bulk RNA seq data
# CRAN
library(readr)
library(dplyr)
library(WGCNA)

# Bioconductor
library(DESeq2)

################################################################################
#
# load data
#
################################################################################

r_count <- read_csv("Documents/GitHub/Revisit_RNAseq_2020_Chow_PNAS/Data/Chow_PNAS_rawcounts.csv")
r_count <- r_count[!duplicated(r_count$gene), ]
cnt_mtx <- as.matrix(r_count[, -1])
rownames(cnt_mtx) <- r_count$gene

sp_meta <- read_csv("Documents/GitHub/Revisit_RNAseq_2020_Chow_PNAS/Data/Chow_PNAS_meta.csv")

