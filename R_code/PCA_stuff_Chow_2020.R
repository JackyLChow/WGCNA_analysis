# Playing around w/ PCA stuff

# CRAN
library(readr)
library(dplyr)
library(WGCNA)
library(scico)
library(circlize)
library(matrixStats)
library(stringr)
library(ggplot2)

# Bioconductor
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ReactomePA)
library(clusterProfiler)

################################################################################
#
# directories
#
################################################################################

setwd("~/Documents/BFX_proj/WGCNA_analysis")

################################################################################
#
# load data
#
################################################################################

### c_mtx: normalized counts data---
count <- read_csv("_input/Chow_PNAS_2020/Chow_PNAS_normcounts.csv")
c_mtx <- as.matrix(count[, 2:ncol(count)])
rownames(c_mtx) <- count$gene

### meta: metadata, ---
meta <- data.frame(read_csv("_input/Chow_PNAS_2020/Chow_PNAS_metashort.csv"))
rownames(meta) <- str_replace(meta$sample, "-", ".")
meta <- meta[, 2:ncol(meta)]

### ensure metadata match ---
meta <- meta[colnames(c_mtx), ]

### gene synonym reference ---
hs <- org.Hs.eg.db
hs <- AnnotationDbi::select(hs,
                            keys = rownames(c_mtx),
                            columns = c("ENTREZID"),
                            keytype = "SYMBOL")
hs <- hs[!duplicated(hs$SYMBOL), ]

# clean up
rm(count)

################################################################################
#
# PCA stuff
#
################################################################################

pca <- prcomp(t(c_mtx)) # run PCA
str(pca) # look at structure

pt_df <- data.frame(meta,
                    PC1 = pca$x[, 1],
                    PC2 = pca$x[, 2])

ggplot(pt_df, aes(PC1, PC2, color = factor(SBRT))) +
  geom_point() +
  geom_hline(yintercept = median(pt_df[pt_df$SBRT == 0, "PC2"]), color = "red") +
  geom_hline(yintercept = median(pt_df[pt_df$SBRT == 1, "PC2"]), color = "blue")

eig <- pca$rotation
eig_PC2 <- eig[, "PC2"]
eig_PC2 <- eig_PC2[order(eig_PC2, decreasing = T)]

head(eig_PC2)
tail(eig_PC2)

pt_df <- cbind(pt_df,
               t(c_mtx[c("LYZ", "PTPRC", "HSPA1A"), ]))

ggplot(pt_df, aes(PC1, PC2, color = factor(SBRT), size = HSPA1A)) + geom_point()











