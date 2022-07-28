library(DESeq2)
library(spatstat)
library(readr)

### raw counts ---
r_count <- read_csv("~/Documents/BFX_proj/_data/Chow_PNAS_2020/Chow_PNAS_rawcounts.csv")
r_count <- r_count[!duplicated(r_count$gene), ] # remove duplicate genes
c_mtx <- as.matrix(r_count[, -1]) # make numerical matrix
rownames(c_mtx) <- r_count$gene # assign gene name to row names

### metadata ---
meta <- read_csv("~/Documents/BFX_proj/_data/Chow_PNAS_2020/Chow_PNAS_meta.csv")

### normalize counts in DESeq2 ---
ds2_ <- DESeqDataSetFromMatrix(countData = c_mtx, colData = meta, design = ~ 1) # make DESeq2 object
ds2_ <- DESeq(ds2_) # run DESeq2

### extract normalized counts ---
n_count <- counts(ds2_, normalized = T)
n_count <- data.frame(gene = rownames(n_count), n_count)

write.csv(n_count, "~/Documents/BFX_proj/WGCNA_analysis/_input/Chow_PNAS_2020/Chow_PNAS_normcounts.csv", row.names = F)

### dummify metadata for WGCNA ---
meta <- meta[, c("sample", "age", "sex", "path_T_stage", "treatment")] # traits of interest
meta <- data.frame(meta[, c("sample", "age")],
                   dummify(meta$sex),
                   dummify(meta$treatment))
                     

write.csv(meta, "~/Documents/BFX_proj/WGCNA_analysis/_input/Chow_PNAS_2020/Chow_PNAS_metashort.csv", row.names = F)
