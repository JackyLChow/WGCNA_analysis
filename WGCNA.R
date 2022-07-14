## Revisit Chow et al., PNAS 2020 bulk RNA seq data
# CRAN
library(readr)
library(dplyr)
library(WGCNA)
library(spatstat)
library(circlize)

# Bioconductor
library(DESeq2)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ReactomePA)
library(clusterProfiler)

################################################################################
#
# load parent data
#
################################################################################

setwd("~/Documents/BFX_proj/Revisit_RNAseq_2020_Chow_PNAS/")

### counts data ---
r_count <- read_csv("Data/Chow_PNAS_rawcounts.csv")
r_count <- r_count[!duplicated(r_count$gene), ] # remove duplicate genes
cmtx <- as.matrix(r_count[, -1]) # make numerical matrix
rownames(cmtx) <- r_count$gene

### metadata ---
meta <- read_csv("Data/Chow_PNAS_meta.csv")

### gene synonym reference ---
hs <- org.Hs.eg.db
hs <- AnnotationDbi::select(hs,
                            keys = rownames(cmtx),
                            columns = c("ENTREZID"),
                            keytype = "SYMBOL")
hs <- hs[!duplicated(hs$SYMBOL), ]

# clean up
rm(r_count)

################################################################################
#
# normalize and prepare for WGCNA
#
################################################################################

### normalize counts in DESeq2 ---
ds2_ <- DESeqDataSetFromMatrix(countData = cmtx, colData = meta, design = ~ 1)
ds2_ <- DESeq(ds2_) # run DESeq2
cmtx_w <- t(counts(ds2_, normalized = T)) # extract normalized counts, transpose for WGCNA

### WGCNA qc for genes and samples ---
gsg_w <- goodSamplesGenes(cmtx_w)
if (!gsg_w$allOK){
  if (sum(!gsg_w$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(colnames(cmtx_w)[!gsg_w$goodGenes], collapse = ", ")));
  if (sum(!gsg_w$goodSamples) > 0) 
    printFlush(paste("Removing samples:", paste(rownames(cmtx_w)[!gsg_w$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  cmtx_w <- cmtx_w[gsg_w$goodSamples, gsg_w$goodGenes]
}

### shorten and dummify metadata
meta_w <- meta[, c("age", "sex", "path_T_stage", "treatment")]
meta_w <- data.frame(meta_w[, "age"],
                     dummify(meta_w$sex),
                     dummify(meta_w$treatment))

# clean up
rm(ds2_, gsg_w)

################################################################################
#
# run modified WGCNA
#
################################################################################

### filter top variant genes for speed ---
cmtx_w <- cmtx_w[, order(colVars(as.matrix(cmtx_w)), decreasing = T)[1:5000]]

### dendrogram of samples ---
#sampleTree <- hclust(dist(cmtx_w), method = "average")
#par(cex = 0.5)
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#dev.off()

### assay possible thresholds ---
sft_w <- pickSoftThreshold(cmtx_w, powerVector = c(c(1:10), seq(from = 12, to = 20, by = 2)))

### calculate gene similarity ---
adj_w <- adjacency(cmtx_w, power = sft_w$powerEstimate) # calculate adjacency matrix; how connected each gene is to other genes
tom_w <- TOMsimilarity(adj_w) # calculate topological overlap
geneTree_w <- hclust(as.dist(1 - tom_w), method = "average") # calculate gene tree by dissimilarity

### assign genes to gene module; produces vector of module assignments for each gene ---
mods_w <- cutreeDynamic(geneTree_w, # clustering
                        distM = 1 - tom_w, # distance matrix
                        deepSplit = 2,
                        pamRespectsDendro = F,
                        minClusterSize = 30) # minimum cluster size is large ~30

### estimate module similarity ---
MEs_w <- moduleEigengenes(cmtx_w, colors = labels2colors(mods_w)) # calculate eigengenes (samples) with un-merged ME color assignment
MEs_ext_w <- MEs_w$eigengenes # extract eigengenes; eigengenes are samples, therefore looking for similarity across all samples
METree_w <- hclust(as.dist(1-cor(MEs_ext_w)), # dissimilarity matrix: col is gene modules, row is samples
                   method = "average") # construct dendrogram
plot(METree_w, sub = "", xlab = "", cex = 0.6)
abline(h = 0.1, col = "red") # dissimilarity threshold e.g. (0.25 is similarity of 0.75)

### merge similar modules ---
mods_merg_w <- mergeCloseModules(cmtx_w, labels2colors(mods_w), cutHeight = 0.1) # merge similar modules

plotDendroAndColors(geneTree_w, cbind(labels2colors(mods_w), mods_merg_w$colors),
                    c("Initial modules", "Merged modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    cex.colorLabels = 0.6,
                    addGuide = TRUE, guideHang = 0.05)

### correlate final modules to traits ---
MEs_w <- moduleEigengenes(cmtx_w, mods_merg_w$colors)$eigengenes # re-calculate eigengenes with merged ME colors
MEs_w <- orderMEs(MEs_w) # cluster MEs by similarity; put grey (unassigned) at end
moduleTraitCor_w <- cor(MEs_w, meta_w, use = "p") # correlation of sample eigenvalue vs sample meta value; pairwiase complete observations
moduleTraitPvalue_w <- corPvalueStudent(moduleTraitCor_w, nrow(cmtx_w)) # calculate P value of correlations
moduleTraitFDR_w <- matrix(p.adjust(moduleTraitPvalue_w, method = "BH"),
                           nrow = nrow(moduleTraitPvalue_w),
                           ncol = ncol(moduleTraitPvalue_w)) # adjust by FDR

### plot module eigenvalue vs trait ---
col_fun <- circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                                c(scico(5, palette = "vanimo")[1], scico(5, palette = "vanimo")[2],
                                  scico(5, palette = "vanimo")[3],
                                  scico(5, palette = "vanimo")[4], scico(5, palette = "vanimo")[5]))

Heatmap(moduleTraitCor_w,
        name = "Correlation",
        col = col_fun,
        width = unit(ncol(moduleTraitCor_w) * 5, "mm"),
        height = unit(nrow(moduleTraitCor_w) * 5, "mm"),
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y, r = abs(-log10(moduleTraitFDR_w)[i, j]/2) * min(unit.c(width, height)), # size is FDR
                      gp = gpar(fill = col_fun(moduleTraitCor_w[i, j]), col = NA)) # color is correlation
          },
        cluster_columns = F)

Heatmap(moduleTraitCor_w,
        name = "Correlation",
        col = col_fun,
        width = unit(ncol(moduleTraitCor_w) * 5, "mm"),
        height = unit(nrow(moduleTraitCor_w) * 5, "mm"),
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(-log10(moduleTraitFDR_w)[i, j]/2) * min(unit.c(width, height)), # size is FDR
                      gp = gpar(fill = col_fun(moduleTraitCor_w[i, j]), col = NA)) # color is correlation
          if(moduleTraitFDR_w[i, j] > 0.5){ # greyed out does not pass FDR threshold
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = "grey30"))
          } else {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          }
        },
        cluster_columns = F)






































