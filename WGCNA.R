## Revisit Chow et al., PNAS 2020 bulk RNA seq data
# CRAN
library(readr)
library(dplyr)
library(WGCNA)
library(spatstat)
library(scico)
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

setwd("~/Documents/BFX_proj/WGCNA_analysis")

### counts data ---
r_count <- read_csv("data/raw/Chow_PNAS_rawcounts.csv")
r_count <- r_count[!duplicated(r_count$gene), ] # remove duplicate genes
c_mtx <- as.matrix(r_count[, -1]) # make numerical matrix
rownames(c_mtx) <- r_count$gene # assign gene name to row names

### metadata ---
meta <- read_csv("data/raw/Chow_PNAS_meta.csv")

### gene synonym reference ---
hs <- org.Hs.eg.db
hs <- AnnotationDbi::select(hs,
                            keys = rownames(c_mtx),
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
ds2_ <- DESeqDataSetFromMatrix(countData = c_mtx, colData = meta, design = ~ 1)
ds2_ <- DESeq(ds2_) # run DESeq2
tnc_mtx <- t(counts(ds2_, normalized = T)) # extract normalized counts, transpose for WGCNA

### WGCNA qc for genes and samples ---
gsg_w <- goodSamplesGenes(tnc_mtx)
if (!gsg_w$allOK){
  if (sum(!gsg_w$goodGenes) > 0) 
    printFlush(paste("Removing genes:", paste(colnames(tnc_mtx)[!gsg_w$goodGenes], collapse = ", ")));
  if (sum(!gsg_w$goodSamples) > 0) 
    printFlush(paste("Removing samples:", paste(rownames(tnc_mtx)[!gsg_w$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  tnc_mtx <- tnc_mtx[gsg_w$goodSamples, gsg_w$goodGenes]
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
tnc_mtx <- tnc_mtx[, order(colVars(as.matrix(tnc_mtx)), decreasing = T)[1:5000]]

### dendrogram of samples ---
#sampleTree <- hclust(dist(tnc_mtx), method = "average")
#par(cex = 0.5)
#plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#dev.off()

### assay possible thresholds ---
sft_w <- pickSoftThreshold(tnc_mtx, powerVector = c(c(1:10), seq(from = 12, to = 20, by = 2)))

### calculate gene similarity ---
adj_w <- adjacency(tnc_mtx, power = sft_w$powerEstimate) # calculate adjacency matrix; how connected each gene is to other genes
tom_w <- TOMsimilarity(adj_w) # calculate topological overlap
geneTree_w <- hclust(as.dist(1 - tom_w), method = "average") # calculate gene tree by dissimilarity

### assign genes to gene module; produces vector of module assignments for each gene ---
mods <- cutreeDynamic(geneTree_w, # clustering
                      distM = 1 - tom_w, # distance matrix
                      deepSplit = 2,
                      pamRespectsDendro = F,
                      minClusterSize = 30) # minimum cluster size is large ~30

### estimate module similarity ---
#MEs_w <- moduleEigengenes(tnc_mtx, colors = labels2colors(mods_w)) # calculate eigengenes (samples) with un-merged ME color assignment
#MEs_ext_w <- MEs_w$eigengenes # extract eigengenes; eigengenes are samples, therefore looking for similarity across all samples
#METree_w <- hclust(as.dist(1-cor(MEs_ext_w)), method = "average") # dissimilarity matrix: col is gene modules, row is samples, construct dendrogram
                    
#plot(METree_w, sub = "", xlab = "", cex = 0.6)
#abline(h = 0.1, col = "red") # dissimilarity threshold e.g. (0.25 is similarity of 0.75)

### merge similar modules ---
mods <- mergeCloseModules(tnc_mtx, labels2colors(mods), cutHeight = 0.1) # merge similar modules

#plotDendroAndColors(geneTree_w, cbind(labels2colors(mods_w), mods_w$colors),
#                   c("Initial modules", "Merged modules"),
#                    dendroLabels = FALSE, hang = 0.03,
#                    cex.colorLabels = 0.6,
#                    addGuide = TRUE, guideHang = 0.05)

### correlate final modules to traits ---
MEs_w <- moduleEigengenes(tnc_mtx, mods$colors)$eigengenes # re-calculate eigengenes with merged ME colors
MEs_w <- orderMEs(MEs_w) # cluster MEs by similarity; put grey (unassigned) at end
moduleTraitCor_w <- cor(MEs_w, meta_w, use = "p") # correlation of sample eigenvalue vs sample meta value; pairwiase complete observations
moduleTraitPvalue_w <- corPvalueStudent(moduleTraitCor_w, nrow(tnc_mtx)) # calculate P value of correlations
moduleTraitFDR_w <- matrix(p.adjust(moduleTraitPvalue_w, method = "BH"),
                           nrow = nrow(moduleTraitPvalue_w),
                           ncol = ncol(moduleTraitPvalue_w)) # adjust by FDR

### clean up ---
mod_trait <- list(cor = moduleTraitCor_w, pval = moduleTraitPvalue_w, fdr = moduleTraitFDR_w)

saveRDS(tnc_mtx, "data/derived/tnc_mtx.rds")
saveRDS(mods, "data/derived/mods.rds")
saveRDS(mod_trait, "data/derived/mod_trait.rds")

rm(list = ls()[grepl("_w$", ls())])

################################################################################
#
# analyze gene modules by GSEA
#
################################################################################

mod_genes <- list()
mod_gsea <- list()

### for loop to pull genes and run three levels of gsea
for(mod_ in unique(mods$colors)){
  genes_ <- colnames(tnc_mtx)[mods$colors == mod_] # get symbols of module genes
  
  entrez_ <- hs$ENTREZID[hs$SYMBOL %in% genes_] # convert symbols to entrezid
  entrez_ <- entrez_[!is.na(entrez_)] # filter missing entrezid
  
  reactome_ <- data.frame(enrichPathway( # run reactome enrichment
    entrez_, hs$ENTREZID,
    organism = "human",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE)); if(nrow(reactome_) > 0){reactome_$source <- "Reactome"}

  kegg_ <- data.frame(enrichKEGG( # run KEGG enrichment
    entrez_, universe = hs$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    use_internal_data = FALSE)); if(nrow(kegg_) > 0){kegg_$source <- "KEGG"}
  
  go_ <- data.frame(enrichGO( # run GO enrichment
    entrez_,
    org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.5,
    pAdjustMethod = "BH",
    minGSSize = 10,
    maxGSSize = 500,
    readable = FALSE,
    pool = FALSE)); if(nrow(go_) > 0){go_$source <- "GO"}
  
  gsea_ <- rbind(reactome_, kegg_, go_)
  
  # add values to list
  mod_genes[[paste0("ME",mod_)]] <- genes_
  mod_gsea[[paste0("ME",mod_)]] <- gsea_
  
  rm(list = ls()[grepl("_$", ls())])
}

saveRDS(mod_genes, "data/derived/mod_genes.rds")
saveRDS(mod_gsea, "data/derived/mod_gsea.rds")

################################################################################
#
# visualize gene module correlation to traits
#
################################################################################

### plot module eigenvalue vs trait ---
# heatmap fill scale
col_fun <- circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                                c(scico(5, palette = "vanimo")[1], scico(5, palette = "vanimo")[2],
                                  scico(5, palette = "vanimo")[3],
                                  scico(5, palette = "vanimo")[4], scico(5, palette = "vanimo")[5]))

# heatmap module name to GSEA
mod_gsea_short <- data.frame(row.names = names(mod_gsea), mod = names(mod_gsea))
mod_gsea_short$paths <- NULL

for(mod_ in names(mod_gsea)){
  if(nrow(mod_gsea[[mod_]]) > 0){
    gsea_ <- mod_gsea[[mod_]]
    gsea_$s_d <- paste0(gsea_$source, ": ", substr(gsea_$Description, 1, 40), "...")
    gsea_ <- gsea_[order(gsea_$pvalue), ]
    
    paths_ <- paste0(gsea_$s_d[1], "\n", gsea_$s_d[2])
    mod_gsea_short[mod_, "paths"] <- paths_
  }
}

mod_gsea_short <- mod_gsea_short[rownames(mod_trait$cor), ]

# heatmap row annotation
anno_row <- HeatmapAnnotation(GSEA = anno_text(mod_gsea_short$paths, gp = gpar(fontsize = 10)),
                              which = "row")

# plot heatmap
png("figs/heat1.png", height = 600, width = 700, res = 100)
Heatmap(mod_trait$cor,
        name = "Correlation",
        col = col_fun,
        right_annotation = anno_row,
        show_row_names = F,
        width = unit(ncol(mod_trait$cor) * 8, "mm"),
        height = unit(nrow(mod_trait$cor) * 8, "mm"),
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y, r = abs(-log10(mod_trait$fdr)[i, j]/1.5) * min(unit.c(width, height)), # size is FDR
                      gp = gpar(fill = col_fun(mod_trait$cor[i, j]), col = NA)) # color is correlation
          },
        cluster_columns = F)
dev.off()

png("figs/heat2.png", height = 600, width = 700, res = 100)
Heatmap(mod_trait$cor,
        name = "Correlation",
        col = col_fun,
        right_annotation = anno_row,
        show_row_names = F,
        width = unit(ncol(mod_trait$cor) * 8, "mm"),
        height = unit(nrow(mod_trait$cor) * 8, "mm"),
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.circle(x = x, y = y, r = abs(-log10(mod_trait$fdr)[i, j]/1.5) * min(unit.c(width, height)), # size is FDR
                      gp = gpar(fill = col_fun(mod_trait$cor[i, j]), col = NA)) # color is correlation
          if(mod_trait$fdr[i, j] > 0.5){ # greyed out does not pass FDR threshold
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = scico(3, palette = "vanimo")[2]))
          } else {
            grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          }
        },
        cluster_columns = F)
dev.off()






































