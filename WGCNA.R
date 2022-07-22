# Generalized WGCNA workflow
## input: normalized counts data, curated metadata for WGCNA
## output: image of WGCNA heatmap, table of WGCNA gene modules, table of WGCNA GSEA results

# CRAN
library(readr)
library(dplyr)
library(WGCNA)
library(scico)
library(circlize)
library(matrixStats)

# Bioconductor
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ReactomePA)
library(clusterProfiler)

################################################################################
#
# data input
#
################################################################################

setwd("~/Documents/BFX_proj/WGCNA_analysis")

out_tab_dir <- "_output/Chowtput_PNAS_2020/tab/"
out_rds_dir <- "_output/Chowtput_PNAS_2020/rds/"
out_fig_dir <- "_output/Chowtput_PNAS_2020/fig/"

### normalized counts data ---
count <- read_csv("_input/Chow_PNAS_2020/Chow_PNAS_normcounts.csv")
c_mtx <- as.matrix(count[, 2:ncol(count)])
rownames(c_mtx) <- count$gene

### trim count matrix for speed, if NULL then no trim ---
trim_c_mtx <- 5000

### metadata ---
meta <- read_csv("_input/Chow_PNAS_2020/Chow_PNAS_metashort.csv")

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
# normalize and prepare for WGCNA
#
################################################################################

tnc_mtx <- t(c_mtx) # transpose counts for WGCNA

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

# clean up
rm(gsg_w)

################################################################################
#
# run gently modified WGCNA workflow
#
################################################################################

### filter top variant genes for speed ---
if(is.numeric(trim_c_mtx)){
  tnc_mtx <- tnc_mtx[, order(colVars(as.matrix(tnc_mtx)), decreasing = T)[1:trim_c_mtx]]
}

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
moduleTraitCor_w <- cor(MEs_w, meta, use = "p") # correlation of sample eigenvalue vs sample meta value; pairwiase complete observations
moduleTraitPvalue_w <- corPvalueStudent(moduleTraitCor_w, nrow(tnc_mtx)) # calculate P value of correlations
moduleTraitFDR_w <- matrix(p.adjust(moduleTraitPvalue_w, method = "BH"),
                           nrow = nrow(moduleTraitPvalue_w),
                           ncol = ncol(moduleTraitPvalue_w)) # adjust by FDR

### correlate genes to final modules ---
geneModuleMembership_w <- cor(tnc_mtx, MEs_w, use = "p")
MMPvalue_w <- corPvalueStudent(as.matrix(geneModuleMembership_w), nrow(tnc_mtx))
MMPFDR_w <- matrix(p.adjust(MMPvalue_w, method = "BH"),
                   nrow = nrow(geneModuleMembership_w),
                   ncol = ncol(geneModuleMembership_w)) # adjust by FDR

#geneTraitSignificance_w <- as.data.frame(cor(tnc_mtx, weight, use = "p"))
#GSPvalue_w <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_w), nrow(tnc_mtx)))

### clean up ---
mod_trait <- list(cor = moduleTraitCor_w, pval = moduleTraitPvalue_w, fdr = moduleTraitFDR_w)
gene_mod <- list(gMM = geneModuleMembership_w, pval = MMPvalue_w, fdr = MMPFDR_w)

#saveRDS(tnc_mtx, paste0(out_rds_dir, "tnc_mtx.rds"))
#saveRDS(mods, paste0(out_rds_dir, "mods.rds"))
#saveRDS(mod_trait, paste0(out_rds_dir, "mod_trait.rds"))
#saveRDS(gene_mod, paste0(out_rds_dir, "gene_mod.rds"))

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
#saveRDS(mod_genes, paste0(out_rds_dir, "mod_genes.rds"))
#saveRDS(mod_gsea, paste0(out_rds_dir, "mod_gsea.rds"))

### export genes and gsea as csv ---
mod_genes_ <- data.frame(row.names = names(mod_genes),
                         mod_name = names(mod_genes))
mod_genes_$mod_genes <- NULL
for(mod_ in mod_genes_$mod_name){
  mod_genes_[mod_, "mod_genes"] <- paste(mod_genes[[mod_]], collapse = ", ")
}
#write.csv(mod_genes_, paste0(out_tab_dir, "mod_genes_long.csv"))
#write.csv(bind_rows(mod_gsea, .id = "column_label"), paste0(out_tab_dir, "mod_gsea_long.csv"))

rm(mod_genes_) # clean up

################################################################################
#
# visualize gene module correlation to traits
#
################################################################################

### load necessary rds files ---
mod_genes <- readRDS(paste0(out_rds_dir, "mod_genes.rds"))
mod_gsea <- readRDS(paste0(out_rds_dir, "mod_gsea.rds"))
mod_trait <- readRDS(paste0(out_rds_dir, "mod_trait.rds"))
gene_mod <- readRDS(paste0(out_rds_dir, "gene_mod.rds"))

### heatmap parameters ---
# heatmap fill scale
col_fun <- circlize::colorRamp2(c(-1, -0.5, 0, 0.5, 1),
                                c(scico(5, palette = "vanimo")[1], scico(5, palette = "vanimo")[2],
                                  scico(5, palette = "vanimo")[3],
                                  scico(5, palette = "vanimo")[4], scico(5, palette = "vanimo")[5]))

# heatmap module annotation with shortened GSEA and gene lists
mod_gene_dat_short <- data.frame(row.names = names(mod_gsea), mod = names(mod_gsea))
mod_gene_dat_short$paths <- NULL
mod_gene_dat_short$genes <- NULL

for(mod_ in names(mod_gsea)){
  if(nrow(mod_gsea[[mod_]]) > 0){
    gsea_ <- mod_gsea[[mod_]]
    gsea_$s_d <- paste0(gsea_$source, ": ", substr(gsea_$Description, 1, 40), "...")
    gsea_ <- gsea_[order(gsea_$pvalue), ]
    
    paths_ <- paste0(gsea_$s_d[1], "\n", gsea_$s_d[2])
    mod_gene_dat_short[mod_, "paths"] <- paths_
  }
  
  gene_ <- gene_mod$gMM[mod_genes[[mod_]], mod_] # trim module membership table to genes assigned to module
  gene_ <- gene_[order(gene_)] # order by membership
  gene_ <- paste0("pos: ", paste(names(head(gene_, 5)), collapse = ", "), "\n",
                  "neg: ", paste(names(tail(gene_, 5)), collapse = ", "))
  mod_gene_dat_short[mod_, "genes"] <- gene_
}

mod_gene_dat_short <- mod_gene_dat_short[rownames(mod_trait$cor), ] # rearrange to match correlation matrix

### plot module eigenvalue vs trait ---
# plot heatmap1; correlation only
anno_row <- HeatmapAnnotation(ME = rownames(mod_trait$cor),
                              col = list(ME = setNames(str_remove(rownames(mod_trait$cor), "ME"), rownames(mod_trait$cor))),
                              which = "row",
                              show_legend = F)

png(paste0(out_fig_dir, "heat1.png"), height = nrow(mod_trait$cor) * 37.5, width = ncol(mod_trait$cor) * 37.5 + 300, res = 100)
set.seed(415); Heatmap(mod_trait$cor,
                       name = "Correlation",
                       col = col_fun,
                       right_annotation = anno_row,
                       show_row_names = T,
                       row_names_side = "right",
                       width = unit(ncol(mod_trait$cor) * 8, "mm"),
                       height = unit(nrow(mod_trait$cor) * 8, "mm"),
                       cluster_columns = F)
dev.off()

# plot heatmap2; correlation, significance, top 2 GSEA
anno_row <- HeatmapAnnotation(GSEA = anno_text(mod_gene_dat_short$paths, gp = gpar(fontsize = 10)),
                              which = "row")

png(paste0(out_fig_dir, "heat2.png"), height = nrow(mod_trait$cor) * 37.5, width = ncol(mod_trait$cor) * 37.5 + 500, res = 100)
set.seed(415); Heatmap(mod_trait$cor,
        name = "Correlation",
        col = col_fun,
        right_annotation = anno_row,
        show_row_names = F,
        show_row_dend = F,
        row_names_side = "right",
        width = unit(ncol(mod_trait$cor) * 8, "mm"),
        height = unit(nrow(mod_trait$cor) * 8, "mm"),
        rect_gp = gpar(type = "none"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
          grid.circle(x = x, y = y, r = abs(-log10(mod_trait$fdr)[i, j]/max(-log10(mod_trait$fdr))) * min(unit.c(width, height) * 1.5), # size is FDR
                      gp = gpar(fill = col_fun(mod_trait$cor[i, j]), col = NA)) # color is correlation
          },
        cluster_columns = F)
dev.off()

# plot heatmap3; correlation, significance, top 5 positive and negative correlated genes for each module
anno_row <- HeatmapAnnotation(GSEA = anno_text(mod_gene_dat_short$genes, gp = gpar(fontsize = 10)),
                              which = "row")

png(paste0(out_fig_dir, "heat3.png"), height = nrow(mod_trait$cor) * 37.5, width = ncol(mod_trait$cor) * 37.5 + 500, res = 100)
Heatmap(mod_trait$cor,
        name = "Correlation",
        col = col_fun,
        right_annotation = anno_row,
        show_row_names = F,
        show_row_dend = F,
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






































