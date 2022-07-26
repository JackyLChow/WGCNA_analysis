library(spatstat)
library(readr)

### normalized counts data ---
count <- read_csv("~/Documents/BFX_proj/_x_x_data/immotion/raw/counts.csv")
c_mtx <- as.matrix(count[, 2:ncol(count)])
rownames(c_mtx) <- count$gene

### load metadata, parent is EGAD00001006617 ---
meta <- data.frame(read_csv("~/Documents/BFX_proj/_x_x_data/immotion/raw/anno.csv"))
rownames(meta) <- meta$RNASEQ_SAMPLE_ID
meta <- meta[colnames(c_mtx), ]

### dummify metadata for WGCNA ---
meta <- meta[, c("RNASEQ_SAMPLE_ID", "ARM", "OBJECTIVE_RESPONSE")]
meta$OBJECTIVE_RESPONSE <- factor(meta$OBJECTIVE_RESPONSE, levels = c("CR", "PR", "SD", "PD", "NE"))
meta$ARM_OR <- factor(paste0(meta$ARM, "_", meta$OBJECTIVE_RESPONSE),
                              levels = c(paste0("sunitinib", "_", c("CR", "PR", "SD", "PD", "NE")),
                                         paste0("atezo_bev", "_", c("CR", "PR", "SD", "PD", "NE"))))
meta <- data.frame(RNASEQ_SAMPLE_ID = meta$RNASEQ_SAMPLE_ID,
                   dummify(meta$OBJECTIVE_RESPONSE),
                   dummify(meta$ARM_OR))

write.csv(meta, "~/Documents/BFX_proj/_x_x_data/immotion/derived/immotion_metashort.csv", row.names = F)
