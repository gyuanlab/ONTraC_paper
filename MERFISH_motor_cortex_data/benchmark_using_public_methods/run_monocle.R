#! /usr/bin/Rscript

library(monocle3)

# ---------- parameters ----------
args <- commandArgs(trailingOnly=TRUE)

cds_file <- args[1]
output_dir <- args[2]

cds <- readRDS(cds_file)
cds <- preprocess_cds(cds, method = "PCA", num_dim = 20)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "PCA")
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition=TRUE)
root.use <- rownames(cds@colData[cds@colData$sim_time==min(cds@colData$sim_time),])
cds <- order_cells(cds, root_cells=root.use)
pseudotime <- pseudotime(cds) 
monocle3_rank <- rank(pseudotime,ties.method = 'random', na.last = 'keep')
monocle3_norm <- scales::rescale(monocle3_rank,to = c(0,1))
cds@colData$monocle3_norm <- monocle3_norm

write.csv(cds@colData, file.path(output_dir, "monocle.csv"))