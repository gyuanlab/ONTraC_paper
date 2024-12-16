#! /usr/bin/Rscript

library(scater)
library(destiny)

# ---------- parameters ----------
args <- commandArgs(trailingOnly=TRUE)

sce_file <- args[1]
output_dir <- args[2]

sce <- readRDS(sce_file)
dm.use <- DiffusionMap(data = sce)
root.use <- rownames(sce@colData)[which(sce@colData$sim_time == 0)[1]]
dpt.use <- DPT(dm.use,tips = match(root.use,rownames(dm.use@eigenvectors)))
dpt_rank <- rank(dpt.use$dpt,ties.method = 'random', na.last = 'keep')
sce@colData$dpt_norm <- scales::rescale(dpt_rank,to = c(0,1))

write.csv(sce@colData, file.path(output_dir, "DPT.csv"))