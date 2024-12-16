#! /usr/bin/Rscript

library(TSCAN)

# ---------- parameters ----------
args <- commandArgs(trailingOnly=TRUE)

counts_file <- args[1]
output_dir <- args[2]

count_matrix <- read.csv(counts_file, row.names=1)
exp <- as.matrix(count_matrix)
procdata <- preprocess(exp, minexpr_value = 0, minexpr_percent = 0.001, cvcutoff = 0.01)
lpsmclust <- exprmclust(procdata)
tscan_order = TSCANorder(lpsmclust)
tscan_rank = rank(tscan_order[colnames(count_matrix),"Pseudotime"], ties.method = 'random', na.last = 'keep')
tscan_norm = scales::rescale(tscan_rank,to = c(0,1))
output <- data.frame(tscan_norm)
rownames(output) <- colnames(count_matrix)

write.csv(output, file.path(output_dir, "tscan.csv"))