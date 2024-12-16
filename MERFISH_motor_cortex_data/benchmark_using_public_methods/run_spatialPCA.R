#! /usr/bin/Rscript

library(SpatialPCA)
library(slingshot)

# ---------- parameters ----------
args <- commandArgs(trailingOnly=TRUE)

exp_file <- args[1]
meta_file <- args[2]
output_dir <- args[3]

exp_df <- read.csv(exp_file, row.names = 1)
rownames(exp_df) <- gsub("_", "-", rownames(exp_df))
meta_info_df <- read.csv(meta_file, row.names = 1)

xy_coods <- as.matrix(meta_info_df[colnames(exp_df),c('x', 'y')])
spatialPCA_obj <- CreateSpatialPCAObject(counts=as.matrix(exp_df), location=xy_coods, project = "SpatialPCA",gene.type="spatial",sparkversion="spark",numCores_spark=5,gene.number=3000, customGenelist=NULL,min.loctions = 2, min.features=2)
spatialPCA_obj <- SpatialPCA_buildKernel(spatialPCA_obj, kerneltype="gaussian", bandwidthtype="SJ",bandwidth.set.by.user=NULL)
spatialPCA_obj <- SpatialPCA_EstimateLoading(spatialPCA_obj,fast=FALSE,SpatialPCnum=10) 
spatialPCA_obj <- SpatialPCA_SpatialPCs(spatialPCA_obj, fast=FALSE)
clusterlabel= walktrap_clustering(clusternum=6,
                                  latent_dat=spatialPCA_obj@SpatialPCs,
                                  knearest=round(sqrt(dim(spatialPCA_obj@SpatialPCs)[2]))
                                 ) 
sce <- SingleCellExperiment(assays = as.matrix(exp_df))
reducedDims(sce) <- SimpleList(DRM = t(spatialPCA_obj@SpatialPCs))
colData(sce)$clusterlabel <- factor(clusterlabel)    
sce  <-slingshot(sce, clusterLabels = 'clusterlabel', reducedDim = 'DRM',start.clus="3" ) 
# in this data we set white matter region as start cluster, one can change to their preferred start region 

# summary(sce@colData@listData)
pseudotime_traj1 <- sce@colData@listData$slingPseudotime_1 # in this data only one trajectory was inferred

write.csv(pseudotime_traj1, file.path(output_dir, "spatial_pca.csv"))