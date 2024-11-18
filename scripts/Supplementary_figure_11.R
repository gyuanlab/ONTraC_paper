# Using https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast
# Published at https://www.nature.com/articles/s41467-023-43458-x

cell_types_table <- readxl::read_xlsx("data/xenium/Cell_Barcode_Type_Matrices.xlsx",
                                      sheet = 4)

download.file(url = "https://cf.10xgenomics.com/samples/xenium/1.0.1/Xenium_FFPE_Human_Breast_Cancer_Rep1/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip",
              "data/xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip")

untar("data_xenium/Xenium_FFPE_Human_Breast_Cancer_Rep1_outs.zip",
      exdir = "data/xenium")

##%######################################################%##
#                                                          #
####                   Create object                    ####
#                                                          #
##%######################################################%##

library(Giotto)

expression_matrix <- get10Xmatrix("data/xenium/outs/cell_feature_matrix/",
                       gene_column_index = 2)

expression_matrix <- expression_matrix$`Gene Expression`

spatial_locs <- data.table::fread("data/xenium/outs/cells.csv.gz")
spatial_locs <- spatial_locs[, c("x_centroid", "y_centroid", "cell_id")]
colnames(spatial_locs) <- c("sdimx", "sdimy", "cell_ID")

instructions <- createGiottoInstructions(python_path = NULL,
                                         save_plot = TRUE,
                                         show_plot = FALSE,
                                         return_plot = FALSE,
                                         save_dir = "results_xenium")

xenium <- createGiottoObject(expression = expression_matrix,
                             spatial_locs = spatial_locs,
                             instructions = instructions)

cell_types <- cell_types_table$Cluster

xenium <- addCellMetadata(xenium,
                          new_metadata = cell_types)

##%######################################################%##
#                                                          #
####                         QC                         ####
#                                                          #
##%######################################################%##

filterDistributions(xenium,
                    detection = "feats")

filterDistributions(xenium,
                    detection = "cells")

##%######################################################%##
#                                                          #
####                      Filter                        ####
#                                                          #
##%######################################################%##

xenium <- filterGiotto(xenium,
                       feat_det_in_min_cells = 100,
                       min_det_feats_per_cell = 10)

# Number of cells removed:  4215  out of  167780 
# Number of feats removed:  0  out of  313 


##%######################################################%##
#                                                          #
####                    Save object                     ####
#                                                          #
##%######################################################%##

saveGiotto(xenium, "xenium_object")

##%######################################################%##
#                                                          #
####                Create ONTraC table                 ####
#                                                          #
##%######################################################%##

metadata <- pDataDT(xenium)

coordinates <- getSpatialLocations(xenium,
                                   output = "data.table")

# create data table
dataset_table <- data.frame(Cell_ID = metadata$cell_ID,
                            Sample = rep("sample1_rep1", 163565),
                            Cell_Type = metadata$cell_types,
                            x = coordinates$sdimx,
                            y = coordinates$sdimy)


# subset region
write.table(dataset_table[dataset_table$x > 3000 & dataset_table$x < 4500  & dataset_table$y > 1000 & dataset_table$y < 2500, ], 
            "dataset.csv", sep = ",", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

##%######################################################%##
#                                                          #
####                     Run ONTraC                     ####
#                                                          #
##%######################################################%##

conda activate ONTraC_v1

## bash code
ONTraC -d dataset.csv \
--preprocessing-dir selection/data_xenium/preprocessing_dir \
--GNN-dir selection/output_xenium/GNN \
--NTScore-dir selection/output_xenium/NTScore \
--device cpu --epochs 300 -s 42 --patience 100 \
--min-delta 0.001 --min-epochs 50 --lr 0.03 --hidden-feats 4 \
-k 4 --n-neighbors 200 --beta 0.3 | tee xenium_final.log

ONTraC_analysis -d dataset.csv \
--preprocessing-dir selection/data_xenium/preprocessing_dir \
--GNN-dir selection/output_xenium/GNN \
--NTScore-dir selection/output_xenium/NTScore \
-o selection/analysis_xenium -r True -l xenium_final.log


##%######################################################%##
#                                                          #
####                  Spatial Plot                      ####
#                                                          #
##%######################################################%##

dataset <- read.csv("dataset.csv")

xenium <- subsetGiotto(xenium,
                       cell_ids = dataset3$Cell_ID)

##%######################################################%##
#                                                          #
####                       Plots                        ####
#                                                          #
##%######################################################%##

spatPlot2D(xenium,
           cell_color = "cell_types",
           point_size = 2)

