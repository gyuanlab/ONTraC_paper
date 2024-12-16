import sys
import numpy as np
import pandas as pd

import anndata as adata
import scanpy as sc
import squidpy as sq
from SpaceFlow import SpaceFlow

exp_file, meta_file, output_dir = sys.argv[1:]

exp_df = pd.read_csv(exp_file, index_col=0)
meta_info_df = pd.read_csv(meta_file, index_col=0)

ad = adata.AnnData(exp_df.T)
ad.obs = ad.obs.join(meta_info_df)
ad.obsm['spatial'] = ad.obs[['x', 'y']].values

sf = SpaceFlow.SpaceFlow(adata=ad)
sf.preprocessing_data(n_top_genes=3000)
sf.train(spatial_regularization_strength=0.1, z_dim=50, lr=1e-3, epochs=1000, max_patience=50, min_stop=100, random_seed=42, gpu=0, regularization_acceleration=True, edge_subset_sz=1000000)
sf.segmentation(domain_label_save_filepath=f"{output_dir}/domains.tsv", n_neighbors=50, resolution=1.0)
sf.pseudo_Spatiotemporal_Map(pSM_values_save_filepath=f"{output_dir}/pSM_values.tsv", n_neighbors=20, resolution=1.0)