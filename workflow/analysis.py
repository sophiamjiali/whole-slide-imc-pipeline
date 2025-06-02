"""
Script:          analysis.py
Purpose:         Downstream analysis of WS-IMC data
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            05-26-2025
"""

import scanpy as sc
import numpy as np
import pandas as pd

def compute_marker_correlation(adata, n_components):
    # Computes correlations between markers across the whole-slide

    # Perform PCA with the specified number of components
    sc.tl.pca(adata, n_components)
    components = adata.obsm['X_pca'][:, :n_components]
    
    # Calculate correlations between markers and each component
    correlations = pd.Dataframe(index = adata.var_names,
                                columns = [f'PC{i+1}' for i in range(n_components)])
    for i, marker in enumerate(adata.var_names):
        for j in range(n_components):
            correlations.iloc[i, j] = np.corrcoef(adata.X[:, i], components[:, j])[0, 1]

    return correlations


def compute_window_centroid(adata):
    # Computes the centroid of each window in the whole-slide image

    return

