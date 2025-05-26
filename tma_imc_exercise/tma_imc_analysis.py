import tifffile as tf
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import os

def load_tiff_image(file_path):
    image = tf.imread(file_path)
    return image


def get_cell_expression(image, mask):
    # return a dataframe (cells x proteins) of the protein counts per cell

    # Loop through each cell label and calculate the expression
    cell_labels = np.unique(mask)
    processed_labels = []
    cell_expression = []

    for label in cell_labels:

        # Skip if background label
        if label == 0:
            continue

        # Calculate expression across all channels for the current cell
        cell_mask = (mask == label)
        cell_image = image[:, cell_mask]
        cell_mean_expression = np.mean(cell_image, axis=1)

        # Append the mean expression to the list
        cell_expression.append(cell_mean_expression)
        processed_labels.append(label)

    # Convert the list to a numpy array then a pandas dataframe
    cell_expression = np.array(cell_expression)
    cell_expression = pd.DataFrame(cell_expression)

    # Set generic marker names for the columns (not provided)
    cell_expression.columns = [f"Marker_{i}" for i in range(cell_expression.shape[1])]

    # Exclude the background label from the index
    cell_expression.index = [str(label) for label in processed_labels]
    return cell_expression


def normalize_expression(cell_expression, cofactors = 5.0):
    # Normalize the expression data using arcsinh transformation

    cell_expression = np.arcsinh(cell_expression / cofactors)
    return cell_expression


def quality_control(cell_expression):
    # Perform quality control on the cell expression data

    # Convert to AnnData object
    adata = sc.AnnData(cell_expression)

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata)
    
    # Filter cells with too few detected markers
    sc.pp.filter_genes(adata, min_cells = 3)

    # Filter markers detected in too few cells
    sc.pp.filter_cells(adata, min_genes = 100)

    # Filter cells with extreme total counts
    total_counts = adata.obs['total_counts']
    adata = adata[(total_counts > np.percentile(total_counts, 1)) &
                  (total_counts < np.percentile(total_counts, 99)), :]
    
    return adata

def compute_marker_correlation(adata):
    # Computes the correlation of each marker with the first 10 PCA components

    # Perform PCA with ten components
    sc.tl.pca(adata, n_comps = 10)
    components = adata.obsm['X_pca'][:, :10]

    # Calculate correlations between markers and each component
    correlations = pd.DataFrame(index = adata.var_names, 
                                columns = [f'PC{i+1}' for i in range(10)])
    for i, marker in enumerate(adata.var_names):
        for j in range(10):
            correlations.iloc[i, j] = np.corrcoef(adata.X[:, i], components[:, j])[0, 1]

    return correlations


def plot_correlation(correlations):
    # Plot a heatmap of the correlations

    plt.figure(figsize = (10, 8))
    sns.heatmap(correlations.astype(float), cmap = 'vlag', center = 0)
    plt.title('Correlation of Markers with First 10 PCA Components')
    plt.xlabel('PCA Component')
    plt.ylabel('Marker')
    plt.tight_layout()
    plt.show()


def cluster_cells(adata):
    # Perform clustering on the cell expression data

    # Reduce dimensionality
    sc.tl.pca(adata, n_comps = 20)

    # Construct a neighbourhood graph
    sc.pp.neighbors(adata, n_neighbors = 10, n_pcs = 20)

    # Cluster via the leiden algorithm
    sc.tl.leiden(adata, resolution = 0.5)

    # Compute UMAP coordinates
    sc.tl.umap(adata)

    return adata

def plot_umap(adata):
    # Plot UMAP results with clusters

    sc.pl.umap(adata, color = 'leiden', title = 'UMAP Clustering', frameon = False)

def compute_differential_expression(adata):
    # Perform differential expression analysis on the clustered data
    
    # Obtain cluster-specific differentially expressed markers
    sc.tl.rank_genes_groups(adata, groupby = 'leiden', method = 'wilcoxon')
    
    return adata

def plot_top_differentials(adata):
    # Plot the top five differentially expressed markers as a dotplot

    sc.pl.rank_genes_groups_dotplot(adata, groupby = 'leiden', standard_scape = 'var', n_genes = 5)



def main():

    # load the IMC image and mask
    imc_image = load_tiff_image("data/BaselTMA_SP41_23.475kx17.66ky_10000x5000_12_20170905_56_62_X7Y7_108_a0_full.tiff")
    imc_mask = load_tiff_image("data/BaselTMA_SP41_23.475kx17.66ky_10000x5000_12_20170905_56_62_X7Y7_108_a0_full_maks.tiff")

    # calculate the cell expression as cells x proteins 
    cell_expression = get_cell_expression(imc_image, imc_mask)

    # preprocess the data
    cell_expression = normalize_expression(cell_expression)
    adata = quality_control(cell_expression)

    # calculate marker correlations and visualize them
    correlation = compute_marker_correlation(adata)
    plot_correlation(correlation)

    # cluster the cells and visualize the UMAP results
    adata = cluster_cells(adata)

    # perform differential expression analysis on the clustered data
    adata = compute_differential_expression(adata)
    plot_top_differentials(adata)


if __name__ == "__main__": 
    main()