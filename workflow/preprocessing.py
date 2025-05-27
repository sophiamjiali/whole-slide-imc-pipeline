"""
Script:          preprocessing.py
Purpose:         Preprocesses WS-IMC image data for downstream analysis
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            05-26-2025
"""

import tifffile as tf
import numpy as np
import pandas as pd
import scanpy as sc
import os

from scipy.stats.mstats import winsorize
from sklearn.preprocessing import MinMaxScaler


def load_imc_image(filepath):
    # Loads an IMC image from a TIFF file

    try:
        image = tf.imread(filepath)
    except Exception as e:
        file_name = os.path.basename(filepath)
        print(f"Error loading image {file_name}: {e}")
        return None

    return image


def load_marker_mapping(filepath):
    # Loads marker mappings from a CSV file



    return


def create_mask(image, threshold):
    # Defines a foreground mask, filtering out blank acquisitions

    # Sum intensities across all channels to create a 2D map
    channel_sums = image.sum(axis = 0)

    # Create a binary mask for foreground pixels based on the threshold
    foreground_mask = channel_sums > threshold
    foreground_mask = np.where(foreground_mask, 1, 0)

    # Assign unique labels to each foreground pixel
    labeled_mask = np.zeros_like(foreground_mask, dtype = int)
    labeled_mask[foreground_mask == 1] = np.arange(1, np.sum(foreground_mask == 1) + 1)

    return labeled_mask

def save_generated_mask(mask, filepath):
    # Saves the generated foreground mask to a file

    return


def extract_pixel_intensities(image, mask, marker_names):
    # Extracts pixel-level intensities from the image using the foreground mask

    # Extract intensities for each labeled pixel in the mask
    pixel_labels = np.unique(mask)
    pixel_expression = []

    for label in pixel_labels:
        
        # Skip if background label
        if label == 0:
            continue

        # Append pixel intensities for the current label
        curr_expression = image[:, mask == label]
        pixel_expression.append(curr_expression)

    # Conver the list to a NumPy array, then to a Pandas DataFrame
    pixel_expression = np.array(pixel_expression)
    pixel_expression = pd.DataFrame(pixel_expression)

    # Annotate the DataFrame with marker names
    pixel_expression.columns = marker_names

    # Exclude the background label from the index
    pixel_expression.index = [str(label) for label in pixel_labels if label != 0]

    return pixel_expression


def normalize_intensities(pixel_expression, quantile):
    # Normalize and scale the image via winsorization and min-max scaling

    winsorized_expr = winsorize(pixel_expression, limits = [quantile, quantile])
    scaler = MinMaxScaler()
    normalized_expr = scaler.fit_transform(winsorized_expr)

    return normalized_expr

def create_anndata(pixel_intensities):
    # Creates an AnnData object from the image and foreground mask

    adata = sc.AnnData(pixel_intensities)
    return adata


def quality_control(adata):
    # Performs quality control on the image data

    # applicable if expression is pixel-level intensities?

    return


def filter_intensities(adata, threshold):
    # Filters image based on threshold value for pixel intensities

    return


def save_anndata(anndata, filepath):
    # Saves the AnnData object to a file

    return

