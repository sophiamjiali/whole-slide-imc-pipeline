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
import scipy as scp
import os

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


def create_foreground_mask(image, threshold):
    # Defines a foreground mask, filtering out blank acquisitions

    # sum intensity across all channels, 2D map of total channel

    # boolean mask of pixels with summed intensity above threshold

    # label foreground as 'cells', background is blank acquisition

    return

def save_generated_mask(mask, filepath):
    # Saves the generated foreground mask to a file

    return



def extract_pixel_intensities(image, mask):
    # Extracts pixel-level intensities from the image using the foreground mask

    return





def normalize_image(image, cofactors):
    # Normalize and scale the image via winsorization and Z-score

    image = 

    # to do
    return


def quality_control(image):
    # Performs quality control on the image data

    return


def filter_image(image, threshold):
    # Filters image based on threshold value for pixel intensities

    return



def create_anndata(image, mask):
    # Creates an AnnData object from the image and foreground mask

    return

def save_anndata(anndata, filepath):
    # Saves the AnnData object to a file

    return

