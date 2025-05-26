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

def load_imc_tiff(filepath):

    # to do
    return

def normalize_image(image, cofactors = 5.0):

    # to do
    return

def sliding_window(image, channel_names, window_size = 200, step = 200):
    """
    Calculate marker correlations in sliding windows over a whole-slide IMC image.
    """
