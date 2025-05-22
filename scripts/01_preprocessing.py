import numpy as np
import imctools
import tifffile
import os

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
