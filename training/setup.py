"""
Script:          feature_extraction.py
Purpose:         Extracts sliding window features from WS-IMC images
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            05-26-2025

PyTorch Version: 2.7.1
"""

# ========== Imports ==========

## Standard Libraries
import os
import numpy as np
import random
from PIL import Image
from types import SimpleNamespace

## Imports for Plotting
import matplotlib.pyplot as plt
from IPython.display import set_matplotlib_formats
import matplotlib as mpl
import seaborn as sns

set_matplotlib_formats('svg', 'pdf') # Set default plotting style
mpl.rcParams['lines.linewidth'] = 2 # Set default line width for plots
sns.reset_orig()

## PyTorch Libraries
import torch
import torch.nn as nn
import torch.utils.data as data
import torch.optim as optim

import torchvision
from torchvision import transforms

import numpy as np

# ========== Configurations ==========

DATASET_PATH = 'data/processed'

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
set_seed(42)  

torch.backends.cudnn.deterministic = True
torch.backends.cudnn.benchmark = False

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")