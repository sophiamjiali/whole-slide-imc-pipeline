"""
Script:          run_pipeline.py
Purpose:         Runs the WS-IMC analysis pipeline
Author:          Sophia Li
Affiliation:     Campbell Lab
Date:            05-26-2025
"""

from utils.setup import load_config, initialize_directories
from workflow.preprocessing import load_imc_tiff, normalize_image, sliding_window



# set up environment with setup.py
# convert MCD files to TIFF with convert_mcd_to_tiff.py
# start workflow