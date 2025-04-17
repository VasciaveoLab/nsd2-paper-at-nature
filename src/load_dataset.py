import numpy as np
import pandas as pd

import scanpy as sc
sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

import scvelo as scv

import matplotlib.pyplot as plt

import os
import getpass
import logging
from datetime import datetime

## These lines are to fix saving PDF with text recognizible by Adobe Illustrator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="matplotlib")
warnings.filterwarnings("ignore", category=UserWarning, module="scanpy")

my_random_seed = 666
n_cells_to_subsample = 3e3
n_cores = 8

### --- MAIN --- ###
def main():
    
    # Set the experiment name
    experiment_name = "nsd2-paper-experiment"
    workspace_absolute_path = "/shares/vasciaveo_lab/programs/nepc_organoids"

    date_time = datetime.now().strftime("%b_%d_%Y_h%Hm%M")
    print("date and time:",date_time)
    experiments_dir = experiment_name + "_" + date_time
    # experiments_path = os.path.join(experiments_dir, "figures")
    experiments_path = os.path.join(workspace_absolute_path, experiments_dir )
    os.makedirs(experiments_path,exist_ok=True)

if __name__ == '__main__':
    main()
