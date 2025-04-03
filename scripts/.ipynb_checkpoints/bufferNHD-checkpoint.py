from argparse import ArgumentParser

import os
import sys
import time

import matplotlib.pyplot as plt

import geopandas as gpd
import pandas as pd

from pandarallel import pandarallel

from reaches import readNHD
from utils import specialBuffer

def bufferNHD(width_set, index, cpus_per_task, data_path, save_dir):
    
    
    
    
    
    
### PARSE ARGUMENTS
parser = ArgumentParser(description='Please specify whether you would\
                        like to use the min, mean, or max predicted\
                        bankfull width for this analysis.')
parser.add_argument('segmented', type=bool, help='True or False')
parser.add_argument('width_set', type=str, help='min, mean, or max')
args=parser.parse_args()
segmented = arge.segmented
width_set = args.width_set

### GET JOB INDEX, CPUS
slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])
cpus = int(os.environ.get('SLURM_CPUS_PER_TASK'))
cpus_per_task = cpus if cpus < 65 else 1

if __name__ == "__main__":

    save_dir = '
    
    bufferNHD(width_set=width_set, index=slurm, cpus_per_task=cpus_per_task, huc2=huc2, data_path=data_path, save_dir=save_dir)