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

def bufferNHD(width_set, index, cpus_per_task):
    # Control flow
    if width_set == 'mean':
        width = 'WidthM'
    elif width_set == 'min':
        width = 'WidthM_Min'
    elif width_set == 'max':
        width = 'WidthM_Max'
    else:
        print('Invalid width option specified, exiting.')
        
    # Initialize parallel
    pandarallel.initialize(nb_workers=cpus_per_task)
    
    ## Prepare data
    # Read segmented NHD
    df, huc4, huc2 = readNHD(index=slurm, segmented=False)

    # Buffer flowlines with an extra 50 m on each side to be safe
    # This is beyond the max distance that the pixels could extend
    # once converted to pseudo pixels
    df['buffers'] = df.parallel_apply(user_defined_function=specialBuffer,
                                                             args=(width,
                                                                   # args are
                                                                   # cap_style, segmented, extra
                                                                   'flat', False, False),
                                                             axis=1)
        
    # Drop original reach geometry column, set buffered geometry as active geometry
    df = df.drop(columns='geometry').set_geometry('buffers').set_crs(crs=df.crs)

        # Write out
    # Set write filepath
    save_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_buffered/'
    save_path = os.path.join(save_path, huc2)
    save_file = huc4 + '_prepped_buffered_' + width_set + '.parquet'

    #Write out gdf as parquet file
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    df.to_parquet(os.path.join(save_path, save_file))
    
    print('Script completed, wrote out results.')
    
### PARSE ARGUMENTS
parser = ArgumentParser(description='Please specify whether you would\
                        like to use the min, mean, or max predicted\
                        bankfull width for this analysis.')

parser.add_argument('width_set', type=str, help='min, mean, or max')
args = parser.parse_args()
width_set = args.width_set

print('width_set: ' + width_set)

### GET JOB INDEX, CPUS
slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])
cpus = int(os.environ.get('SLURM_CPUS_PER_TASK'))
cpus_per_task = cpus if cpus < 65 else 1

if __name__ == "__main__":
    bufferNHD(width_set=width_set, index=slurm, cpus_per_task=cpus_per_task)