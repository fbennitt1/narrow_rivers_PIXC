'''
A script to write out GeoJSONs of buffered NHD reaches in the crs of the NLCD.

Input are the the netcdf files which contain gridded daily forcings. (Coordinates are "Lat", "Lon")

The output are netcdf files which contain zonal daily forcings for each lake (Coordinates are "lake_id" and "date")

Input: 1980_01.nc, 1981_02.nc, 1983_03.nc, ..., 2022_12.nc
Output: zonal_1980_01.nc, zonal_1981_02.nc, ..., zonal_2022_12.nc

'''
from argparse import ArgumentParser

import os
import sys
import time

import matplotlib.pyplot as plt

import geopandas as gpd
import pandas as pd
import rasterio

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
    # print('Original crs: ' + str(df.crs))

    # Buffer flowlines with an extra 32 m on each side to be safe
    # This is just beyond the max distance that the pseudo-pixels
    # could extend beyond the pixel centroid.
    df['buffers'] = df.parallel_apply(user_defined_function=specialBuffer,
                                                             args=(width,
                                                                   # args are
                                                                   # cap_style, segmented, extra
                                                                   'flat', False, True),
                                                             axis=1)
        
    # Drop original reach geometry column, set buffered geometry as active geometry
    df = df.drop(columns='geometry').set_geometry('buffers').set_crs(crs=df.crs)
    df = df[['NHDPlusID', 'buffers']]
    
#     # Open NLCD to get proper crs
#     nlcd_path = '/nas/cee-water/cjgleason/data/NLCD/Annual_NLCD_LndCov_2023_CU_C1V0.tif'
    
#     with rasterio.open(nlcd_path) as src:
#         nlcd_crs = src.crs

    # Reproject data to match NLCD
    # df = df.to_crs(nlcd_crs)
    # print('New crs: ' + str(df.crs))

    ## Write out
    # Set write filepath
    save_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_buffered_json/'
    save_path = os.path.join(save_path, huc2)
    save_file = huc4 + '_prepped_buffered_' + width_set + '.json'

    #Write out gdf as parquet file
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    df.to_file(filename=os.path.join(save_path, save_file), driver='GeoJSON')
    
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