'''
To create the files needed for this, first run bufferFlowlinesNLCD.py using
runBufferFlowNLCD.sh to generate the buffered NHD flowlines as GeoJSON with
the correct CRS.
'''

import os

import pandas as pd

from argparse import ArgumentParser
from exactextract import exact_extract

def intersectNLCD(width_set, index):
    # Control flow
    if width_set == 'mean':
        width = 'WidthM'
    elif width_set == 'min':
        width = 'WidthM_Min'
    elif width_set == 'max':
        width = 'WidthM_Max'
    else:
        print('Invalid width option specified, exiting.')

    nlcd_path = '/nas/cee-water/cjgleason/data/NLCD/Annual_NLCD_LndCov_2023_CU_C1V0.tif'
    
    ## Set-up
    prep_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_buffered_json/'
    
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    # Define dtypes for lookup tables to preserve leading zeros
    dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
    # Read in HUC lookup table
    lookup = pd.read_csv(os.path.join(mdata_path,
                                      'HUC4_lookup_no_great_lakes.csv'),
                         dtype=dtype_dic)

    # Get current HUC2 and HUC4 IDs
    huc2 = 'HUC2_' + lookup.loc[index,'HUC4'][0:2]
    huc4 = 'NHDPLUS_H_' + lookup.loc[index,'HUC4'] + '_HU4_GDB'

    file_path = os.path.join(prep_path, huc2, huc4 + \
                             '_prepped_buffered_' + width_set + '.json')
    
    ## Intersect
    # IGNORE THE MISMATCHED SPATIAL REF WARNING, IT IS NOT TRUE IN
    # THIS CASE. Tested in prepNLCD.ipynb, we are good.
    df = exact_extract(rast=nlcd_path, vec=file_path,
                       ops=['frac', 'variety', 'mode', 'minority'],
                       include_cols=["NHDPlusID"], output='pandas')
   
    ## Write out
    # Set write filepath
    save_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_NLCD/'
    save_path = os.path.join(save_path, huc2)
    save_file = huc4 + '_NLCD_intersected_' + width_set + '.parquet'

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

if __name__ == "__main__":
    intersectNLCD(width_set=width_set, index=slurm)