from argparse import ArgumentParser
import os
import sys
import time

import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr

from pandarallel import pandarallel
# from shapely.geometry import box

from reaches import *
from utils import *

### PARSE ARGUMENTS
parser = ArgumentParser(description='Please specify whether you would\
                        like to use the min, mean, or max predicted\
                        bankfull width for this analysis.')
parser.add_argument('width_set', type=str, help='min, mean, or max')
args=parser.parse_args()
width_set = args.width_set

# # FOR NOW, SET
# width = 'WidthM'

# Set correct width
if width_set == 'mean':
    width = 'WidthM'
elif width_set == 'min':
    width = 'WidthM_Min'
elif width_set == 'max':
    width = 'WidthM_Max'
else:
    print('Invalid width option specified, exiting.')
    sys.exit()


### PIXEL CLOUD
# Get PIXC index metadata
mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
dtype_dic= {'cycle': str, 'pass': str, 'tile': str, 'version': str}

# Read in HUC lookup table
pixc_lookup = pd.read_csv(os.path.join(mdata_path,
                                       'PIXC_v2_0_HUC2_01_best_files_no_exits.csv'),
                          dtype=dtype_dic).drop(columns='index')

# Get job index
# slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])
slurm = 0

# Get filepath and other things
file_name = pixc_lookup.loc[slurm, 'files']
granule_name = file_name[:-3]
tile_name = file_name[20:28]
pass_num = int(file_name[20:23]) 

## Read in PIXC
# PIXC datapath
data_path = '/nas/cee-water/cjgleason/fiona/data/PIXC_v2_0_HUC2_01/'
pixc_path = os.path.join(data_path, file_name)

# Read in pixel group
ds_PIXC = xr.open_mfdataset(paths=pixc_path, group = 'pixel_cloud',
                            engine='h5netcdf')

# Make mask
mask = bitwiseMask(ds_PIXC)

# Check that we didn't lose all pixels
if mask.shape[0] == 0:
    print('This granule has no pixels after masking, exiting.')
    sys.exit()

# Set desired data vars
variables = ['azimuth_index', 'range_index', 'cross_track',
             'pixel_area', 'height', 'geoid',
             'dlatitude_dphase', 'dlongitude_dphase',
             'dheight_dphase', 'classification']
                   
# Convert PIXC to GeoDataFrame
gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_vars=variables)

### NHDPlus HR
## Find correct HUC4s
# Read in tile and HUC4 intersection data
dtype_dic= {'tile': str, 'huc4': str, 'coverage': float}
tile_huc4 = pd.read_csv(os.path.join(mdata_path,
                                    'huc4_swot_science_tiles.csv'),
                        dtype=dtype_dic)
# Make list of HUC4s that intersect our tile
hucs = list(tile_huc4[tile_huc4['tile'] == tile_name]['huc4'])

## Get NHD index metadata
# Define dtypes for lookup tables to preserve leading zeros
dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
# Read in HUC lookup table
huc_lookup = pd.read_csv(os.path.join(mdata_path,
                                  'HUC4_lookup_no_great_lakes.csv'),
                     dtype=dtype_dic)
# Extract indices for read-in
indices = list(huc_lookup[huc_lookup['HUC4'].isin(hucs)]['slurm_index'])
                 
## READ IN HUC4 BASINS
# Create merged dataframe of all basins intersected
if len(indices) == 1:
    # Read prepped NHD
    flowlines, huc4_list, huc2_list = readNHD(index=indices[0])
else:
    # Initialize lists
    d = []
    huc4_list = []
    huc2_list = []
    # Loop through indices and store in lists
    for i in indices:
        # Read prepped NHD
        flowlines, huc4, huc2 = readNHD(index=i)
        # Append to lists
        d.append(flowlines)
        huc4_list.append(huc4) # I DON'T DO ANYTHING WITH THIS
        huc2_list.append(huc2) # I DON'T DO ANYTHING WITH THIS
    # Merge GeoDataFrames
    flowlines = pd.concat(d)
    
# Project CRS (currently to WGS 84 / UTM zone 18N)
flowlines = flowlines.to_crs(epsg=32618)

# Initialize panarallel
pandarallel.initialize()
                   
# Buffer flowlines with an extra 50 m on each side to be safe
# This is beyond the max distance that the pixels could
# extend once converted to pseudo pixels
flowlines['buffer'] = flowlines.parallel_apply(user_defined_function=specialBuffer,
                                                         args=(width,
                                                               'flat',True),
                                                         axis=1)          
# Set geometry to buffered reaches
flowlines = flowlines.set_geometry('buffer').set_crs(epsg=32618)

## Clip masked pixels to buffered reaches
gdf_PIXC_clip = gpd.sjoin(gdf_PIXC, flowlines, how='inner', predicate='within')
# Check that there are pixels that intersect
if gdf_PIXC_clip.shape[0] == 0:
    print('This granule has no pixels that intersect reaches, exiting.')
    sys.exit() 
# Drop unneeded cols
gdf_PIXC_clip = gdf_PIXC_clip.drop(columns=['index_right',
                                            'Bin', 'GNIS_Name',
                                            'LengthKM', 'NHDPlusID',
                                            'WidthM', 'geometry_right'])
                   
### NADIR TRACK
# Get single pixel for selecting correct nadir segment
pixel_pt = gdf_PIXC_clip.iloc[0].geometry      
# Find correct nadir segment and return its geometry
nadir_segment_ln = findNadir(pass_num=pass_num, pixel_pt=pixel_pt)
                   
### MAKE PSEUDO-PIXELS
# Set along-track pixel resolution
azimuth_res = 21 # meters
# Make pseudo pixels
gdf_PIXC_clip['pseudo_pixel'] = gdf_PIXC_clip.parallel_apply(user_defined_function=makePseudoPixels,
                                                         args=(nadir_segment_ln,
                                                               azimuth_res),
                                                         axis=1)
# Clean-up
gdf_PIXC_clip = gdf_PIXC_clip.rename(columns={'geometry': 'pixel_centroid'}).set_geometry('pseudo_pixel')
#Get bounds of PIXC tile
pseudo_bounds = gdf_PIXC_clip.total_bounds
# Copy geometry column as sjoin will discard it
gdf_PIXC_clip['pseudo_geom'] = gdf_PIXC_clip.geometry
                   
### READ IN SEGMENTS
# Create merged dataframe of all basins intersected
if len(indices) == 1:
    # Read prepped NHD
    segments, _, _ = readNHD(index=indices[0]) # DO I WANT TO ALSO EXTRACT HUC4
else:
    # Initialize lists
    d = []
    # Loop through indices and store in lists
    for i in indices:
        # Read prepped NHD
        segments, huc4, _ = readSegments(index=i)
        # Make column with HUC4 id
        segments['huc4_long'] = huc4
        segments['huc4'] = segments['huc4_long'].str[10:14] # DO I USE THESE
        # Rename segments to geometry
        segments = segments.rename(columns={'segments': 'geometry'}).set_geometry('geometry')
        # Append to list
        d.append(segments)
    # Merge GeoDataFrames
    segments = pd.concat(d)

# Project CRS (currently to WGS 84 / UTM zone 18N)
segments = segments.to_crs(epsg='32618')
# Clean-up
segments = segments.reset_index().rename(columns={'index': 'index_old'})
print("Shape of segments after reset: " + str(segments.shape))

# Assign a unique counter within each index group
segments['counter'] = segments.groupby('NHDPlusID').cumcount()    
print("Shape of segments after counter: " + str(segments.shape))

# Keep only first ten segments (some reaches repeat)
segments = segments[segments['counter'] < 10]
print("Shape of segments after dropping extra: " + str(segments.shape))

# Clip the segments to the bounds of the PIXC with pseudo-pixels           
segments = segments.clip(pseudo_bounds)
print("Shape of segments after clipping: " + str(segments.shape))

# Keep only reaches that are fully contained in PIXC granule
segments = segments.groupby('NHDPlusID').filter(lambda x: len(x) == 10)
print("Shape of segments after groupby: " + str(segments.shape))
                   
# Buffer segments
segments['buffer'] = segments.parallel_apply(user_defined_function=specialBuffer, args=(width,'flat', False), axis=1)        
# Set active geometry col to buffered segments
segments = segments.set_geometry('buffer')                   
# Calculate segment area
segments['segment_area'] = segments.geometry.area

# Merge the segments and pseudo-puxels by intersection
sj = gpd.sjoin(segments, gdf_PIXC_clip, how='left', predicate='intersects')
# Drop unneeded columns
sj = sj.drop(columns=['index_right', 'points', 'azimuth_index',
                      'range_index', 'cross_track', 'pixel_area',
                      'height', 'geoid', 'dlatitude_dphase',
                      'dlongitude_dphase', 'dheight_dphase',
                      'klass', 'latitude', 'longitude', ])
# Set active geometry column for dissolve
sj = sj.set_geometry('pseudo_geom')

# Dissolve
sj = sj.groupby('NHDPlusID', as_index=False).parallel_apply(user_defined_function=specialDissolve)
# Drop multi-index
sj = sj.reset_index().drop(columns=['level_0', 'level_1'])
# Clip dissolved pseudo-pixels to node areas
sj['pseudo_geom_clip'] = sj.parallel_apply(user_defined_function=specialClip,
                                                         axis=1)
# Calculate the pseudo-pixel area within each node
sj['pseudo_area'] = sj.pseudo_geom_clip.area
# Calculate coverage
sj['coverage'] = sj.pseudo_area/sj.segment_area
# Copy to fill zeros for stats
sj_w_zero = sj.copy()
sj_w_zero['coverage'] = sj_w_zero['coverage'].fillna(0)

### DO STATS
bins = sj.Bin.unique()
## Nodes
# Get descriptive stats
node_desc = sj.groupby('Bin')['coverage'].describe().reset_index()
node_desc['with_zero'] = 0

node_quant = pd.DataFrame(sj.groupby('Bin')['coverage'].quantile(q=[x / 100.0 for x in range(0,100,1)])).reset_index().rename(columns={'level_1': 'quantile'})
node_quant['with_zero'] = 0

## Nodes with zeros
node_desc_w_zero = sj_w_zero.groupby('Bin')['coverage'].describe().reset_index()
node_desc_w_zero['with_zero'] = 1

node_quant_w_zero = pd.DataFrame(sj_w_zero.groupby('Bin')['coverage'].quantile(q=[x / 100.0 for x in range(0,100,1)])).reset_index().rename(columns={'level_1': 'quantile'})
node_quant_w_zero['with_zero'] = 1

# nodes_mean = sj.groupby('Bin')['coverage'].mean().to_list()
# nodes_std = sj.groupby('Bin')['coverage'].std().to_list()
# d = {'mean': nodes_mean, 'std': nodes_std}
# nodes = pd.DataFrame(data=d).T
# nodes.columns = bins

## Reaches
d = {}
for i in range(1, 10):
    threshold = i/10
    # print(threshold)
    
    detected = sj.groupby(['Bin', 'NHDPlusID'])['coverage'].apply(lambda x: (x > threshold).sum()) / 10
    reach = detected.groupby('Bin').mean().to_list()
    #quantiles
    #std dev
    #n
    # size of the dataframe (are we including the reaches with nothing detected?????)
    
    d[threshold] = reach
    
reaches = pd.DataFrame(data=d).T
reaches.columns = bins
                   
### WRITE OUT
save_path = os.path.join('/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/', 'PIXC_v2_0_HUC2_01_testing')

if not os.path.isdir(save_path):
    os.makedirs(save_path)
    
# sj.to_csv(os.path.join(save_path, granule_name + '_coverage.csv'))
nodes.to_csv(os.path.join(save_path, granule_name + '_nodes.csv'))
reaches.to_csv(os.path.join(save_path, granule_name + '_reaches.csv'))