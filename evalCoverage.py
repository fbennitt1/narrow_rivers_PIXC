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

def evalCoverage(width_set, index, cpus_per_task, huc2, data_path, save_dir):
    ## SET UP
    if width_set == 'mean':
        width = 'WidthM'
        binn = 'Bin'
    elif width_set == 'min':
        width = 'WidthM_Min'
        binn = 'Bin_Min'
    elif width_set == 'max':
        width = 'WidthM_Max'
        binn = 'Bin_Max'
    else:
        print('Invalid width option specified, exiting.')

    ### PIXEL CLOUD
    # Get PIXC filepath
    file_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/PIXC_v2_0_HUC2_' + huc2 + '_filtered.json'
    data = open_json(file_path)
    file_name = data[index]
    
    # Get PIXC metadata
    granule_name = file_name[:-3]
    tile_name = file_name[20:28]
    pass_num = int(file_name[20:23])

    print(granule_name)

    ## Read in PIXC
    # PIXC datapath
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
                 'prior_water_prob', 'classification']

    # Convert PIXC to GeoDataFrame
    gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_vars=variables)

    ### NHDPlus HR
    ## Find correct HUC4s
    # Read in tile and HUC4 intersection data
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    dtype_dic= {'tile': str, 'huc4': str, 'coverage': float}
    tile_huc4 = pd.read_csv(os.path.join(mdata_path,
                                        'huc4_swot_science_tiles.csv'),
                            dtype=dtype_dic)
    # Make list of HUC4s that intersect the tile
    hucs = list(tile_huc4[tile_huc4['tile'] == tile_name]['huc4'])
    # Limit to the current HUC2
    hucs = [x for x in hucs if x.startswith(huc2)]

    ## Get NHD index metadata
    # Define dtypes for lookup tables to preserve leading zeros
    dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
    # Read in HUC lookup table
    huc_lookup = pd.read_csv(os.path.join(mdata_path,
                                      'HUC4_lookup_no_great_lakes.csv'),
                         dtype=dtype_dic)
    # Extract indices for read-in
    indices = list(huc_lookup[huc_lookup['HUC4'].isin(hucs)]['slurm_index'])

    ## READ IN HUC4 FLOWLINES
    # Create merged dataframe of all flowlines intersected
    if len(indices) == 1:
        # Read prepped NHD
        flowlines, _, _,  = readNHD(index=indices[0])
        # huc4_list, huc2_list = readNHD(index=indices[0])
    else:
        # Initialize lists
        d = []
        # huc4_list = []
        # huc2_list = []
        # Loop through indices and store in lists
        for i in indices:
            # Read prepped NHD
            flowlines, _, _ = readNHD(index=i)
            # huc4, huc2 = readNHD(index=i)
            # Append to lists
            d.append(flowlines)
            # huc4_list.append(huc4) # I DON'T DO ANYTHING WITH THIS
            # huc2_list.append(huc2) # I DON'T DO ANYTHING WITH THIS
        # Merge GeoDataFrames
        flowlines = pd.concat(d)

    # Project CRS (currently to WGS 84 / UTM zone 18N)
    flowlines = flowlines.to_crs(epsg=3857)

    # Initialize panarallel
    pandarallel.initialize(nb_workers=cpus_per_task)

    # Buffer flowlines with an extra 50 m on each side to be safe
    # This is beyond the max distance that the pixels could
    # extend once converted to pseudo pixels
    flowlines['buffer'] = flowlines.parallel_apply(user_defined_function=specialBuffer,
                                                             args=(width, 'flat', False, True),
                                                             axis=1)          
    # Set geometry to buffered reaches
    flowlines = flowlines.set_geometry('buffer').set_crs(epsg=3857)

    ## Clip masked pixels to buffered reaches
    gdf_PIXC_clip = gpd.sjoin(gdf_PIXC, flowlines, how='inner', predicate='within')
    # Check that there are pixels that intersect
    if gdf_PIXC_clip.shape[0] == 0:
        print('This granule has no pixels that intersect reaches, exiting.')
        sys.exit()
    # Drop unneeded cols
    gdf_PIXC_clip = gdf_PIXC_clip.drop(columns=['index_right', 'NHDPlusID',
                                                'GNIS_Name', 'LengthKM',
                                                'WidthM', 'WidthM_Min',
                                                'WidthM_Max', 'Bin', 'Bin_Min',
                                                'Bin_Max', 'geometry_right'])

    ### NADIR TRACK
    # Get single pixel for selecting correct nadir segment
    pixel_pt = gdf_PIXC_clip.iloc[0].geometry      
    # Find correct nadir segment and return its geometry
    nadir_segment_ln = findNadir(pass_num=pass_num, pixel_pt=pixel_pt)

    ### MAKE PSEUDO-PIXELS
    # Set along-track pixel resolution
    azimuth_res = 22 # meters
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

    print('Indices: ' + str(indices))

    ### READ IN SEGMENTS
    # Create merged dataframe of all basins intersected
    if len(indices) == 1:
        # Read prepped NHD
        segments, _, _ = readNHD(index=indices[0], segmented=True) # DO I WANT TO ALSO EXTRACT HUC4
    else:
        # Initialize lists
        d = []
        # Loop through indices and store in lists
        for i in indices:
            # Read prepped NHD
            segments, huc4, _ = readNHD(index=i, segmented=True)
            # Make column with HUC4 id
            segments['huc4_long'] = huc4
            segments['huc4'] = segments['huc4_long'].str[10:14] # DO I USE THESE
            # Rename segments to geometry
            # segments = segments.rename(columns={'segments': 'geometry'}).set_geometry('geometry')
            # Append to list
            d.append(segments)
        # Merge GeoDataFrames
        segments = pd.concat(d)

    # print(segments.columns)
    # sys.exit()

    # Project CRS (currently to WGS 84 / UTM zone 18N)
    segments = segments.to_crs(epsg='3857')

    ## Clean-up
    segments = segments.reset_index().rename(columns={'index': 'index_old'})
    # Assign a unique counter within each index group
    segments['counter'] = segments.groupby('NHDPlusID').cumcount()    
    # Keep only first ten segments (some reaches repeat)
    segments = segments[segments['counter'] < 10]
    # Clip the segments to the bounds of the PIXC with pseudo-pixels           
    segments = segments.clip(pseudo_bounds)
    # Keep only reaches that are fully contained in PIXC granule
    segments = segments.groupby('NHDPlusID').filter(lambda x: len(x) == 10)
    
    # Get number of reaches per bin
    counts = pd.DataFrame(segments[binn].value_counts()).reset_index()

    ## Buffer segments
    segments['buffer'] = segments.parallel_apply(user_defined_function=specialBuffer, args=(width,'flat', True, False), axis=1)        
    # Set active geometry col to buffered segments
    segments = segments.set_geometry('buffer')                   
    # Calculate segment area
    segments['segment_area'] = segments.geometry.area

    ## Merge the segments and pseudo-puxels by intersection
    sj = gpd.sjoin(segments, gdf_PIXC_clip, how='left', predicate='intersects')
    # Drop unneeded columns
    sj = sj.drop(columns=['index_right', 'points', 'azimuth_index',
                          'range_index',
                          'height', 'geoid',
                          'klass', 'latitude', 'longitude'])
    # Set active geometry column for dissolve
    sj = sj.set_geometry('pseudo_geom')

    ## Dissolve
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
    # Fill zeros for stats
    sj['coverage'] = sj['coverage'].fillna(0)

    ### DO STATS
    bins = sj.Bin.unique()
    
    # reaches_cent, reaches_min = summarizeCoverage(df=sj, binn=binn,
    #                                                 bins=bins, counts=counts)
    reaches_cent, reaches_thresh, reaches_min = summarizeCoverage(df=sj, binn=binn,
                                                bins=bins, counts=counts)

    ### WRITE OUT
    save_path = os.path.join('/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_output/',
                             save_dir)

    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    # sj.to_csv(os.path.join(save_path, granule_name + '_coverage.csv'))

    ### MAKE PARQUET
    reaches_cent.to_parquet(os.path.join(save_path, granule_name + '_reaches_cent.parquet'))
    reaches_thresh.to_parquet(os.path.join(save_path, granule_name + '_reaches_thresh.parquet'))
    reaches_min.to_parquet(os.path.join(save_path, granule_name + '_reaches_min.parquet'))
    print('Script completed, wrote out results.')
    
### PARSE ARGUMENTS
parser = ArgumentParser(description='Please specify whether you would\
                        like to use the min, mean, or max predicted\
                        bankfull width for this analysis.')
parser.add_argument('width_set', type=str, help='min, mean, or max')
args=parser.parse_args()
width_set = args.width_set

### GET JOB INDEX, CPUS
slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])
cpus = int(os.environ.get('SLURM_CPUS_PER_TASK'))
cpus_per_task = cpus if cpus < 65 else 1

if __name__ == "__main__":
    huc2 = '15'
    data_path = '/nas/cee-water/cjgleason/fiona/data/PIXC_v2_0_HUC2_' + huc2
    # pixc_ref = 'PIXC_v2_0_HUC2_01_best_files_no_exits.csv' ## CHANGE THIS
    save_dir = 'PIXC_v2_0_HUC2_' + huc2 + '_2025_03_06_'+ width_set
    
    evalCoverage(width_set=width_set, index=slurm, cpus_per_task=cpus_per_task, huc2=huc2, data_path=data_path, save_dir=save_dir)