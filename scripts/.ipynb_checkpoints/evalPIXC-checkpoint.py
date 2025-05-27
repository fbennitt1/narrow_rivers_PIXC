from argparse import ArgumentParser
import os
import sys
import time

import dask.dataframe as dd
import dask_geopandas
import geopandas as gpd
import numpy as np
import pandas as pd
import shapely
import xarray as xr

from pandarallel import pandarallel

from reaches import *
from utils import *

def evalCoverage(width_set, suffix, index, cpus_per_task, huc2, data_path, save_dir):
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
        
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'

    ### PIXEL CLOUD
    # Get PIXC filepath
    file_path = os.path.join(mdata_path, 'PIXC_v2_0_HUC2_' + \
                             huc2 + '_filtered' + suffix + '.json') 
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
    ds_PIXC = xr.open_mfdataset(paths=pixc_path,
                                group = 'pixel_cloud',
                                engine='h5netcdf')

    # Make mask
    mask = bitwiseMask(ds_PIXC)

    # Check that we didn't lose all pixels
    if mask.shape[0] == 0:
        print('This granule has no pixels after masking, exiting.')
        sys.exit()

    # Set desired data vars
    variables = ['cross_track', 'water_frac',
                 'pixel_area', 'height', 'geoid', 'solid_earth_tide',
                 'load_tide_fes', 'pole_tide', 'prior_water_prob',
                 'classification_qual', 'geolocation_qual',
                 'sig0_qual', 'classification']

    # Convert PIXC to GeoDataFrame
    gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_vars=variables)
    
    del ds_PIXC
    
    ## Flag as-in RiverSP ## HERE
    gdf_PIXC['geo_qual_wse_suspect'] = 0
    gdf_PIXC['geo_qual_wse_suspect'] = np.where((gdf_PIXC['geolocation_qual'] >=2**7),
                                                1, gdf_PIXC['geo_qual_wse_suspect'])

    gdf_PIXC['class_qual_area_suspect'] = 0
    gdf_PIXC['class_qual_area_suspect'] = np.where((gdf_PIXC['classification_qual'] >=2**7),
                                                   1, gdf_PIXC['geo_qual_wse_suspect'])

    gdf_PIXC = gdf_PIXC.drop(columns=['geolocation_qual', 'classification_qual',
                                      'sig0_qual'])
    
    ## Make indicator for wse aggregation
    # (don't use lane_near_water or dark_water)
    wse_klass = [3.0, 4.0, 6.0, 7.0]

    gdf_PIXC['agg_wse'] = 0
    gdf_PIXC.loc[gdf_PIXC['klass'].isin(wse_klass), 'agg_wse'] = 1
    
    # Calculate wse
    gdf_PIXC['wse'] = gdf_PIXC.height - gdf_PIXC.geoid - \
                      gdf_PIXC.solid_earth_tide - gdf_PIXC.load_tide_fes -\
                      gdf_PIXC.pole_tide
    gdf_PIXC = gdf_PIXC.drop(columns=['height', 'geoid', 'solid_earth_tide',
                                      'load_tide_fes', 'pole_tide'])
    
    ### NHDPlus HR
    ## Find correct HUC4s
    # Read in tile and HUC4 intersection data
    dtype_dic= {'tile': str, 'huc4': str, 'coverage': float}
    tile_huc4 = pd.read_csv(os.path.join(mdata_path,
                                        'huc4_swot_science_tiles.csv'),
                            dtype=dtype_dic)
    
    # Make list of HUC4s that intersect the tile
    huc4s = list(tile_huc4[tile_huc4['tile'] == tile_name]['huc4'])
    # Limit to the current HUC2
    huc4s = [x for x in huc4s if x.startswith(huc2)]


    ## READ IN HUC4 BUFFERED FLOWLINES
    # (with extra 32 m on each side to capture full pseudo-pixels)
    # Create merged dataframe of all polygons

    data_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_buffered/HUC2_' + \
                 huc2 + '/'

    file_paths = []

    for huc in huc4s:
        file_path = data_path + 'NHDPLUS_H_' + huc + '_HU4_GDB_prepped_buffered_' + width_set + '.parquet'
        file_paths.append(file_path)
        
    reach_mask = dask_geopandas.read_parquet(path=file_paths,
                                             columns=['NHDPlusID', 'Slope',
                                                      'LengthKM', 'buffers'])
    
    reach_mask = reach_mask.compute()
    
    # Clip masked pixels to buffered reaches with extra width
    gdf_PIXC = gpd.sjoin(gdf_PIXC, reach_mask, how='inner',
                         predicate='within').reset_index().drop(columns=['index',
                                                                         'index_right'])
    if gdf_PIXC.shape[0] == 0:
        print('This granule has no pixels that intersect reaches, exiting.')
        sys.exit() 

    ### NADIR TRACK
    # Get single pixel for selecting correct nadir segment
    pixel_pt = gdf_PIXC.iloc[0].geometry      
    # Find correct nadir segment and return its geometry
    nadir_segment_ln = findNadir(pass_num=pass_num, pixel_pt=pixel_pt)
    
    del pixel_pt
    
    ## FIND AZIMUTH OF NADIR
    az_nadir = calcAzimuth(line=nadir_segment_ln)

    ### EVALUATE COVERAGE
    ## Make pseudo-pixels
    # Set along-track pixel resolution
    azimuth_res = 22 # meters
    
    pandarallel.initialize(nb_workers=int(os.environ.get('SLURM_CPUS_PER_TASK')))
    
    gdf_PIXC_cov = gdf_PIXC[['points', 'pixel_area', 'geometry']].copy()
    
    # Make pseudo pixels
    # Make pseudo pixels
    gdf_PIXC_cov['pseudo_pixel'] = gdf_PIXC_cov.parallel_apply(user_defined_function=makePseudoPixels,
                                                               args=(nadir_segment_ln, azimuth_res),
                                                               axis=1)
    
    del nadir_segment_ln
    
    # Clean-up
    gdf_PIXC_cov = gdf_PIXC_cov.rename(columns={'geometry': 'pixel_centroid'}).set_geometry('pseudo_pixel').set_crs(epsg=3857)
    # Get bounds of PIXC tile
    pseudo_bounds = gdf_PIXC_cov.total_bounds
    # Copy geometry column as sjoin will discard it
    gdf_PIXC_cov['pseudo_geom'] = gdf_PIXC_cov.geometry # DO I NEED THIS

    ## Read in buffered segments
    data_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_segmented_buffered/HUC2_' + \
    huc2 + '/'
    
    file_paths = []

    for huc in huc4s:
        file_path = data_path + 'NHDPLUS_H_' + huc + \
        '_HU4_GDB_prepped_segmented_buffered_' + width_set + '.parquet'
        file_paths.append(file_path)

    fields = ['NHDPlusID', width, binn, 'counter', 'buffers']
    segments = dask_geopandas.read_parquet(path=file_paths, columns=fields)
    segments = segments.compute()
    segments = segments.clip(pseudo_bounds)
    
    # Keep only reaches that are fully contained in PIXC granule
    segments = segments.groupby('NHDPlusID').filter(lambda x: len(x) == 10)
    segments = segments.sort_values(by=['NHDPlusID',
                                        'counter']).reset_index()
    
    segments = segments.drop(columns='index')
    
    # Calculate segment area
    segments['segment_area'] = segments.geometry.area
    
    ## Join and calculate coverage
    # Merge the segments and pseudo-puxels by intersection
    sj = gpd.sjoin(segments, gdf_PIXC_cov, how='left',
                   predicate='intersects').reset_index()
    
    del segments, gdf_PIXC_cov
    
    sj = sj.drop(columns=['index', 'points', 'pixel_centroid',
                          'pixel_area', 'index_right'])
    sj = sj.set_geometry('pseudo_geom')
    sj = sj.groupby('NHDPlusID',
                    as_index=False).parallel_apply(user_defined_function=specialDissolve)
    
    sj = sj.reset_index().drop(columns=['level_0', 'level_1'])
    sj['pseudo_geom_clip'] = sj.parallel_apply(user_defined_function=specialClip,
                                               axis=1)
    # Calculate the pseudo-pixel area within each node
    sj['pseudo_area'] = sj.pseudo_geom_clip.area
    sj['coverage'] = sj.pseudo_area/sj.segment_area
    sj['coverage'] = sj['coverage'].fillna(0)
    
    # Drop geometry columns
    sj = sj.drop(columns=['pseudo_geom', 'buffers', 'pseudo_geom_clip',
                          'pseudo_area', 'segment_area'])
    sj['nodes'] = 'node_' + sj['counter'].astype(str)
    sj = sj.drop(columns='counter')
    
    df = sj.pivot(index='NHDPlusID', columns='nodes', values='coverage').reset_index()
    df.index.name = 'index'
    
    df = pd.merge(left=df, right=sj[['NHDPlusID', width, binn]].drop_duplicates(subset='NHDPlusID'),
                  how='left', on='NHDPlusID')
    del sj
    
    ### EVALUATE HEIGHTS
    ## Read in flowlines
    data_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped/HUC2_' +\
                 huc2 + '/'
    
    file_paths = []
    for huc in huc4s:
        file_path = data_path + 'NHDPLUS_H_' + huc + '_HU4_GDB_prepped.parquet'
        file_paths.append(file_path)
        
    fields = ['NHDPlusID', 'GNIS_Name', 'LengthKM', 'geometry']
    flowlines = dask_geopandas.read_parquet(path=file_paths,
                                            columns=fields)
    flowlines = flowlines.compute()
    flowlines.loc[:,'geometry'] = flowlines.geometry.explode().force_2d()
    
    ## Calculate azimuth and sinuosity (should move sinuosity to static)
    flowlines['az_nadir'] = az_nadir
    flowlines['temp'] = flowlines.parallel_apply(user_defined_function=calcAzSin,
                                                 axis=1)
    flowlines[['alignment', 'sinuosity']] = pd.DataFrame(flowlines['temp'].tolist(),
                                                         index=flowlines.index)
    flowlines = flowlines.drop(columns='temp')
    # Rename geometry to avoid conflicts with heights
    flowlines = flowlines.rename_geometry('flowlines')
    
    gdf_PIXC = pd.merge(left=gdf_PIXC,
                    right=flowlines[['NHDPlusID', 'alignment',
                                     'sinuosity', 'flowlines']],
                    on='NHDPlusID')
    del flowlines
    
    ## Calculate distance along reaches
    gdf_PIXC['distance'] = gdf_PIXC.apply(lambda x: \
                                          project_point(x['flowlines'],
                                                        x['geometry']), axis=1)
    # Drop geometry columns
    df_PIXC = gdf_PIXC.drop(columns=['geometry', 'flowlines'])
    del gdf_PIXC
    
    ## Would calculate weights here
    
    # Calculate slopes
    df_PIXC_for_agg = df_PIXC[df_PIXC['agg_wse'] == 1]
    slopes = df_PIXC_for_agg.groupby('NHDPlusID').apply(func=group_regress,
                                                        include_groups=False)
    del df_PIXC_for_agg
    
    df_PIXC = pd.merge(left=df_PIXC, right=slopes, how='left', on='NHDPlusID')
    del slopes
    
    df_PIXC = df_PIXC.drop(columns=['pixel_area', 'klass', 'distance'])
    
    df_PIXC['DropM'] = df_PIXC['Slope'] * df_PIXC['LengthKM'] * 1000
    
    num_pixels = df_PIXC.groupby('NHDPlusID').count()['points'].reset_index()
    
    temp = df_PIXC.groupby('NHDPlusID')['agg_wse'].value_counts().reset_index()
    temp = temp[temp['agg_wse'] == 1]
    temp = temp.reset_index().drop(columns=['agg_wse', 'index'])
    
    df_PIXC = df_PIXC.drop(columns=['points', 'agg_wse', 'wse'])
    
    df_PIXC = pd.merge(df_PIXC, num_pixels, how='left', on='NHDPlusID')
    df_PIXC = pd.merge(df_PIXC, temp, how='left', on='NHDPlusID')
    
    df_PIXC = df_PIXC.rename(columns={'points': 'num_pixels', 'count': 'num_wse'})
    
    df_PIXC = df_PIXC.drop_duplicates(subset='NHDPlusID').reset_index()
    
    df = pd.merge(left=df, right=df_PIXC, how='left', on='NHDPlusID')
    
    del df_PIXC
    
    ### WRITE OUT
    save_path = os.path.join('/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_output/',
                             save_dir)

    if not os.path.isdir(save_path):
        os.makedirs(save_path)

    # sj.to_csv(os.path.join(save_path, granule_name + '_coverage.csv'))

    ### MAKE PARQUET
    df.to_parquet(os.path.join(save_path, granule_name + '_reaches.parquet'))
    # reaches_cent.to_parquet(os.path.join(save_path, granule_name + '_reaches_cent.parquet'))
    # reaches_thresh.to_parquet(os.path.join(save_path, granule_name + '_reaches_thresh.parquet'))
    # reaches_min.to_parquet(os.path.join(save_path, granule_name + '_reaches_min.parquet'))
    print('Script completed, wrote out results.')
    
    # end
    
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
    huc2 = '17'
    suffix = '_1'#'_2'
    data_path = '/nas/cee-water/cjgleason/data/SWOT/PIXC_v2_0_HUC2_' + huc2
    # pixc_ref = 'PIXC_v2_0_HUC2_01_best_files_no_exits.csv' ## CHANGE THIS
    save_dir = 'PIXC_v2_0_HUC2_' + huc2 + '_2025_05_24_'+ width_set
    
    evalCoverage(width_set=width_set, suffix=suffix, index=slurm, cpus_per_task=cpus_per_task, huc2=huc2, data_path=data_path, save_dir=save_dir)