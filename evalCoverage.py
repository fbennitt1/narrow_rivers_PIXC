import os
import sys
import time

# import contextily as ctx
import geopandas as gpd
# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapely
import xarray as xr

# from matplotlib import colors
from pandarallel import pandarallel
from shapely.geometry import box

from reaches import readNHD
from reaches import readSegments
from reaches import bitwiseMask
from reaches import makeGDF
from reaches import makePseudoPixels
from reaches import specialDissolve
from reaches import specialClip

### PIXEL CLOUD
# Get PIXC index metadata
mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
dtype_dic= {'cycle': str, 'pass': str, 'tile': str, 'version': str}
# Read in HUC lookup table
pixc_lookup = pd.read_csv(os.path.join(mdata_path,
                                       'PIXC_v2_0_HUC2_01_best_files.csv'),
                          dtype=dtype_dic).drop(columns='index')

# Get job index
slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])

# Get filepath
file_name = pixc_lookup.loc[slurm, 'files']

# PIXC datapath
data_path = '/nas/cee-water/cjgleason/fiona/data/PIXC_v2_0_HUC2_01/'

pixc_path = os.path.join(data_path, file_name)

tile_name = pixc_path[-71:-3]

## Check if tile intersects NHD
# Read in xarray, global
ds_GLOB = xr.open_mfdataset(paths=pixc_path, engine='h5netcdf')

# Get bounding coordinates for SWOT tile
west_lon = ds_GLOB.geospatial_lon_min
south_lat = ds_GLOB.geospatial_lat_min
east_lon = ds_GLOB.geospatial_lon_max
north_lat = ds_GLOB.geospatial_lat_max

# Clip width polygons to current points
bbox = box(west_lon, south_lat, east_lon, north_lat)
bbox = gpd.GeoDataFrame({'geometry': [bbox]}, crs="EPSG:4326")
bbox = bbox.to_crs(epsg='32618')
                   
# Read in HUC4 boundaries (no Great Lakes)
data_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/all_wbd_no_great_lakes.parquet'
wbd = gpd.read_parquet(path=data_path)
                   
# Project CRS
wbd = wbd.to_crs(epsg=32618)

# Check if granule intersects NHD 
test = gpd.sjoin(wbd, bbox, how='inner', predicate='intersects')
if test.shape[0] == 0:
    print('This granule does not intersect the NHD, exiting.')
    sys.exit()
                   
## If so, proceed
pass_num = ds_GLOB.pass_number
pass_num

# Read in xarray, pixel group
ds_PIXC = xr.open_mfdataset(paths=pixc_path,
                            group='pixel_cloud', engine='h5netcdf')

# Make mask
mask = bitwiseMask(ds_PIXC)

if mask.shape[0] == 0:
    print('This granule has no pixels after masking, exiting.')
    sys.exit()

variables = ['azimuth_index', 'range_index', 'cross_track',
             'pixel_area', 'height', 'geoid',
             'dlatitude_dphase', 'dlongitude_dphase',
             'dheight_dphase', 'classification']
                   
# Make PIXC into GeoDataFrame
gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_vars=variables)
                   
### FIND CORRECT HUC4
# Get NHD index metadata
# Define dtypes for lookup tables to preserve leading zeros
dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
# Read in HUC lookup table
lookup = pd.read_csv(os.path.join(mdata_path,
                                  'HUC4_lookup_no_great_lakes.csv'),
                     dtype=dtype_dic)

# Merge on index metadata
wbd = pd.merge(left=wbd, right=lookup, how='inner', left_on='huc4',
               right_on='HUC4').drop(columns=['HUC4', 'HUC2',
                                              'toBasin', 'level'])
        
# Get bounds of PIXC tile
bounds_PIXC = gdf_PIXC.union_all().convex_hull
gdf_bounds = gpd.GeoDataFrame({'geometry': [bounds_PIXC]}, crs=wbd.crs)
                   
# Get slurm indices (from mdata) of basins that intersect PIXC tile
indices = gpd.sjoin(wbd, gdf_bounds, how='inner', predicate='intersects')['slurm_index'].to_list()
                   
### READ IN HUC4 BASINS
# Create merged dataframe of all basins intersected
if len(indices) == 1:
    # Read prepped NHD
    basin, huc4_list, huc2_list = readNHD(index=indices[0])

else:
    # Initialize lists
    d = []
    huc4_list = []
    huc2_list = []
    
    # Loop through indices and store in lists
    for idx in indices:

        # Read prepped NHD
        basin, huc4, huc2 = readNHD(index=idx)

        # Append to lists
        d.append(basin)
        huc4_list.append(huc4)
        huc2_list.append(huc2)
        
    # Merge GeoDataFrames
    basin = pd.concat(d)

# Project crs
basin = basin.to_crs(epsg=32618)
                   
# Buffer with an extra 50 m on each side to be safe
basin['buffer'] = basin.buffer(distance=((basin.WidthM/2)+50), cap_style='flat')
                   
# Set geometry to buffered reaches
basin = basin.set_geometry('buffer')
                   
# Get only pixels within buffered reaches
gdf_PIXC_clip = gpd.sjoin(gdf_PIXC, basin, how='inner', predicate='within')

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
                   
# Read in nadir (science orbit)
nadir = gpd.read_file('/nas/cee-water/cjgleason/data/SWOT/swath/swot_science_hr_Aug2021-v05_shapefile_nadir/swot_science_hr_2.0s_4.0s_Aug2021-v5_nadir.shp')

# Convert CRS to WGS 84 / UTM zone 18N
nadir = nadir.to_crs(epsg=32618)
                   
# Find candidate nadir segments
candidates = nadir[nadir['ID_PASS'] == pass_num]
                   
# Find distance from each candidate to single pixel
candidates['dist'] = candidates.loc[:,'geometry'].distance(pixel_pt)
                   
# Get nadir segment closest to single pixel
nadir_segment = candidates[candidates.dist == candidates.dist.min()]
                   
# Get nadir segment geoemtry
nadir_segment_ln = nadir_segment.geometry[nadir_segment.index[0]]
                   
### MAKE PSEUDO-PIXRLS
# Set along-track pixel resolution
azimuth_res = 22 # meters
                   
pandarallel.initialize()
                   
# Make pseudo pixels
gdf_PIXC_clip['pseudo_pixel'] = gdf_PIXC_clip.parallel_apply(user_defined_function=makePseudoPixels,
                                                             args=(nadir_segment_ln, azimuth_res),
                                                             axis=1)
            
pseudo = gdf_PIXC_clip.drop(columns='geometry').set_geometry('pseudo_pixel').set_crs(crs=gdf_PIXC_clip.crs)
                   
# Get bounds of PIXC tile
# pseudo_bounds = pseudo.union_all().convex_hull
pseudo_bounds = pseudo.total_bounds
pseudo_poly = box(pseudo_bounds[0], pseudo_bounds[1],
                      pseudo_bounds[2], pseudo_bounds[3])
                   
### READ IN SEGMENTS
# Create merged dataframe of all basins intersected
if len(indices) == 1:
    # Read prepped NHD
    segments, _, _ = readSegments(index=indices[0])

else:
    # Initialize lists
    d = []
    
    # Loop through indices and store in lists
    for idx in indices:

        # Read prepped NHD
        segments, huc4, _ = readSegments(index=idx)
        # Make column with HUC4 id
        segments['huc4_long'] = huc4
        segments['huc4'] = segments['huc4_long'].str[10:14]
        # Append to list
        d.append(segments)
        
    # Merge GeoDataFrames
    segments = pd.concat(d)
                   
segments = segments.to_crs(epsg='32618')
segments = segments.reset_index().rename(columns={'index': 'index_old'})
                   
# Assign a unique counter within each index group
segments['counter'] = segments.groupby('NHDPlusID').cumcount()
                   
# Keep only first ten segments (some reaches repeat)
segments = segments[segments['counter'] < 10]
                   
segments_clip = segments.clip(pseudo_bounds)
                   
# Keep only reaches that are fully contained in PIXC granule
segments_clip = segments_clip.groupby('NHDPlusID').filter(lambda x: len(x) == 10)
                   
# Buffer segments
segments_clip['buffered'] = segments_clip.buffer(distance=(segments_clip.WidthM/2), cap_style='flat')
                   
segments_clip = segments_clip.set_geometry('buffered')
                   
# Calculate segment area
segments_clip['segment_area'] = segments_clip.geometry.area
                   
# Copy geometry column as sjoin will discard it
pseudo['pseudo_geom'] = pseudo.geometry
                   
# Merge the segments and pseudo-puxels by intersection
sj = gpd.sjoin(segments_clip, pseudo, how='left', predicate='intersects')

# Drop unneeded columns
sj = sj.drop(columns=['index_right', 'points', 'azimuth_index',
                      'range_index', 'cross_track', 'pixel_area',
                      'height', 'geoid', 'dlatitude_dphase',
                      'dlongitude_dphase', 'dheight_dphase',
                      'klass', 'latitude', 'longitude', ])

# Set active geometry column for dissolve
sj = sj.set_geometry('pseudo_geom')

# Dissolve
dissolved = sj.groupby('NHDPlusID', as_index=False).parallel_apply(user_defined_function=specialDissolve)

# Drop multi-index
dissolved = dissolved.reset_index().drop(columns=['level_0', 'level_1'])

# Clip dissolved pseudo-pixels to node areas
dissolved['pseudo_geom_clip'] = dissolved.parallel_apply(user_defined_function=specialClip,
                                                         axis=1)
# Calculate the pseudo-pixel area within each node
dissolved['pseudo_area'] = dissolved.pseudo_geom_clip.area

# Calculate coverage
dissolved['coverage'] = dissolved.pseudo_area/dissolved.segment_area

# Copy to fill zeros for stats
dissolved_w_zero = dissolved.copy()
dissolved_w_zero['coverage'] = dissolved_w_zero['coverage'].fillna(0)

### DO STATS
bins = dissolved.Bin.unique()

## Nodes
node_desc = dissolved.groupby('Bin')['coverage'].describe().reset_index()
node_desc['with_zero'] = 0

node_quant = pd.DataFrame(dissolved.groupby('Bin')['coverage'].quantile(q=[x / 100.0 for x in range(0,100,1)])).reset_index().rename(columns={'level_1': 'quantile'})
node_quant['with_zero'] = 0

## Nodes with zeros
node_desc_w_zero = dissolved_w_zero.groupby('Bin')['coverage'].describe().reset_index()
node_desc_w_zero['with_zero'] = 1

node_quant_w_zero = pd.DataFrame(dissolved_w_zero.groupby('Bin')['coverage'].quantile(q=[x / 100.0 for x in range(0,100,1)])).reset_index().rename(columns={'level_1': 'quantile'})
node_quant_w_zero['with_zero'] = 1

# nodes_mean = dissolved.groupby('Bin')['coverage'].mean().to_list()
# nodes_std = dissolved.groupby('Bin')['coverage'].std().to_list()
# d = {'mean': nodes_mean, 'std': nodes_std}
# nodes = pd.DataFrame(data=d).T
# nodes.columns = bins

## Reaches
d = {}
for i in range(1, 10):
    threshold = i/10
    # print(threshold)
    
    detected = dissolved.groupby(['Bin', 'NHDPlusID'])['coverage'].apply(lambda x: (x > threshold).sum()) / 10
    reach = detected.groupby('Bin').mean().to_list()
    #quantiles
    #std dev
    #n
    # size of the dataframe (are we including the reaches with nothing detected?????)
    
    d[threshold] = reach
    
reaches = pd.DataFrame(data=d).T
reaches.columns = bins
                   
### WRITE OUT
save_path = os.path.join('/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/', 'PIXC_v2_0_HUC2_01')

if not os.path.isdir(save_path):
    os.makedirs(save_path)
    
# sj.to_csv(os.path.join(save_path, tile_name + '_coverage.csv'))
nodes.to_csv(os.path.join(save_path, tile_name + '_nodes.csv'))
reaches.to_csv(os.path.join(save_path, tile_name + '_reaches.csv'))