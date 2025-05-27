import os

import geopandas as gpd
import math
import numpy as np
import pandas as pd
import shapely as shp
import xarray as xr

from scipy.stats import linregress
from shapely.geometry import *

def readNHD(index, segmented=False):
    '''
    This function takes the index for an NHD HUC4 basin (see
    ./data/HUC4_lookup_no_great_lakes_PIXC.csv), reads it in,
    forces the geometry to 2D, and returns the basin dataframe.
    '''
    
    ## Set-up
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    if segmented == False:
        print('type: normal')
        prep_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped/'
    else:
        print('type: segmented')
        prep_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_segmented/'

    # Define dtypes for lookup tables to preserve leading zeros
    dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
    # Read in HUC lookup table
    lookup = pd.read_csv(os.path.join(mdata_path,
                                      'HUC4_lookup_no_great_lakes.csv'),
                         dtype=dtype_dic)

    # Get current HUC2 and HUC4 IDs
    huc2 = 'HUC2_' + lookup.loc[index,'HUC4'][0:2]
    huc4 = 'NHDPLUS_H_' + lookup.loc[index,'HUC4'] + '_HU4_GDB'
    print(str(huc4))

    if segmented == False:
        # Set data filepath
        file_path = os.path.join(prep_path, huc2, huc4 + '_prepped.parquet')

        ## Read in prepped NHD flowlines
        features = ['NHDPlusID', 'GNIS_Name', 'LengthKM', 'WidthM', 'WidthM_Min',
                    'WidthM_Max', 'Bin', 'Bin_Min', 'Bin_Max', 'StreamOrde',
                    'Slope', 'geometry']
        basin = gpd.read_parquet(path=file_path, columns=features)
        print('flowlines read-in')

        # Make geometry 2D LineStrings
        basin['geometry'] = basin.geometry.explode().force_2d()
        print('exploded')

    else:
        # Set data filepath
        file_path = os.path.join(prep_path, huc2, huc4 + '_prepped_segmented.parquet')
        print(file_path)

        ## Read in prepped NHD flowlines
        basin = gpd.read_parquet(path=file_path)
        print('segments read-in')
        
    return basin, huc4, huc2

def findNadir(pass_num, pixel_pt):
    '''
    XXX
    '''
    # Read in nadir (science orbit)
    nadir = gpd.read_file('/nas/cee-water/cjgleason/data/SWOT/swath/swot_science_hr_Aug2021-v05_shapefile_nadir/swot_science_hr_2.0s_4.0s_Aug2021-v5_nadir.shp')
    # Project CRS (currently to WGS 84 / UTM zone 18N)
    nadir = nadir.to_crs(epsg=3857)
    # Find candidate nadir segments
    candidates = nadir[nadir['ID_PASS'] == pass_num]
    # Find distance from each candidate to single pixel
    candidates.loc[:,'dist'] = candidates.loc[:,'geometry'].distance(pixel_pt)
    # Get nadir segment closest to single pixel
    nadir_segment = candidates[candidates.dist == candidates.dist.min()]
    # Get nadir segment geometry
    nadir_segment_ln = nadir_segment.geometry[nadir_segment.index[0]]
    
    return nadir_segment_ln

def cut(line, distance):
    '''
    This function takes one reach centerline and 1/10th of the reach length
    and cuts the line in two at the distance from its starting point. It then
    returns the trimmed piece and the remainder of the reach centerline.
    '''
    if distance <= 0.0 or distance >= line.length:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p))
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]
        
def bitwiseMask(ds): # HERE
    '''
    This function masks a PIXC granules: for now, it ony remove pixels
    with land classification and those with bad geolocation_qual.
    # See page 65 of PIXC PDD: https://podaac.jpl.nasa.gov/SWOT?tab=datasets-information&sections=about%2Bdata
    '''
    # Fow now, eliminate the really bad stuff
    # mask = np.where((ds.classification > 1) & (ds.geolocation_qual < 2**16) &
    #                 (np.abs(ds.cross_track) > 10000) & (np.abs(ds.cross_track) < 60000))[0]
    mask = np.where((ds.classification > 1) &
                    (ds.interferogram_qual < 2**16) &
                    (ds.classification_qual < 2**16) &
                    (ds.geolocation_qual < 2**16) &
                    (ds.sig0_qual < 2**16) &
                    (np.abs(ds.cross_track) > 10000) &
                    (np.abs(ds.cross_track) < 60000))[0]
    
    print(mask.shape)
    return mask

def makeGDF(ds, mask, data_vars):
    '''
    This function takes the pixel cloud xarray object, makes a masked
    GeoDataFrame, renames columns as needed, set the CRS, reprojects
    the CRS, and returns a GeoDataFrame.
    '''

    # Subset xarray, convert to masked DataFrame
    xarr = ds[data_vars]
    df = xarr.to_dataframe().loc[mask].reset_index()

    # Create GDF
    gdf_PIXC = gpd.GeoDataFrame(df,
                                geometry=gpd.points_from_xy(df.longitude,
                                                            df.latitude),
                                crs="EPSG:4326") # PIXC has no native CRS, setting same as River_SP

    if 'classification' in gdf_PIXC.columns:
        gdf_PIXC.rename(columns={'classification': 'klass'}, inplace=True)
    
    # Convert the crs to WGS 84 / UTM zone 18N
    gdf_PIXC = gdf_PIXC.to_crs(epsg='3857')
    
    return gdf_PIXC

def calcAzimuth(line):
    # Set-up
    x_coords = [coord[0] for coord in line.coords]
    y_coords = [coord[1] for coord in line.coords]
    
    # Find deltas
    dx = x_coords[-1] - x_coords[0]
    dy = y_coords[-1] - y_coords[0]
    
    # Find azimuth
    azimuth = math.degrees(math.atan2(dy, dx))
        
    if azimuth < 0:
        azimuth += 360
    
    return azimuth

def calcAzSin(df):
    '''
    '''
    # Set-up
    line = df.geometry
    x_coords = [coord[0] for coord in line.coords]
    y_coords = [coord[1] for coord in line.coords]
    
    # Regress
    result = linregress(x_coords, y_coords)
    slope = result.slope
    intercept = result.intercept
    
    # Find deltas
    dx = x_coords[-1] - x_coords[0]
    dy = y_coords[-1] - y_coords[0]
    dy_reg = (slope*x_coords[-1] + intercept) - y_coords[0]
    
    # Find azimuth
    azimuth = math.degrees(math.atan2(dy_reg, dx))
        
    if azimuth < 0:
        azimuth += 360
        
    alignment = azimuth - df.az_nadir
    
    # Find sinuosity
    distance = np.sqrt(dx**2 + dy**2)
    
    sinuosity = line.length/distance
    
    return alignment, sinuosity

def segmentReach(reach):
    '''
    Ths function takes a reach center line, cuts it into ten segments of
    equal length, and returns those ten segments.
    '''
    # Get linestring
    line = reach.geometry
    # Find length of segments
    dist = line.length/10

    # Make empty list for segments
    segments = []
    
    # # Set failed to default
    failed = 0

    # Segment the reach
    for i in range(9):
        try:
            # Chop the reach, store remainder
            new, line = cut(line=line, distance=dist)
            # Append new segment to GeoSeries of all segments
            segments.append(new)
            # segments[i]: new
        except:
            failed = 1
            # print('Failed to segment reach (#' + str(i) +'): ' + str(reach['NHDPlusID']))
    
    # Append final segment to list
    segments.append(line)

    return segments, failed

def makePseudoPixels(pixel, segment_ln, azimuth_res):
    '''
    XXX
    '''
    # Get pixel geometry
    pixel_pt = pixel.geometry #[pixel.index[0]]
    
    # Calculate width
    width = pixel.pixel_area/azimuth_res #[pixel.index[0]]/azimuth_res
    
    # Create linestring from pixel to closest point on nadir track
    orthogonal = shp.shortest_line(pixel_pt, segment_ln)
    
    # return orthogonal
    
    # Make line parallel to orthogonal at correct azimuth_res
    up = orthogonal.parallel_offset(distance=azimuth_res/2, side='right')
    down = orthogonal.parallel_offset(distance=azimuth_res/2, side='left')
    
    # Get coords for inner edge
    one_coord = shp.ops.substring(up, start_dist=0, end_dist=width/2).coords[1]
    two_coord = shp.ops.substring(down, start_dist=0, end_dist=width/2).coords[1]
    
    # Get inner and outer edges of polygon
    inner_edge = LineString([two_coord, one_coord])
    outer_edge = inner_edge.parallel_offset(distance=width, side='right')
    
    # Get coords for outer edge
    three_coord = outer_edge.coords[0]
    four_coord = outer_edge.coords[1]
    
    # Make pseudo pixel
    pseudo_pixel = Polygon((one_coord, two_coord, three_coord, four_coord, one_coord))
    
    return pseudo_pixel


def group_regress(group):
    
    # group = group.filter(lambda x: x['agg_wse'] == 1)
    
    if (len(group) > 2) and (len(group['distance'].unique()) > 1):
        slope, _, r, _, se = (linregress(x=group['distance'], y=group['wse']))
        return pd.Series({'swot_slope': slope, 'swot_r': r, 'swot_se': se})
    else:
        return pd.Series({'swot_slope': np.nan, 'swot_r': np.nan, 'swot_se': np.nan})
    

def summarizeCoverage(df, binn, bins, counts):
    '''
    '''
    ### REACHES
    ## Centiles, thresholds
    d_c = {}
    d = {}

    for i in range(1, 10):
        threshold = i/10

        detected = df.groupby([binn, 'NHDPlusID'])['coverage'].apply(lambda x: (x > threshold).sum()) / 10

        reach_c = detected.groupby(binn).quantile(q=[x / 100.0 for x in range(0,100,1)]).reset_index()
        d_c[threshold] = reach_c
        
        reach = detected.reset_index()
        d[threshold] = reach
    
    # Centiles
    for threshold, data in d_c.items():
        data['threshold'] = threshold
        
    reaches_cent  = pd.concat(d.values()).rename(columns={'level_1': 'quantile'})    
    # Merge on reach counts
    reaches_cent = pd.merge(left=reaches_cent, right=counts, how='left', on=binn)
    
    # All
    for threshold, data in d.items():
        data['threshold'] = threshold
        
    reaches_thresh = pd.concat(d.values()).rename(columns={'level_1': 'quantile'})
    
    # Minimum coverage
    reaches_min = pd.DataFrame(df.groupby('NHDPlusID')['coverage'].min()).reset_index()
    # Merge on bins
    reaches_min = pd.merge(left=reaches_min, right=df[['NHDPlusID', binn]], how='left', on='NHDPlusID')
    # Take every tenth row to get reach-level results
    reaches_min = reaches_min.sort_values(by=['NHDPlusID'])[::10].reset_index()

    # return reaches_cent, reaches_min
    return reaches_cent, reaches_thresh, reaches_min