from shapely.geometry import *
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import xarray as xr

def readNHD(index):
    '''
    This function takes the index for an NHD HUC4 basin (see
    ./data/HUC4_lookup_no_great_lakes_PIXC.csv), reads it in,
    forces the geometry to 2D, and returns the basin dataframe.
    '''
    
    ## Set-up
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    prep_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped/' # _with_waterbody

    # Define dtypes for lookup tables to preserve leading zeros
    dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
    # Read in HUC lookup table
    lookup = pd.read_csv(os.path.join(mdata_path, 'HUC4_lookup_no_great_lakes.csv'), dtype=dtype_dic)

    # Get current HUC2 and HUC4 IDs
    huc2 = 'HUC2_' + lookup.loc[index,'HUC4'][0:2]
    huc4 = 'NHDPLUS_H_' + lookup.loc[index,'HUC4'] + '_HU4_GDB'
    
    # Set data filepath
    file_path = os.path.join(prep_path, huc2, huc4 + '_prepped.gpkg') # _with_waterbody

    ## Read in prepped NHD flowlines
    features = ['NHDPlusID', 'GNIS_Name', 'LengthKM', 'WidthM', 'Bin', 'geometry']
    basin = gpd.read_file(filename=file_path, columns=features, engine='pyogrio')
    
    # Make geometry 2D LineStrings
    basin['geometry'] = basin.geometry.explode().force_2d()
    
    return basin, huc4

def cut(line, distance):
    '''
    This function takes one reach centerline and 1/10th of the reach length
    and cuts the lie in two at the distance from its starting point. It then
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
        
def bitwiseMask(ds):
    '''
    This function masks a PIXC granules: for now, it ony remove pixels
    with land classification and those with bad geolocation_qual.
    # See page 65 of PIXC PDD: https://podaac.jpl.nasa.gov/SWOT?tab=datasets-information&sections=about%2Bdata
    '''
    mask = np.where(np.logical_and(ds.classification > 1, ds.geolocation_qual < 2**12))[0]
    # print(mask.shape)
    return mask

def makeGDF(ds, mask, data_var):
    '''
    This function takes the pixel cloud xarray object, makes a masked
    GeoDataFrame, renames columns as needed, set the CRS, reprojects
    the CRS, and returns a GeoDataFrame.
    '''

    # Create GDF
    gdf_PIXC = gpd.GeoDataFrame(ds[data_var][mask],
                                geometry=gpd.points_from_xy(
                                    ds.longitude[mask],
                                    ds.latitude[mask]),
                                crs="EPSG:4326") # PIXC has no native CRS, setting same as River_SP

    if data_var == 'classification':
        gdf_PIXC.rename(columns={gdf_PIXC.columns[0]: 'klass'}, inplace=True)
    else:
        gdf_PIXC.rename(columns={gdf_PIXC.columns[0]: data_var}, inplace=True)
    
    # Convert the crs to WGS 84 / Pseudo-Mercator
    gdf_PIXC = gdf_PIXC.to_crs(epsg=3857)
    
    return gdf_PIXC

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

    # Segment the reach
    for i in range(9):
        try:
            # Chop the reach, store remainder
            new, line = cut(line=line, distance=dist)
            # Append new segment to GeoSeries of all segments
            segments.append(new)
            # segments[i]: new
        except:
            print(reach['NHDPlusID'])
    
    # Append final segment to list
    segments.append(line)

    return segments

def getCoverage(reach, basin_crs, gdf_PIXC):
    '''
    Ths function takes a segmented reach, explodes it to get one row
    per reach, buffers the reach segments with 1/2 the calculated width
    on either of the reach segment centerlines, intersects the buffered
    segments with the SWOT PIXC, and calculates the percent coverage for
    the whole reach, and returns this proportion.
    
    Right now, a segment only needs to have one pixel in it to be
    considered "detected".
    '''
    
    # Make GeoSeries with just segments
    segments = gpd.GeoSeries(data=reach['segments'], crs=basin_crs)

    # Buffer segments by 1/2 the calculated width
    # CITE BEIGER
    segments = segments.buffer(distance=(reach.WidthM/2),cap_style='flat')
    
    # Make GeoDataFrame of segments
    segments = gpd.GeoDataFrame(geometry=segments)
    
    # Get pixels in reach
    inside = gpd.sjoin(gdf_PIXC, segments, how='inner',
                       predicate='within').rename(columns={'index_right': 'segment'})

    if inside.empty:
        prop = 0
    else:    
        # Get number of pixels in each reach segment
        counts = pd.DataFrame(inside.segment.value_counts().reset_index().sort_index())
        # Calculate coverage (proportion of reaches with >= 1 pixel in them)
        prop = len(counts.loc[counts['count'] != 0])/10
    return prop