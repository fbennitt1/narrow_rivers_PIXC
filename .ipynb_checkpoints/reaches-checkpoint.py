from shapely.geometry import *
import geopandas as gpd
import numpy as np
import os
import pandas as pd
import shapely as shp
import xarray as xr

def readNHD(index):
    '''
    This function takes the index for an NHD HUC4 basin (see
    ./data/HUC4_lookup_no_great_lakes_PIXC.csv), reads it in,
    forces the geometry to 2D, and returns the basin dataframe.
    '''
    
    ## Set-up
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    prep_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped/'

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
    
    # Set data filepath
    file_path = os.path.join(prep_path, huc2, huc4 + '_prepped.gpkg')

    ## Read in prepped NHD flowlines
    features = ['NHDPlusID', 'GNIS_Name', 'LengthKM', 'WidthM', 'Bin', 'geometry']
    basin = gpd.read_file(filename=file_path, columns=features, engine='pyogrio')
    print('read in')
    
    # Make geometry 2D LineStrings
    basin['geometry'] = basin.geometry.explode().force_2d()
    print('exploded')
    
    return basin, huc4, huc2

def readSegments(index):
    '''
    This function takes the index for an NHD HUC4 basin (see
    ./data/HUC4_lookup_no_great_lakes_PIXC.csv), reads in the
    segmented and exploded version, and returns the basin
    dataframe.
    '''
    
    ## Set-up
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
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
    
    # Set data filepath
    file_path = os.path.join(prep_path, huc2, huc4 + '_prepped_segmented.parquet')

    ## Read in prepped NHD flowlines
    basin = gpd.read_parquet(path=file_path)
    print('read in')
    
    return basin, huc4, huc2

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
        
def bitwiseMask(ds):
    '''
    This function masks a PIXC granules: for now, it ony remove pixels
    with land classification and those with bad geolocation_qual.
    # See page 65 of PIXC PDD: https://podaac.jpl.nasa.gov/SWOT?tab=datasets-information&sections=about%2Bdata
    '''
    # Fow now, eliminate the really bad stuff
    mask = np.where((ds.classification > 2) & (ds.geolocation_qual < 2**16) &
                    (np.abs(ds.cross_track) > 10000) & (np.abs(ds.cross_track) < 60000))[0]
    
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
    gdf_PIXC = gdf_PIXC.to_crs(epsg='32618')
    
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

def specialDissolve(reach):
    reach = reach.dissolve(by='counter', as_index=False)
    return reach

def specialClip(sj):
    left = gpd.GeoSeries(sj.pseudo_geom)
    right = gpd.GeoSeries(sj.buffered)
    pseudo_geom_clip = left.clip(right)
    return pseudo_geom_clip

# OLD VERSION
# def getCoverage(reach, basin_crs, gdf_PIXC, num):
#     '''
#     Ths function takes a segmented reach, explodes it to get one row
#     per reach, buffers the reach segments with 1/2 the calculated width
#     on either of the reach segment centerlines, intersects the buffered
#     segments with the SWOT PIXC, and calculates the percent coverage for
#     the whole reach, and returns this proportion.
    
#     Right now, a segment only needs to have one pixel in it to be
#     considered "detected".
#     '''
    
#     # Make GeoSeries with just segments
#     segments = gpd.GeoSeries(data=reach['segments'], crs=basin_crs)

#     # Buffer segments by 1/2 the calculated width
#     # CITE BEIGER
#     segments = segments.buffer(distance=(reach.WidthM/2),cap_style='flat')
    
#     # Make GeoDataFrame of segments
#     segments = gpd.GeoDataFrame(geometry=segments)
    
#     # Get pixels in reach
#     inside = gpd.sjoin(gdf_PIXC, segments, how='inner',
#                        predicate='within').rename(columns={'index_right': 'segment'})

#     if inside.empty:
#         prop = 0
#     else:    
#         # Get number of pixels in each reach segment
#         counts = pd.DataFrame(inside.segment.value_counts().reset_index().sort_index())
        
#         # Calculate coverage (proportion of reaches with >= 1 pixel in them)
#         prop = len(counts.loc[counts['count'] >= 1])/10
        
#     return prop