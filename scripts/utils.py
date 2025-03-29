import json

import geopandas as gpd

from shapely.geometry import *

def specialBuffer(df, width_col, cap_style, segmented=False, extra=False):
    '''
    XXX
    '''
    # print(df.index)
    if segmented == True:
        reach = df.segments
    else:
        reach = df.geometry
    
    half_width = df[width_col]/2
    
    if extra == True:
        buffer = reach.buffer(distance=(half_width+50), cap_style=cap_style)
    else:
        buffer = reach.buffer(distance=half_width, cap_style=cap_style)
    return buffer
    
def specialClip(df):
    '''
    XXX
    '''
    left = gpd.GeoSeries(df.pseudo_geom)
    right = gpd.GeoSeries(df.buffer)
    pseudo_geom_clip = left.clip(right)
    return pseudo_geom_clip

def specialDissolve(reach):
    '''
    XXX
    '''
    reach = reach.dissolve(by='counter', as_index=False)
    return reach

def open_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
        return data