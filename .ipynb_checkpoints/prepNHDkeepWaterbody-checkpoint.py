import geopandas as gpd
import os
import pandas as pd

'''
This function takes a filepath to the NHD HR Plus (currently located at
/nas/cee-water/cjgleason/craig/CONUS_ephemeral_data, merges on the VAA
and EROMMA tables, intersects the flowlines with the waterbodies and 
discards those that are within, finds the physiographic division for
each reach, and calculates the width for each reach. It writes out
the merged HUC4 files as gpkgs.
'''

slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])

def prepNHDwaterbody(data_path, save_path):
    
    ## Set-up
    mdata_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/'
    # Max binsize of 1000 is plenty for CONUS w/o lakes
    # Mississippi is ~472 m wide at mouth as calculated here
    bins = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 500, 1000]
    
    # Define dtypes for lookup tables to preserve leading zeros
    dtype_dic= {'HUC4': str, 'HUC2': str, 'toBasin': str, 'level': str}
    # Read in HUC lookup table
    lookup = pd.read_csv(os.path.join(mdata_path, 'HUC4_lookup_no_great_lakes.csv'), dtype=dtype_dic)
    
    # Get slurm job index
    i = slurm
    
    # Get current HUC2 and HUC4 IDs
    huc2 = 'HUC2_' + lookup.loc[i,'HUC4'][0:2]
    huc4 = 'NHDPLUS_H_' + lookup.loc[i,'HUC4'] + '_HU4_GDB'

    # Set data filepath
    file_path = os.path.join(data_path, huc2, huc4, huc4 + '.gdb')
    
    # Set write filepath
    save_path = os.path.join(save_path, huc2)
    save_file = huc4 + '_prepped_with_waterbody.gpkg'
    
    ## Prep Physiographic Regions
    # https://www.sciencebase.gov/catalog/item/631405bbd34e36012efa304e
    physio = gpd.read_file(filename=os.path.join(data_path,
                                                 'other_shapefiles/physio.shp'),
                           engine='pyogrio')
    # Set CRS to Web Mercator
    physio = physio.to_crs(epsg=3857)
    # Dissolve provinces by division
    physio = physio.dissolve(by='DIVISION').reset_index()
    # Drop all columns besides division and geometry
    physio = physio[['DIVISION', 'geometry']]
    
    ## Get bankfull width coefficients from Bieber et al. 2015, Table 3
    bankfull = pd.read_csv(os.path.join(mdata_path, 'bieger_2015_bankfull_width.csv'))

    ## Merging
    # Read in NHD flowlines
    basin = gpd.read_file(filename=file_path, layer='NHDFlowline', engine='pyogrio')
    # Set CRS to Pseudo-Mercator https://epsg.io/3857
    basin = basin.to_crs(epsg=3857)

    # Read in VAA
    vaa = gpd.read_file(filename=file_path, layer='NHDPlusFlowlineVAA', engine='pyogrio')
    # Merge on VAA
    basin = basin.merge(right=vaa, how='inner', on=['NHDPlusID', 'VPUID', 'ReachCode'])
    # Read in EROMMA
    eromma = gpd.read_file(filename=file_path, layer='NHDPlusEROMMA', engine='pyogrio')
    # Merge on EROMMA
    basin = basin.merge(right=eromma, how='inner', on=['NHDPlusID', 'VPUID'])

    ## Filtering
    # Drop reaches that aren't stream types or artificial path
    basin = basin.loc[(basin.FType == 460) | (basin.FType == 558)]
    # Drop reaches that are terminal paths
    basin = basin.loc[basin.TerminalFl == 0]
    # Drop reaches with discharge of zero
    basin = basin.loc[basin.QBMA > 0]
    # Drop reaches with stream order of zero
    basin = basin.loc[basin.StreamOrde > 0]

    ## Find the physiographic division each reach is within
    # Using intersects to foil the broken topology even after the dissolve
    # and neither shapely nor sf fully repaired it
    basin = basin.sjoin(df=physio, how='left',
                        predicate='intersects').drop(columns='index_right')
    # Drop all reaches where DIVISION == NaN (in Canada and off the coast)
    basin = basin[~basin.DIVISION.isnull()]

    ## Get bankfull widths
    # Merge on bankfull width coefficient
    basin = basin.merge(bankfull, on='DIVISION', how='left')
    # Calculate width from cumulative drainage area
    basin['WidthM'] = basin.a*basin.TotDASqKm**basin.b
        
    # Drop reaches that are shorter than their width
    basin = basin[basin['LengthKM']*1000 > basin['WidthM']]
    
    ## Bin reaches by width, set to string for parquet
    basin['Bin'] = pd.cut(basin['WidthM'], bins).astype(str)
    basin['Bin_Min'] = pd.cut(basin['WidthM_Min'], bins).astype(str)
    basin['Bin_Max'] = pd.cut(basin['WidthM_Max'], bins).astype(str)

    ## Write out gdf as gpkg file
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    basin.to_file(os.path.join(save_path, save_file), driver='GPKG')
            
data_path = '/nas/cee-water/cjgleason/craig/CONUS_ephemeral_data/'
save_path = '../narrow_rivers_PIXC_data/NHD_prepped_with_waterbody/'

prepNHDwaterbody(data_path=data_path, save_path=save_path)