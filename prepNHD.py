import geopandas as gpd
import os
import pandas as pd

# This function takes some things
# merges on certain columns from the VAA and EROMMA tables,
# finds the physiographic division of each reach,
# and calculates the width for each reach. It writes out
# the merged HUC4 files as gpkgs.

def prepNHD(data_path):
    
    ## Set-up
    codes_huc2 = ['01','02','03','04','05','06','07','08','09',
              '10','11','12','13','14','15','16','17','18']
    fieldsF = ['GNIS_ID', 'GNIS_Name', 'LengthKM',  'FlowDir',
               'WBArea_Permanent_Identifier', 'FType', 'FCode',
               'NHDPlusID', 'VPUID', 'geometry']
    fieldsVAA = ['NHDPlusID', 'StreamOrde', 'FromNode', 'ToNode',
                'LevelPathI', 'TerminalFl', 'TotDASqKm', 'VPUID']
    fieldsEROMMA = ['NHDPlusID', 'QBMA', 'VPUID']
    
    ## Prep Physiographic Regions
    # https://www.sciencebase.gov/catalog/item/631405bbd34e36012efa304e
    physio = gpd.read_file('/nas/ceewater/cjgleason/craig/CONUS_ephemeral_data/other_shapefiles/physio.shp')
    # Dissolve provinces by division
    physio = physio.dissolve(by='DIVISION').reset_index()
    # Set CRS to Web Mercator
    physio = physio.to_crs(epsg=3857)
    # Drop all columns besides division and geometry
    physio = physio[['DIVISION', 'geometry']]
    
    ## Get bankfull width coefficients from Bieber et al. 2015, Table 3
    bankfull = pd.read_csv('/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC/data/bieger_2015_bankfull_width.csv')
    

    ## Loop through HUC4 basins, prep the data, and write out new file
    for i in range(len(codes_huc2)):
        
    # Get all HUC4 paths for current HUC2 (excluding WBD)
    sub_paths = [fn for fn in os.listdir(os.path.join(datapath, 'HUC2_' +       codes_huc2[i])) if fn.startswith('NHD')]
    
        for j in sub_paths:
            path = os.path.join(datapath, 'HUC2_' + codes_huc2[i],
                            j, j + '.gdb')
            ## Merging
            # Read in NHD flowlines
            basin = gpd.read_file(filename=data_path, layer='NHDFlowline',
                                  columns=fieldsF)
            # Read in VAA
            vaa = gpd.read_file(filename=data_path, layer='NHDPlusFlowlineVAA',
                                columns=fieldsVAA)
            # Merge on VAA
            basin = basin.merge(vaa, on=['NHDPlusID', 'VPUID'])
            # Read in EROMMA
            eromma = gpd.read_file(filename=data_path, layer='NHDPlusEROMMA',
                        columns=fieldsEROMMA)
            # Merge on EROMMA
            basin = basin.merge(eromma, on=['NHDPlusID', 'VPUID'])
            # Set CRS to Pseudo-Mercator https://epsg.io/3857
            basin = basin.to_crs(epsg=3857)
            
            ## Filtering
            # Keep only reaches that are stream types or artificial path
            basin = basin.loc[(basin.FType == 460 | basin.FType == 558)]
            # Keep only reaches that are not terminal paths
            basin = basin.loc[basin.TerminalFl == 0]
            # Keep only reaches with non-zero discharge
            basin = basin.loc[basin.QBMA > 0]
            # Keep only reaches with non-zero stream order
            basin = basin.loc[basin.StreamOrde > 0]
            
            ## Find the physiographic division each reach is within
            basin = gpd.sjoin(left_df=basin, right_df=physio,
                             how='left', predicate='within') # do I want within
            
            ## Merge on bankfull width coefficient
            basin = basin.merge(bankfull, on='DIVISION', how='left')
            # Calculate width from cumulative drainage area
            basin['WidthM'] = basin.a*basin.TotDASqKm**basin.b
            
            ## Write out gdf as gpkg file
            save_path = '../narrow_rivers_PIXC_data/'
            save_file = j + '.gpkg'
            if not os.path.isdir(save_path):
                os.makedirs(save_path)
            basins.to_file(os.path.join(save_path, save_file), driver='GPKG')
            
            