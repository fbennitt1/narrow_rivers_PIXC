import os

import geopandas as gpd
import pandas as pd
from pandarallel import pandarallel

from reaches import readNHD
from reaches import cut
from reaches import segmentReach

'''
'''
def segmentReaches(index, cpus_per_task, save_path):
    ## Prepare data
    # Read prepped NHD
    gdf, huc4, huc2 = readNHD(index=index, segmented=False)

    # Initialize panarallel
    pandarallel.initialize(nb_workers=cpus_per_task)

    # Segment the reaches
    gdf[['segments', 'failed']] = gdf.parallel_apply(user_defined_function=segmentReach, axis=1).apply(pd.Series)

    # Keep only reaches that were sucessfully segmented
    gdf = gdf[gdf['failed'] == 0]

    # Drop failure column
    gdf = gdf.drop(columns='failed')

    # Explode the segments to give each own row
    gdf =  gdf.explode(column='segments', index_parts=True)

    # Drop original reach geometry column, set segments as active geometry
    gdf = gdf.drop(columns='geometry').set_geometry('segments').set_crs(crs=gdf.crs)

    # Write out
    # Set write filepath
    # save_path = '../narrow_rivers_PIXC_data/NHD_prepped_segmented/'
    save_path = os.path.join(save_path, huc2)
    save_file = huc4 + '_prepped_segmented.parquet'

    #Write out gdf as parquet file
    if not os.path.isdir(save_path):
        os.makedirs(save_path)
    gdf.to_parquet(os.path.join(save_path, save_file))
    
    print('Script completed.')

if __name__ == "__main__":
    ### GET JOB INDEX, CPUS
    slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])
    cpus = int(os.environ.get('SLURM_CPUS_PER_TASK'))
    cpus_per_task = cpus if cpus < 65 else 1
    save_path = '/nas/cee-water/cjgleason/fiona/narrow_rivers_PIXC_data/NHD_prepped_segmented/'
    segmentReaches(index=slurm, cpus_per_task=cpus_per_task, save_path=save_path)