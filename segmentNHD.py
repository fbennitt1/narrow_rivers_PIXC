import os

import geopandas as gpd
import pandas as pd
from pandarallel import pandarallel

from reaches import readNHD
from reaches import cut
from reaches import segmentReach

'''
'''
# Get slurm job index
slurm = int(os.environ['SLURM_ARRAY_TASK_ID'])

## Prepare data
# Read prepped NHD
basin, huc4, huc2 = readNHD(index=slurm)

pandarallel.initialize()

# Segment the reaches
basin[['segments', 'failed']] = basin.parallel_apply(user_defined_function=segmentReach, axis=1).apply(pd.Series)

# Keep only reaches that were sucessfully segmented
basin = basin[basin['failed'] == 0]

# Drop failure column
basin = basin.drop(columns='failed')

# Explode the segments to give each own row
basin =  basin.explode(column='segments', index_parts=True)

# Drop original reach geometry column, set segments as active geometry
basin = basin.drop(columns='geometry').set_geometry('segments').set_crs(crs=basin.crs)

# Write out
# Set write filepath
save_path = '../narrow_rivers_PIXC_data/NHD_prepped_segmented/'
save_path = os.path.join(save_path, huc2)
save_file = huc4 + '_prepped_segmented.parquet'

#Write out gdf as parquet file
if not os.path.isdir(save_path):
    os.makedirs(save_path)
basin.to_parquet(os.path.join(save_path, save_file))