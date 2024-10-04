import numpy as np
import os
import pandas as pd
import xarray as xr

from reaches import readNHD
from reaches import cut
from reaches import bitwiseMask
from reaches import makeGDF
from reaches import segmentReach
from reaches import getCoverage

## Setup
# Slurm indices for one HUC4 basin per HUC2 basin (randomly sampled)
# indices = [120, 122, 13, 126, 33, 37, 134, 52, 54, 142, 73, 155, 83, 92, 95, 99, 110, 113]
# Order of all possible bins for plotting
order = ['(0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]',
         '(50, 60]', '(60, 70]', '(70, 80]', '(80, 90]', '(90, 100]',
         '(100, 150]', '(150, 200]', '(200, 500]', '(500, 1000]']

data_path = '/nas/cee-water/cjgleason/fiona/data/small_rivers/mar_2024_ver_c/'

## For testing
# Set to 0 or 1
i = 0 # 1

if i == 0:
    index = 4 # HUC4 0108, Connecticut
    pixc_path = os.path.join(data_path, 'leaf_off/SWOT_L2_HR_PIXC_014_341_229R_20240429T152954_20240429T153005_PIC0_01.nc')
elif i == 1:
    index = 109 # HUC4 1711, Snoqualmie
    pixc_path = os.path.join(data_path, 'SWOT_L2_HR_PIXC_011_345_239R_20240227T044209_20240227T044220_PIC0_01.nc')

tile_name = pixc_path[-71:-3]

## Prepare data
# Read prepped NHD
basin, huc4 = readNHD(index=index)

# Read in PIXC granule
ds_PIXC = xr.open_mfdataset(paths=pixc_path, group='pixel_cloud', engine='h5netcdf')

# Make mask for PIXC quality flags
mask = bitwiseMask(ds=ds_PIXC)

# Make PIXC GDF
gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_var='classification')

# Get bounds of PIXC granule
bounds_PIXC = gdf_PIXC.union_all().convex_hull

# Crop NHD basin to bounds of PIXC granule
basin_crop = basin.clip(bounds_PIXC)

# Drop all geometries that aren't LineStrings
# (If a river flows across the PIXC boundary then back in, becomes
# MultiLineString and needs to be dropped for now).
basin_crop = basin_crop[basin_crop.geometry.geometry.type=='LineString']

# Segment the reaches
basin_crop['segments'] = basin_crop.apply(func=segmentReach, axis=1)

## Analyze coverage
# Set number of pixels per segment needed for detection
num = 1

# Get the coverage for each reach (at least one pixel/segment)
basin_crop['coverage'] = basin_crop.apply(func=getCoverage,
                                          args=(basin_crop.crs,
                                                gdf_PIXC, num), axis=1)

print(pd.DataFrame(basin_crop.coverage.value_counts().reset_index().sort_index()))

# Get mean coverage by width bin
mean_cov = basin_crop.groupby('Bin')['coverage'].mean()
print(mean_cov)