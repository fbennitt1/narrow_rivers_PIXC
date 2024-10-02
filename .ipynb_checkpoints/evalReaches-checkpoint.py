import numpy as np
import pandas as pd
import xarray as xr

from reaches import readNHD
from reaches import cut
from reaches import bitwiseMask
from reaches import makeGDF
from reaches import segmentReach
from reaches import getCoverage

## Setup
# SLurm indices for one HUC4 basin per HUC2 basin (randomly sampled)
indices = [120, 122, 13, 126, 33, 37, 134, 52, 54, 142, 73, 155, 83, 92, 95, 99, 110, 113]
# Order of all possible bins for plotting
order = ['(0, 10]', '(10, 20]', '(20, 30]', '(30, 40]', '(40, 50]',
         '(50, 60]', '(60, 70]', '(70, 80]', '(80, 90]', '(90, 100]',
         '(100, 150]', '(150, 200]', '(200, 500]', '(500, 1000]']
index = 109
pixc_path = '../data/small_rivers/mar_2024_ver_c/SWOT_L2_HR_PIXC_011_345_239R_20240227T044209_20240227T044220_PIC0_01.nc'
tile_name = pixc_path[-71:-3]

# Read prepped NHD
basin, huc4 = readNHD(index=index)

# Read in xarray
ds_PIXC = xr.open_mfdataset(paths=pixc_path, group='pixel_cloud', engine='h5netcdf')

# Make dict for legend labels
flags = ds_PIXC.classification.flag_meanings.split() # extract each flag meaning
codes = {idx:k for idx, k in enumerate(flags, start=1)}

# Make mask
mask = bitwiseMask(ds=ds_PIXC)

# Make PIXC
gdf_PIXC = makeGDF(ds=ds_PIXC, mask=mask, data_var='classification')

bounds_PIXC = gdf_PIXC.union_all().convex_hull
basin_crop = basin.clip(bounds_PIXC)
basin_crop = basin_crop[basin_crop.geometry.geometry.type=='LineString']

# Find width bins present in cropped hydrography
bins = list(np.unique(basin_crop.Bin))

# Make sorted list of bins present
ordered = []
for binn in order:
    if binn in bins:
        ordered.append(binn)
ordered

# Segment the reaches
basin_crop['segments'] = basin_crop.apply(func=segmentReach, axis=1)

# Get the coverage for each reach (at least one pixel/segment)
basin_crop['coverage'] = basin_crop.apply(func=getCoverage,
                                          args=(basin_crop.crs, gdf_PIXC), axis=1)

print(pd.DataFrame(basin_crop.coverage.value_counts().reset_index().sort_index()))

mean_cov = basin_crop.groupby('Bin')['coverage'].mean()
print(mean_cov)



# plt.scatter(x=ordered, y=coverage, c='k')
fig, ax = plt.subplots(figsize=(12,8))
plt.scatter(x=mean_cov.keys(), y=mean_cov.values, c='k')
plt.title('Mean percent of reach detected in ' + huc4 +' by \n'
          + tile_name)
plt.xticks(rotation=-45);

# Set number of ticks for x-axis
ax.set_xticks(range(len(ordered)))
# Set ticks labels for x-axis
ax.set_xticklabels(ordered)

plt.xlabel('Bankfull Width (m)')
plt.ylabel('Percent Detected')
plt.savefig(fname='./figures_for_unit_test/one_tile_'+ huc4 +'.png', bbox_inches='tight')