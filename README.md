# fwi-wrfout

A set of python scripts that calculate the Fire Weather Index (FWI), from the Canadian Forest Fire Weather Index (FWI) System (more info about FWI: https://cwfis.cfs.nrcan.gc.ca/background/summary/fwi), in an efficient way in order to convert wrfout files into output netcdf files with the FWI index and its sub-indices (FFMC, DMC, DC, BUI, ISI). The source code for the FWI calculation is based on the pyfwi project (https://code.google.com/archive/p/pyfwi/).

fwi-wrf is distributed under the 3-Clause BSD License (see the LICENSE.txt file).

# Prerequisites

Python 3.10 or posterior versions are recommended.

Tests were performed with these python modules:
* numpy: 1.23.5 
* netCDF4: 1.6.3 
* matplotlib: 3.7.1 
* geopandas: 0.12.2 
* xarray: 2023.4.2 
* cartopy: 0.21.1 
* shapely: 2.0.1 

# Installation

Clone the git repository:
```
git clone https://github.com/nunomrm/fwi-wrfout.git
```
Install required Python modules:
```
pip install -r requirements.txt
```

# Examples
```
generate_fwi_nc.py
```
Usage:
1. Go to ```examples/```
2. Run ```python generate_fwi_nc.py```

Input: ```data/wrfout_files/wrfout_test.nc```
Output: ```output/nc_files/fwi_test.nc```
Description: Computes and generates a NetCDF FWI file with the same time, lon, lat format as the wrfout. Uses an ```xarray``` dataset to open the wrfout and extract its climate/weather information. For very large wrfouts (e.g. >8GB), using ```xarray.open_dataset``` for data extraction made the runtime decrease dramatically, due to parallelization routines with ```dask``` (```xarray``` dependency), compared to using ```netCDF4.Dataset```. Functions in ```custom_functions``` contain all required tools for the computation of FWI. Emphasys is given to ```compute_fwi```, which calls ```fwi_functions``` (from the pyfwi project, mentioned above) and makes calculation of FWI with numpy arrays. Lastly, a NetCDF4 file is created, storing FWI and its' sub-indices.

```
plots_fwi.py
```
Usage:
1. Go to ```examples/```
2. Run ```python generate_fwi_nc.py```

Description: Makes (and saves in ```.PNG``` format) three simple plots with Matplotlib and Cartopy that can be seen below.

```

```





