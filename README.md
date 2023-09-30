# fwi-wrfout

fwi-wrfout (version 1.0) is a utility library that performs calculations of the [Fire Weather Index](https://cwfis.cfs.nrcan.gc.ca/background/summary/fwi) (FWI), among related operations. The main goal of this utility is to convert wrfout files into output [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) files containing the FWI index and its sub-indices (FFMC, DMC, DC, BUI, ISI). The source code for the FWI calculation is based on the [pyfwi project](https://code.google.com/archive/p/pyfwi/). 

The functions to compute FWI (```fwi_functions.py```) use the same calculations as the original *pyfwi* FWI functions. However, they were adapted for calculations in numpy arrays, which is much more efficient than calculation of individual numbers, as it was originally designed in *pyfwi*. This makes a significant difference in terms of efficiency when using data-heavy ```(LON,LAT,XTIME)``` numpy arrays, which is often the case of *wrfout* files (*wrfouts* are NetCDF output files from the [WRF](https://www.mmm.ucar.edu/models/wrf)model).

fwi-wrfout is distributed under the 3-Clause BSD License (see the LICENSE.txt file).

# Prerequisites

[Python 3.10](https://www.python.org/downloads/release/python-3100/) or posterior versions are recommended.

Tests were performed with these python modules:
* numpy: 1.23.5 
* netCDF4: 1.6.3 
* matplotlib: 3.7.1 
* geopandas: 0.12.2 
* xarray: 2023.4.2 
* cartopy: 0.21.1 
* shapely: 2.0.1 

# Setting up

Clone the git repository:
```
git clone https://github.com/nunomrm/fwi-wrfout.git
```
Install required Python modules:
```
pip install -r requirements.txt
```

# Run examples
```
generate_fwi_nc.py
```
Usage:
1. Go to ```examples/```
2. Run ```python generate_fwi_nc.py```

Input: ```data/wrfout_files/wrfout_test.nc```
Output: ```output/nc_files/fwi_test.nc```
Description: Computes and generates a NetCDF FWI file with the same time, lon, lat format as the wrfout. Uses an ```xarray``` dataset to open the wrfout and extract its climate/weather information. For very large wrfouts (e.g. >8GB), using ```xarray.open_dataset``` for data extraction made the runtime decrease dramatically, due to parallelization routines with ```dask``` (```xarray``` dependency), compared to using ```netCDF4.Dataset```. In the utilities folder (```utils/```) has all tools for the computation of FWI. The ```compute_fwi``` function is the allows for computation of FWI, and which calls ```fwi_functions``` (which was adapted from the *pyfwi* project) and makes calculation of FWI with numpy arrays. Lastly, a NetCDF4 file is created, storing FWI and its' sub-indices.

```
plots_fwi.py
```
Usage:
1. Go to ```examples/```
2. Run ```python plots_fwi.py```

Description: Creates three simple map plots with Matplotlib and Cartopy. Applies masking outside Portugal's Vila Real, Bragança and Guarda districts. These plots are saved in the ```PNG``` format in ````output/plots/```. Below we see two example plots.


Map of the FWI mean without masking (```fwi_mean_nomask.png```):

<img src="https://github.com/nunomrm/fwi-wrfout/blob/main/output/plots/fwi_mean_nomask.png" width="350"/>

Map of the FWI mean with masking (```fwi_mean.png```):

<img src="https://github.com/nunomrm/fwi-wrfout/blob/main/output/plots/fwi_mean.png" width="350"/>

