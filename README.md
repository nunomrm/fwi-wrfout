# fwi-wrfout

fwi-wrfout (version 1.1) is a utility library that performs calculations of the [Fire Weather Index](https://cwfis.cfs.nrcan.gc.ca/background/summary/fwi) (FWI), among related operations. The main goal of this utility is to convert *wrfout* files (output files from the [WRF](https://www.mmm.ucar.edu/models/wrf) atmospheric model) into output [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) files containing the FWI index and its sub-indices (FFMC, DMC, DC, BUI, ISI). The source code for the FWI calculation is based on the [pyfwi project](https://code.google.com/archive/p/pyfwi/). 

The functions to compute FWI (```fwi_functions.py```) use the same calculations as the original *pyfwi* FWI functions. However, they were adapted for calculations in numpy arrays instead of individual numbers, as it was originally designed in *pyfwi*. This makes a significant difference in terms of efficiency when using data-heavy ```(LON,LAT,XTIME)``` numpy arrays, which is often the case of *wrfout* files (in the original version several for-loops would have to be performed).

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

# Release Notes

## Catalog of versions
* [v1.0](https://github.com/nunomrm/fwi-wrfout/tree/v1.0) (30/09/2023) [Deprecated]
* v1.1 (24/06/2024) [Current: read more in https://github.com/nunomrm/fwi-wrfout/releases/tag/v1.1]

## Potential future improvements and enhancements until version 2.0
- each of the functions in ```main_utils.py``` and in ```fwi_fuctions.py``` in their own individual files inside `utils/`
- adapt the chunk ``XTIME` size in the `xr.open_dataset` line of `generate_fwi_nc.py` depending on file size (e.g., very large wrfout files of ~8 GB need 20 chunks in order to )
- there's a few Python library incompatibilities with Python 3.11 (to fix in order to have a more generalizable and robust product)
- more generalizable initialization approach for the computing of FWI (instead of narrowing down months specific months, as described in the summary of changes made in v0.2), such as looking for the best day (in terms of weather conditions and time of the year) to initialize FWI automatically, depending on geographic location
- more and better documentation
- major refactorization

# Setting up

Clone the git repository:
```
git clone https://github.com/nunomrm/fwi-wrfout.git
```
Install required Python modules:
```
pip install -r requirements.txt
```

# Utilities description
The functions present in ```utils/main_utils.py``` contain the core utilities of fwi-wrfout. These utilities are listed and described below:
* ```fwi_idx```: filters the original time array indices into indices for calculating the FWI;
* ```extract_climate_vars```: imports the variables necessary to calculate the FWI and performs some calculations (e.g., obtain relative humidity, obtain hourly precipitation) and extracts a dictionary data structure with the attributes ```t2```, ```wind```,```rain_cumulative```, ```wind```, ```rh``` (for all time instants), ```t2_fwi```, ```wind_fwi```, ```rain_fwi```, ```rh_fwi``` (for FWI indices);
* ```compute_fwi```: allows for computation of FWI by calling ```fwi_functions``` (adapted from the *pyfwi* project to calculate FWI with numpy arrays), and attributes initial values of FWI sub-indices dependent on previous day values (FFMC, DMC and DC);
* ```calc_rh```: calculates relative humidity from pressure (```psfc```), the ratio of saturation mixture (```q2```), and temperature (```t2```);
* ```make_map```: creates the figure and draws the map with cartopy;
* ```plot_fwi_vars```: plots FWI variables, and currently adapted for FWI and n_days of FWI (in terms of custom colormaps and labels), not yet for FWI sub-indices.

# Run examples
## ```generate_fwi_nc.py```
Usage:
1. Go to ```scripts/```
2. Run ```python generate_fwi_nc.py```

Input: ```data/wrfout_files/wrfout_test.nc```

Output: ```output/nc_files/fwi_test.nc```

Description: Computes and generates a NetCDF FWI file with the same time, lon, lat format as the wrfout. Uses an ```xarray``` dataset to open the wrfout and extract its climate/weather information. For very large wrfouts (e.g. >8GB), using ```xarray.open_dataset``` for data extraction made the runtime decrease dramatically, due to parallelization routines with ```dask``` (```xarray``` dependency), compared to using ```netCDF4.Dataset```. In the utilities folder (```utils/```) has all tools for the computation of FWI. The ```compute_fwi``` function is the allows for computation of FWI, and which calls ```fwi_functions``` (which was adapted from the *pyfwi* project) and makes calculation of FWI with numpy arrays. Lastly, a NetCDF file is created, storing FWI and its' sub-indices.

## ```plots_fwi.py```
Usage:
1. Go to ```scripts/```
2. Run ```python plots_fwi.py```

Description: Creates three simple map plots with Matplotlib and Cartopy. Applies masking outside Portugal's Vila Real, Bragança and Guarda districts. These plots are saved in the ```PNG``` format in ```output/plots/```. Below we see two example plots.

Map of the FWI mean without masking (```fwi_mean_nomask.png```):

<img src="https://github.com/nunomrm/fwi-wrfout/blob/main/output/plots/fwi_mean_nomask.png" width="350"/>

Map of the FWI mean with masking (```fwi_mean.png```):

<img src="https://github.com/nunomrm/fwi-wrfout/blob/main/output/plots/fwi_mean.png" width="350"/>

