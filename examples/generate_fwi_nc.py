#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# directories and file names
######################################

data_dir = '../data/'
wrfout_dir = data_dir + 'wrfout_files/'
output_dir = '../output/'
out_nc_files = output_dir + 'nc_files/'
fname = 'wrfout_test.nc'
fname_out = 'fwi_test.nc'

######################################
# import required python modules
######################################

import sys, os
sys.path.append('../utils/') # allowing custom module imports from the utils folder
from fwi_functions import *
from custom_functions import *
import fwi_functions
import numpy as np
import netCDF4 as nc
import xarray as xr

######################################
# extract and pre-process wrfout data
######################################

# open dataset with xarray
ds = xr.open_dataset(wrfout_dir+fname, decode_times=False).chunk({'XTIME':6})

# extract latitude and longitude
lat = np.array(ds['XLAT'])
lon = np.array(ds['XLONG'])
Nx = lat.shape[0]
Ny = lat.shape[1]

# process time
time_input = np.array(ds['XTIME'])
t_units = ds['XTIME'].units
t_cal = ds['XTIME'].calendar
t = nc.num2date(time_input,units=t_units,calendar=t_cal)

# get FWI time indices
i_fwi, time_fwi = fwi_idx(t)
Nt_fwi = len(i_fwi)

# extract vars
climate_vars = extract_climate_vars(ds, t, i_fwi=i_fwi)

ds.close()

######################################
# compute FWI
######################################

ffmc_0, dmc_0, dc_0 = [85, 6, 15]

fwi, bui, isi, dc, dmc, ffmc = compute_fwi(t, i_fwi, Nx,
                                           Ny, lat, 
                                           climate_vars, 
                                           ffmc_0, dmc_0, 
                                           dc_0, fwi_functions)

######################################
# write to netCDF
######################################

ds = nc.Dataset(out_nc_files+fname_out,'w',format='NETCDF4')

# create dimensions
time_dim = ds.createDimension('XTIME', Nt_fwi)
south_north = ds.createDimension('south_north', Nx)
west_east = ds.createDimension('west_east', Ny)

# create variables
time_nc = ds.createVariable('XTIME', 'f4', ('XTIME',))
time_nc.units = t_units
time_nc.calendar = t_cal
lats = ds.createVariable('XLAT', 'f4', ('south_north', 'west_east'))
lons = ds.createVariable('XLONG', 'f4', ('south_north', 'west_east'))
fwi_nc = ds.createVariable('FWI', 'f4', ('XTIME', 'south_north', 'west_east'))
ffmc_nc = ds.createVariable('FFMC', 'f4', ('XTIME', 'south_north', 'west_east'))
dmc_nc = ds.createVariable('DMC', 'f4', ('XTIME', 'south_north', 'west_east'))
dc_nc = ds.createVariable('DC', 'f4', ('XTIME', 'south_north', 'west_east'))
isi_nc = ds.createVariable('ISI', 'f4', ('XTIME', 'south_north', 'west_east'))
bui_nc = ds.createVariable('BUI', 'f4', ('XTIME', 'south_north', 'west_east'))

# assign data to variables
time = time_input[i_fwi]
lats[:]=lat
lons[:]=lon
time_nc[:]=time
fwi_nc[:,:,:]=fwi
ffmc_nc[:,:,:]=ffmc
dmc_nc[:,:,:]=dmc
dc_nc[:,:,:]=dc
isi_nc[:,:,:]=isi
bui_nc[:,:,:]=bui
ds.close()
