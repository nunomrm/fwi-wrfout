#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# directories and file names
######################################

utils_dir = '../utils/'
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
sys.path.append(utils_dir) # importing functions from ../utils/
from fwi_functions import *
from custom_functions import *
import numpy as np
import netCDF4 as nc
import xarray as xr

######################################
# extract and pre-process wrfout data
######################################

ds_nc = nc.Dataset(wrfout_dir+fname) # use only for time variable processing
nctime = ds_nc['XTIME']
time_input = np.array(nctime)
t_units = nctime.units
t_cal = nctime.calendar
t = nc.num2date(nctime,units=t_units,calendar=t_cal)
str_time = [i.strftime("%Y-%m-%d %H:%M") for i in t]
Nt = len(time_input)

ds = xr.open_dataset(wrfout_dir+fname).chunk({'XTIME':6})
lat = np.array(ds['XLAT'])
lon = np.array(ds['XLONG'])
Nx = lat.shape[0]
Ny = lat.shape[1]

# get FWI indices
i_fwi, time_fwi = fwi_idx(t)
Nt_fwi = len(i_fwi)

# extract vars
climate_vars = extract_climate_vars(ds, t, i_fwi=i_fwi)
t2, rh, wind, rain = [np.array(climate_vars['t2_fwi']),
                      np.array(climate_vars['rh_fwi']),
                      np.array(climate_vars['wind_fwi']),
                      np.array(climate_vars['rain_fwi'])]

ds.close()
ds_nc.close()

######################################
# compute FWI
######################################



ffmc_0, dmc_0, dc_0 = [85, 6, 15]
[ffmc, dmc, dc, isi, bui, fwi] = [np.zeros((Nt_fwi,Nx,Ny)),
                                  np.zeros((Nt_fwi,Nx,Ny)),
                                  np.zeros((Nt_fwi,Nx,Ny)),
                                  np.zeros((Nt_fwi,Nx,Ny)),
                                  np.zeros((Nt_fwi,Nx,Ny)),
                                  np.zeros((Nt_fwi,Nx,Ny))]

y_prev = t[0].year 
for i, idx in enumerate(i_fwi):
    mm = t[idx].month
    y = t[idx].year
    if i == 0 or y!=y_prev:
        y_prev = t[idx].year   
        dc_prev = dc_0 * np.ones((Nx,Ny))
        dmc_prev = dmc_0 * np.ones((Nx,Ny))
        ffmc_prev = ffmc_0 * np.ones((Nx,Ny))
        dc[i,:,:] = dc_prev
        dmc[i,:,:] = dmc_prev
        ffmc[i,:,:] = ffmc_prev
    else:
        ffmc[i,:,:] = FFMC(t2[i,:,:],
                       rh[i,:,:],
                       wind[i,:,:],
                       rain[i,:,:],
                       ffmc_prev)
        dmc[i,:,:] = DMC(t2[i,:,:],
                         rh[i,:,:],
                         rain[i,:,:],
                         dmc_prev,
                         lat,
                         mm)
        dc[i,:,:] = DC(t2[i,:,:],
                       rain[i,:,:],
                       dc_prev,
                       lat,
                       mm)         
        ffmc_prev = ffmc[i,:,:]
        dmc_prev = dmc[i,:,:]
        dc_prev = dc[i,:,:]
    isi[i,:,:] = ISI(wind[i,:,:],
                     ffmc[i,:,:])
    bui[i,:,:] = BUI(dmc[i,:,:],
                    dc[i,:,:])
    fwi[i,:,:] = FWI(isi[i,:,:],
                    bui[i,:,:])

ds.close()

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
