# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7

@author: monte
"""

import datetime
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
import geopandas as gpd
from shapely.geometry import Point
import xarray as xr
import sys
import os
src_dir = '../utils'
sys.path.append(src_dir)
import cartopy.crs as ccrs


fwi_dir = '../output/nc_files/'
output_dir = '../output/plots/fwi_vars_map/'
my_proj = ccrs.PlateCarree() # map projection
shapefile = gpd.read_file('../data/shapefiles/distritos.shp')
fwi_lims = [38, 50] # high and very high FWI danger levels
d_lims = [-10, -6,  36.7, 42.3]

fname = 'fwi_test.nc'

# defining regions of Portugal
reg_pt = ['Litoral Norte','Interior Norte e Centro','Litoral Centro','Sul']
colors_reg_pt = ['b','g','orange','r']
regions = {}
regions[reg_pt[0]] = ['Viana do Castelo','Braga','Porto']
regions[reg_pt[1]] = ['Bragança','Vila Real','Guarda','Viseu','Castelo Branco']
regions[reg_pt[2]] = ['Aveiro','Coimbra','Leiria','Lisboa']
regions[reg_pt[3]] = ['Santarém','Setúbal','Beja','Portalegre','Évora','Faro'] 
c = 0
fwi_lims = [38, 50] # very high and extreme FWI danger levels
    
        
fwi_districts = {}
ds = nc.Dataset(fwi_dir+fname)

time_input = np.array(ds['XTIME'])
t_units = ds['XTIME'].units
t_cal = ds['XTIME'].calendar
t = nc.num2date(time_input,units=t_units,calendar=t_cal)
N_t = len(t)
lat = (ds['XLAT'])
lon = (ds['XLONG'])
Nx = lat.shape[0]
Ny = lat.shape[1]
spec_name = fname[fname.find('_'):-3]

# get mean fwi and number of days for very high and extreme fwi conditions
fwi = np.array(ds['FWI'])
fwi_timeseries_reg = {}
for r in reg_pt:
    mask = np.zeros((Nx, Ny))
    s=0
    for i_s, name in enumerate(shapefile.NAME_1):
        if name in regions[r]:
            for i in range(Nx):
                for j in range(Ny):
                    this_point = Point(lon[i, j], lat[i, j])
                    if shapefile.geometry[i_s].contains(this_point):
                        mask[i, j] = 1

    # ds = xr.open_dataset(fwi_dir+fname).chunk({'XTIME': 4})
    fwi = np.array(ds['FWI'])
    fwi_t = np.zeros((N_t))
    for i in range(N_t):
        fwi_i = fwi[i,:,:]
        fwi_i[mask==False] = np.nan
        fwi_i = np.nanmean(fwi_i)
        fwi_t[i] = fwi_i
    fwi_timeseries_reg[r] = fwi_t

if 'MPI' in fname:
    tag = 'MPI'
elif 'ssp245' in fname:
    tag = 'SSP2-4.5'
elif 'ssp370' in fname:
    tag = 'SSP3-7.0'
elif 'ssp585' in fname:
    tag = 'SSP5-8.5'
elif 'ERA5' in fname:
    tag = 'ERA5'

fig, ax = plt.subplots(4,1, figsize=(16,9), sharey=True)
for i, r in enumerate(reg_pt):
    fwi_t = fwi_timeseries_reg[r]
    ax[i].plot_date(t,fwi_t,markersize=.7,
                    color=colors_reg_pt[i],
                    linestyle='-',
                    marker=None)
    ax[i].set_xlim(t[0],t[-1])
    ax[i].set_ylim(0,100)
    ax[i].set_yticks(range(0,100,25))
    ax[i].grid()
ax[i].set_xlabel('time (years)')
plt.rcParams.update({'font.size': 22})
plt.subplots_adjust(hspace = 0.01)
plt.show()

