# -*- coding: utf-8 -*-

######################################
# directories and file names
######################################

data_dir = '../output/' # in this case data is the computed fwi .nc file
fwi_dir = data_dir + 'nc_files/'
output_dir = '../output/plots/'
fname = 'fwi_test.nc'
plot_name1 = 'fwi_mean.png'
plot_name2 = 'ndays_extreme.png'

######################################
# import required python modules
######################################


import sys
sys.path.append('../utils/') # allowing custom module imports from the utils folder
from custom_functions import *
import numpy as np
import xarray as xr
import geopandas as gpd
from shapely.geometry import Point
import cartopy.crs as ccrs

######################################
# definition of map limits and FWI thresholds
######################################

fwi_lims = [38, 50] # high and very high FWI danger levels
d_lims = [-10, -6,  36.7, 42.3] # geographical domain limits

######################################
# reading FWI netcdf files
######################################

ds = xr.open_dataset(fwi_dir+fname).chunk({'XTIME':4})
lat = np.array(ds['XLAT'])
lon = np.array(ds['XLONG'])
Nx = lat.shape[0]
Ny = lat.shape[1]
t = ds['XTIME']
Nt = len(t)
fwi = np.array(ds['FWI'])

######################################
# creating mask of study region
######################################

shapefile = gpd.read_file('../data/extras/distritos.shp')
# "Interior Norte e Centro" region's districts (for masking in the shapefile)
region = ['Bragança','Vila Real','Guarda','Viseu','Castelo Branco']
mask = np.zeros((Nx, Ny))
for i_s, name in enumerate(shapefile.NAME_1):
    if name in region:
        for i in range(Nx):
            for j in range(Ny):
                this_point = Point(lon[i, j], lat[i, j])
                if shapefile.geometry[i_s].contains(this_point):
                    mask[i, j] = 1

######################################
# Mean FWI and n_days calculation
# with masking implementation
######################################

fwi_mean = np.nanmean(fwi,axis=0)
fwi_mean[mask==False] = np.nan

where_fwi_threshold = np.zeros((len(fwi_lims),Nt,Nx,Ny)) 
for i in range(2): # goes through very high and extreme FWI thresholds
    where_fwi_threshold[i,:,:,:] = (fwi>=fwi_lims[i])
    fwi_threshold_ndays = np.zeros((len(fwi_lims),Nx,Ny)) 
for i in range(2):
    fwi_threshold_ndays[i,:,:] = np.sum(where_fwi_threshold[i,:,:,:],axis=0)
for i in range(2):
    nd = fwi_threshold_ndays[i,:,:]
    nd[mask==False] = np.nan
    fwi_threshold_ndays[i,:,:] = nd
    

fwi_fields = {} # dictionary which stores arrays of the calculations on FWI
fwi_fields['fwi_mean'] = fwi_mean
fwi_fields['ndays_veryhigh'] = fwi_threshold_ndays[0,:,:]
fwi_fields['ndays_extreme'] = fwi_threshold_ndays[1,:,:]

######################################
# plotting
######################################

my_proj = ccrs.PlateCarree() # map projection
label_name = 'FWI (mean)'
var = 'fwi_mean'
fig_p = output_dir+plot_name2
cb_lims, delta_c = [[15, 40], 5] # colorbar limits and its' delta
fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
               var_dict=fwi_fields, lon=lon, lat=lat, top_label=label_name,
               colorbar_lims=cb_lims, colorbar_delta=delta_c)

label_name = 'FWI ≥ 38'
danger_lvl = 'veryhigh'
var = 'ndays_'+danger_lvl
fig_p = output_dir+plot_name2
cb_lims, delta_c = [[0, 20], 4]
fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
               var_dict=fwi_fields, lon=lon, lat=lat, top_label=label_name,
               colorbar_lims=cb_lims, colorbar_delta=delta_c)

label_name = 'FWI ≥ 50'
danger_lvl = 'extreme'
var = 'ndays_'+danger_lvl
fig_p = output_dir+plot_name2
cb_lims, delta_c = [[0, 10], 2]
fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
               var_dict=fwi_fields, lon=lon, lat=lat, top_label=label_name,
               colorbar_lims=cb_lims, colorbar_delta=delta_c)
