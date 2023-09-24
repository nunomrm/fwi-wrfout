#!/usr/bin/env python3
# -*- coding: utf-8 -*-

######################################
# directories and file names
######################################

utils_dir = '../utils/'
data_dir = '../output/' # in this case data is the computed fwi .nc file
fwi_dir = data_dir + 'nc_files/'
output_dir = '../output/plots/'
fname = 'fwi_test.nc'
plot_name1 = 'fwi_mean.png'
plot_name2 = 'ndays_extreme.png'

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
import geopandas as gpd
from shapely.geometry import Point
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

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
years_name = '2014'
var = 'fwi_mean'
fig_p = output_dir+plot_name2
cb_lims, delta_c = [[20, 50], 5]
fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
               var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
               colorbar_lims=cb_lims, colorbar_delta=delta_c)


danger_lvl = 'extreme'
var = 'ndays_'+danger_lvl
fig_p = output_dir+plot_name2
cb_lims, delta_c = [[0, 120], 20]
fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
               var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
               colorbar_lims=cb_lims, colorbar_delta=delta_c)


'''
from sklearn.metrics import mean_squared_error
utils_dir = '../utils/'
sys.path.append(utils_dir)
from custom_functions import *

fwi_dir = '../output_data/nc_files/'
output_dir = '../output_data/plots/fwi_vars_map/'
my_proj = ccrs.PlateCarree() # map projection
shapefile = gpd.read_file('../extras/CNTR_RG_01M_2020_4326.shp')
for i, name in enumerate(shapefile.CNTR_NAME):
    if name == 'Portugal':
        pt_polygon = shapefile.geometry[i]
fwi_lims = [38, 50] # high and very high FWI danger levels
d_lims = [-10, -6,  36.7, 42.3]
my_proj = ccrs.PlateCarree()
c = 0
files_list = os.listdir(fwi_dir)
for idx, i in enumerate(files_list):
    if 'MPI' in i:
        idx_f = idx
temp = files_list[idx_f]
files_list[idx_f] = files_list[0]
files_list[0] = temp
print(files_list)
        
for fname in files_list:
    spec_name = fname[fname.find('_'):-3]
    print(fname)
    ds_fwi = nc.Dataset(fwi_dir+fname)
    if c == 0:
        lat = (ds_fwi['XLAT'])
        lon = (ds_fwi['XLON'])
        Nx = lat.shape[0]
        Ny = lat.shape[1]
        mask = np.zeros((Nx, Ny))
        for i in range(Nx):
            for j in range(Ny):
                this_point = Point(lon[i, j], lat[i, j])
                mask[i, j] = pt_polygon.contains(this_point)
    c += 1

    # pre-process time
    nctime = ds_fwi.variables['time'][:]
    time = np.array(nctime)
    t_units = ds_fwi.variables['time'].units
    t_cal = ds_fwi.variables['time'].calendar
    tvalue = nc.num2date(nctime, units=t_units, calendar=t_cal)
    str_time = [i.strftime("%Y-%m-%d %H:%M") for i in tvalue]
    Nt = ds_fwi.dimensions['time'].size

    # get mean fwi and number of days for very high and extreme fwi conditions
    ds = xr.open_dataset(fwi_dir+fname).chunk({'time': 4})
    fwi = np.array(ds['FWI'])
    fwi_mean = np.nanmean(fwi,axis=0)
    fwi_mean[mask==False] = np.nan
    fwi_fields = {}
    fwi_fields['fwi_mean'] = fwi_mean
    where_fwi_extreme = np.zeros((len(fwi_lims),Nt,Nx,Ny)) 
    for i in range(2):
        where_fwi_extreme[i,:,:,:] = (fwi>=fwi_lims[i])
        fwi_extreme_ndays = np.zeros((len(fwi_lims),Nx,Ny)) 
    for i in range(2):
        fwi_extreme_ndays[i,:,:] = np.sum(where_fwi_extreme[i,:,:,:],axis=0)
    for i in range(2):
        nd = fwi_extreme_ndays[i,:,:]
        nd[mask==False] = np.nan
        fwi_extreme_ndays[i,:,:] = nd
    
    years_n = 20  # number of years to calculate ndays/year
    fwi_extreme_ndays = fwi_extreme_ndays/years_n
    fwi_fields['ndays_veryhigh'] = fwi_extreme_ndays[0,:,:]
    fwi_fields['ndays_extreme'] = fwi_extreme_ndays[1,:,:]
    if 'MPI' in fname:
        fwi_fields_mpi = {}
        fwi_fields_mpi['ndays_veryhigh'] = fwi_extreme_ndays[0,:,:]
        fwi_fields_mpi['ndays_extreme'] = fwi_extreme_ndays[1,:,:]
        fwi_fields_mpi['fwi'] = fwi
    if 'MPI' not in fname:
        fwi_fields['ndays_diff_veryhigh'] = fwi_fields['ndays_veryhigh'] - \
                                            fwi_fields_mpi['ndays_veryhigh']
        fwi_fields['ndays_diff_extreme'] = fwi_fields['ndays_extreme'] - \
                                            fwi_fields_mpi['ndays_extreme']
        if 'ERA5' in fname:  
            fwi_rmse = np.zeros((Nx,Ny))
            for i in range(Nx):
                for j in range(Ny):
                    fwi_rmse[i,j] = mean_squared_error(fwi[:,i,j],fwi_fields_mpi['fwi'][:,i,j],squared=False)
            fwi_rmse[mask==False] = np.nan
            fwi_fields['rmse'] = fwi_rmse
    ds.close()

    # plotting
    
    years_name = fname[len(fname)-find_nth(fname[::-1], '_', 2):][:4]+'-' +\
        fname[len(fname)-find_nth(fname[::-1], '_', 2):][5:-3]
    var = 'fwi_mean'
    fig_p = output_dir+'fwimean'+spec_name+'.png'
    cb_lims, delta_c = [[20, 50], 5]
    fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
    plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                   var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                   colorbar_lims=cb_lims, colorbar_delta=delta_c)
    
    danger_lvl = 'veryhigh'
    var = 'ndays_'+danger_lvl
    fig_p = output_dir+'ndays_'+danger_lvl+spec_name+'.png'
    cb_lims, delta_c = [[25,175], 20]
    fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
    plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                   var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                   colorbar_lims=cb_lims, colorbar_delta=delta_c)
    
    danger_lvl = 'extreme'
    var = 'ndays_'+danger_lvl
    fig_p = output_dir+'ndays_'+danger_lvl+spec_name+'.png'
    cb_lims, delta_c = [[0, 120], 20]
    fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])

    plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                   var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                   colorbar_lims=cb_lims, colorbar_delta=delta_c)
    
    if 'MPI' not in fname:
        danger_lvl = 'veryhigh'
        var = 'ndays_diff_'+danger_lvl
        fig_p = output_dir+'ndays_diff_'+danger_lvl+spec_name+'.png'
        cb_lims, delta_c = [[-40,40], 10]
        fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
        years_name = fname[len(fname)-find_nth(fname[::-1], '_', 2):][:4]+'-' +\
            fname[len(fname)-find_nth(fname[::-1], '_', 2):][5:-3]
        plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                       var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                       colorbar_lims=cb_lims, colorbar_delta=delta_c)
        
        danger_lvl = 'extreme'
        var = 'ndays_diff_'+danger_lvl
        fig_p = output_dir+'ndays_diff_'+danger_lvl+spec_name+'.png'
        cb_lims, delta_c = [[-40,40], 10]
        fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
        plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                       var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                       colorbar_lims=cb_lims, colorbar_delta=delta_c)
        
    if 'ERA5' in fname:
        var = 'rmse'
        fig_p = output_dir+'rmse_'+spec_name+'.png'
        cb_lims, delta_c = [[10,30], 5]
        fig, ax = make_map(my_proj, fnt_size=22, fig_size=[10, 14])
        plot_fwi_vars(fig=fig, ax=ax, fig_path=fig_p, domain_lims=d_lims, var=var,
                       var_dict=fwi_fields, lon=lon, lat=lat, years=years_name,
                       colorbar_lims=cb_lims, colorbar_delta=delta_c)
        
    ds.close()
'''