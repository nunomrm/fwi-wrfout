# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 17:34:46 2023

@author: monte
"""

def make_map(proj,fnt_size,fig_size):
    from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(fig_size[0],fig_size[1]),
                            subplot_kw=dict(projection=proj))
    gl = ax.gridlines(draw_labels=True)
    gl.top_labels = gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': fnt_size}
    gl.ylabel_style = {'size': fnt_size}
    return fig, ax

def fwi_idx(t):
    import numpy as np
    N_t = len(t)
    i_t_fwi = []
    t_fwi = np.array([])
    for i in range(N_t):
        mm, hr = [t[i].month, t[i].hour]
        if mm>=4 and mm<=10 and hr==12:
            i_t_fwi.append(i)
            t_fwi = np.append(t_fwi, t[i])
    return i_t_fwi, t_fwi

def extract_climate_vars(ds, t, **kwargs):
    import xarray as xr
    import numpy as np
    import netCDF4 as nc
    
    climate_vars = {}
    climate_vars['t2'] = xr.DataArray(ds['T2']) - 273.15
    climate_vars['q2'] = xr.DataArray(ds['Q2'])
    u10 = xr.DataArray(ds['U10']) * 3.6
    v10 = xr.DataArray(ds['V10']) * 3.6
    climate_vars['wind'] = np.sqrt(u10**2 + v10**2)
    climate_vars['psfc'] = xr.DataArray(ds['PSFC'])
    climate_vars['rh'] = calc_rh(climate_vars['psfc'],
                                 climate_vars['q2'],
                                 climate_vars['t2']+273.15)
    rnc = xr.DataArray(ds['RAINNC'])
    rc = xr.DataArray(ds['RAINC'])
    i_rc = xr.DataArray(ds['I_RAINC'])
    i_rnc = xr.DataArray(ds['I_RAINNC'])
    rain =  (100*i_rnc+rnc)+(100*i_rc+rc)
    climate_vars['rain_cumulative'] = rain
    
    if kwargs!={} and 'i_fwi' in kwargs:
        i_fwi = kwargs['i_fwi']
        climate_vars['t2_fwi'] = climate_vars['t2'][i_fwi,:,:]
        climate_vars['q2_fwi'] = climate_vars['q2'][i_fwi,:,:]
        climate_vars['psfc_fwi'] = climate_vars['psfc'][i_fwi,:,:]
        climate_vars['rh_fwi'] = climate_vars['rh'][i_fwi,:,:]
        climate_vars['wind_fwi'] = climate_vars['wind'][i_fwi,:,:]
        rain_n = xr.zeros_like(climate_vars['t2_fwi'])
        for i in range(len(i_fwi)):
            if i > 0 and t[i_fwi[i]].year == t[i_fwi[i-1]].year:
                rain_n[i,:,:] = rain[i_fwi[i],:,:] - rain[i_fwi[i-1]+1,:,:]
        climate_vars['rain_fwi'] = rain_n
    
    if kwargs!={} and 'i_summer' in kwargs:
        i_s = kwargs['i_summer']
        climate_vars['t2_summer'] = climate_vars['t2'][i_s,:,:]
        climate_vars['q2_summer'] = climate_vars['q2'][i_s,:,:]
        climate_vars['psfc_summer'] = climate_vars['psfc'][i_s,:,:]
        climate_vars['rh_summer'] = climate_vars['rh'][i_s,:,:]
        climate_vars['wind_summer'] = climate_vars['wind'][i_s,:,:]
        rain_n = rain[i_s,:,:]
        for i in range(len(i_s)):
            if t[i_s[i]].year != t[i_s[i-1]].year :
                rain_n[-1,:,:] = rain_n[-1,:,:] - (rain[i_s[i],:,:] - rain[i_s[i-1],:,:])
        climate_vars['rain_cumulative_summer'] = rain_n
    
    return climate_vars


def calc_rh(p,q,t):
    import numpy as np
    
    # constants
    svp1 = 611.2
    svp2 = 17.67
    svp3 = 29.65
    svpt0 = 273.15
    eps = .622   
    # RH formula
    rh = 100 * (p*q/(q*(1.-eps) + eps))/(svp1*np.exp(svp2*(t-svpt0)/(t-svp3)))
    
    return rh

def find_nth(haystack, needle, n):
    start = haystack.find(needle)
    while start >= 0 and n > 1:
        start = haystack.find(needle, start+1)
        n -= 1
    return start
        
def plot_clim_vars(fig, ax, fig_path, domain_lims, var, var_dict,
                   lon, lat, years, colorbar_lims, colorbar_delta):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.feature as cf
    
    var_plot=var_dict[var]
    if 't2' in var:
        labelname = 'T-2 m ('+u'\N{DEGREE SIGN}'+'C)'
        cmap = plt.cm.jet
    elif 'rh' in var:
        labelname = 'RH (%)'
        cmap = plt.cm.Blues
    elif 'rain' in var:
        labelname = 'Precipitation (mm)'
        cmap = plt.cm.YlGnBu
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cb_ticks = np.arange(colorbar_lims[0],colorbar_lims[1]+colorbar_delta,colorbar_delta)
    bounds = np.arange(colorbar_lims[0],colorbar_lims[1]+0.01,colorbar_delta/4)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax.set_extent([domain_lims[0], domain_lims[1], domain_lims[2], domain_lims[3]])
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.BORDERS)
    ax.pcolormesh(lon,lat,var_plot,cmap=cmap,norm=norm,zorder=0)
    ax.text(-9.95,42.05,years,fontsize=24)
    ax2 = fig.add_axes([0.15, 0.06, 0.7, 0.03])
    cb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cmap, norm=norm,
                                    spacing='proportional', boundaries=bounds, 
                                    ticks=cb_ticks, orientation='horizontal',
                                    extend='both')
    cb.set_label(label=labelname, fontsize=24)
    cb.ax.tick_params(labelsize=24)
    plt.savefig(fig_path)
    plt.show()
    
def plot_fwi_vars(fig, ax, fig_path, domain_lims, var, var_dict,
                   lon, lat, years, colorbar_lims, colorbar_delta):
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.feature as cf
    
    var_plot=var_dict[var]
    if 'fwi' in var:
        labelname = 'FWI'
        cmap = plt.cm.hot_r
        labelname = 'FWI'
    elif 'ndays' in var and not 'diff' in var:
        if 'extreme' in var:
            cmap = plt.cm.YlOrRd
        elif 'veryhigh':
            cmap = plt.cm.Oranges
        labelname = 'n_days/year'
    elif 'diff' in var:
        labelname = 'n_days/year anomaly'
        cmap = plt.cm.bwr
    elif 'rmse' in var:
        labelname = 'RMSE'
        cmap = plt.cm.plasma
        
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cb_ticks = np.arange(colorbar_lims[0],colorbar_lims[1]+colorbar_delta,colorbar_delta)
    bounds = np.arange(colorbar_lims[0],colorbar_lims[1]+0.01,colorbar_delta/4)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax.set_extent([domain_lims[0], domain_lims[1], domain_lims[2], domain_lims[3]])
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.BORDERS)
    ax.pcolormesh(lon,lat,var_plot,cmap=cmap,norm=norm,zorder=0)
    ax.text(-9.95,42.05,years,fontsize=24)
    ax2 = fig.add_axes([0.15, 0.06, 0.7, 0.03])
    cb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cmap, norm=norm,
                                    spacing='proportional', boundaries=bounds, 
                                    ticks=cb_ticks, orientation='horizontal',
                                    extend='both')
    cb.set_label(label=labelname, fontsize=24)
    cb.ax.tick_params(labelsize=24)
    plt.savefig(fig_path)
    plt.show()