# -*- coding: utf-8 -*-

import numba as nb

def make_map(proj,fnt_size,fig_size):
    '''
    Creates a figure with matplotlib and initializes the map plot with cartopy 
    using lon/lat grid and label formatting configurations 

    Parameters
    ----------
    proj : cartopy.crs.PlateCarree
        Cartopy projection for visualization of the map.
    fnt_size : int
        Font size of lon (x) and lat (y) labels.
    fig_size : list
        List with width and height of the matplotlib figure.

    Returns
    -------
    fig : matplotlib.figure.Figure
        Figure.
    ax : cartopy.mpl.geoaxes.GeoAxes
        Plot axes.

    '''
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
    '''
    Extracts the time indices that match with the FWI calculation times (at 
    12:00, between Apr-Oct), and the new time array

    Parameters
    ----------
    t : numpy.ndarray
        Time array from the wrfout dataset.

    Returns
    -------
    i_t_fwi : numpy.ndarray
        Indices of the t array that correspond to the FWI calculation instants.
    t_fwi : numpy.ndarray
        Time array for FWI calculation.

    '''
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
    '''
    extracts climate/meteorologic variables that are needed for FWI calculation
    (from the wrfout dataset) and performs unit conversion calculation (wind 
    and temperature) or calculations from the rain variables

    Parameters
    ----------
    ds : xarray.core.dataset.Dataset
        wrfout dataset (opened with xarray).
    t : numpy.ndarray
        Time array from the wrfout dataset.
    **kwargs : str
        The only keyword argument available is the i_fwi array,
        which allows for data filtering of FWI calculation 
        times after extraction of the climate variables.

    Returns
    -------
    climate_vars : dictionary
        Dataset with the extracted variables.

    '''
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
        for i in range(1,len(i_fwi)):
            if t[i_fwi[i]].year == t[i_fwi[i-1]].year:
                rain_n[i,:,:] = rain[i_fwi[i],:,:] - rain[i_fwi[i-1],:,:]
            else:
                rain_n[i,:,:] = 0
        climate_vars['rain_fwi'] = rain_n
        # print('coomputing rain...')
        # @nb.njit(parallel=True)
        # def compute_rain_n(i_fwi, t, rain, rain_n):
        #     for i in nb.prange(1, len(i_fwi)):
        #         if t[i_fwi[i]].year == t[i_fwi[i-1]].year:
        #             rain_n[i,:,:] = rain[i_fwi[i],:,:] - rain[i_fwi[i-1]+1,:,:]
        # climate_vars['rain_fwi'] = compute_rain_n(i_fwi, t, rain, rain_n)
        # print('done')
        
    return climate_vars

def calc_rh(p, q, t):
    '''
    Calculates the relativity humidity with a known formula, from pressure,
    mixing ratio and temperature

    Parameters
    ----------
    p : numpy.ndarray
        Surface pressure.
    q : numpy.ndarray
        Water vapor mixing ratio at 2-meters height.
    t : numpy.ndarray
        Temperature  at 2-meters height.

    Returns
    -------
    rh : numpy.ndarray
        Relative humidity (units: %).

    '''
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
    
def plot_fwi_vars(fig, ax, fig_path, domain_lims, var, var_dict,
                   lon, lat, top_label, colorbar_lims, colorbar_delta):
    '''
    
    plots the FWI variable maps (TO DO: adapt for the FWI sub-indices
    in the "if" conditions)

    Returns
    -------
    None.

    '''
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import cartopy.feature as cf
    
    var_plot=var_dict[var]
    if 'fwi' in var:
        labelname = 'FWI'
        cmap = plt.cm.hot_r
    elif 'ndays' in var:
        if 'extreme' in var:
            cmap = plt.cm.YlOrRd
        elif 'veryhigh':
            cmap = plt.cm.Oranges
        labelname = 'n_days'
        
    cmaplist = [cmap(i) for i in range(cmap.N)]
    cb_ticks = np.arange(colorbar_lims[0],colorbar_lims[1]+colorbar_delta,colorbar_delta)
    bounds = np.arange(colorbar_lims[0],colorbar_lims[1]+0.01,colorbar_delta/4)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    ax.set_extent([domain_lims[0], domain_lims[1], domain_lims[2], domain_lims[3]])
    ax.add_feature(cf.COASTLINE)
    ax.add_feature(cf.BORDERS)
    ax.pcolormesh(lon,lat,var_plot,cmap=cmap,norm=norm,zorder=0)
    ax.text(-9.95,42.05,top_label,fontsize=24)
    ax2 = fig.add_axes([0.15, 0.06, 0.7, 0.03])
    cb = mpl.colorbar.ColorbarBase(ax=ax2, cmap=cmap, norm=norm,
                                    spacing='proportional', boundaries=bounds, 
                                    ticks=cb_ticks, orientation='horizontal',
                                    extend='both')
    cb.set_label(label=labelname, fontsize=24)
    cb.ax.tick_params(labelsize=24)
    plt.savefig(fig_path)
    plt.show()
    
def compute_fwi(t, i_fwi, Nx, Ny, lat, climate_vars, ffmc_0, dmc_0, dc_0, fwi_functions):
    '''
    Computes the Fire-Weather Index (FWI)
    
    '''
    import numpy as np
    
    t2, rh, wind, rain = [np.array(climate_vars['t2_fwi']),
                          np.array(climate_vars['rh_fwi']),
                          np.array(climate_vars['wind_fwi']),
                          np.array(climate_vars['rain_fwi'])]
    
    Nt = len(i_fwi)
    
    [ffmc, dmc, dc, isi, bui, fwi] = [np.zeros((Nt,Nx,Ny)),
                                      np.zeros((Nt,Nx,Ny)),
                                      np.zeros((Nt,Nx,Ny)),
                                      np.zeros((Nt,Nx,Ny)),
                                      np.zeros((Nt,Nx,Ny)),
                                      np.zeros((Nt,Nx,Ny))]
    y_prev = t[i_fwi][0].year 
    mm=0 # to delete
    for i, idx in enumerate(i_fwi):
        y = t[idx].year
        mm = t[idx].month
        if i == 0 or y!=y_prev:
            y_prev = t[idx].year   
            dc_prev = dc_0 * np.ones((Nx,Ny))
            dmc_prev = dmc_0 * np.ones((Nx,Ny))
            ffmc_prev = ffmc_0 * np.ones((Nx,Ny))
            dc[i,:,:] = dc_prev
            dmc[i,:,:] = dmc_prev
            ffmc[i,:,:] = ffmc_prev
        else:
            ffmc[i,:,:] = fwi_functions.FFMC(t2[i,:,:],
                                               rh[i,:,:],
                                               wind[i,:,:],
                                               rain[i,:,:],
                                               ffmc_prev)
            dmc[i,:,:] = fwi_functions.DMC(t2[i,:,:],
                                             rh[i,:,:],
                                             rain[i,:,:],
                                             dmc_prev,
                                             lat,
                                             mm)
            dc[i,:,:] = fwi_functions.DC(t2[i,:,:],
                                           rain[i,:,:],
                                           dc_prev,
                                           lat,
                                           mm)         
            ffmc_prev = ffmc[i,:,:]
            dmc_prev = dmc[i,:,:]
            dc_prev = dc[i,:,:]
        isi[i,:,:] = fwi_functions.ISI(wind[i,:,:],
                                       ffmc[i,:,:])
        bui[i,:,:] = fwi_functions.BUI(dmc[i,:,:],
                                       dc[i,:,:])
        fwi[i,:,:] = fwi_functions.FWI(isi[i,:,:],
                                       bui[i,:,:])
        
    return fwi, bui, isi, dc, dmc, ffmc