#%%
def compute_mean_atmos_vort(date):
    t0        = date - delta_day + delta_hr
    t1        = date
    t_list    = pd.date_range(t0, t1, freq='H')
    
    vort      = ds_sector.sel(time=t_list).vo.values
    vort_mean = vort.mean(axis = 0)
 
    #process_name = current_process().name
    #print(f"Process name: {process_name}.")

    return vort_mean
#%%
def compute_mean_ice_vort(date):
    
    fname = dir_ice + f'{year}/{str(date.month).zfill(2)}/ice_drift_sh_ease2-750_cdr-v1p0_24h-{year}{str(date.month).zfill(2)}{str(date.day).zfill(2)}1200.nc'
    try:
        ds        = xr.open_dataset(fname)
        flag      = ds.status_flag.values[0]
        flag      = (flag>=flag_threshold)
        #extract drift values: original units in kilometers per 24hrs, modified units in meters per second
        dX        = ds.dX[0].values
        dY        = ds.dY[0].values
        dX[~flag] = np.nan
        dY[~flag] = np.nan
        U         = ((dX)*1000)/(24*60*60)
        V         = ((dY)*1000)/(24*60*60)
        vort      = mpcalc.vorticity(u=U*units('m/s'), v=V*units('m/s'), dx=75000*units('m'), dy=-75000*units('m'))

    except(OSError):
        print('[X] OSError "File not found": ' + fname)
        blank    = np.empty((144,144))
        blank[:] = np.nan
        vort     = blank

    #process_name = current_process().name
    #print(f"Process name: {process_name}.")

    return vort
#%%
def swath_remapping(field):
    
    in_data = field
    #resample consideration radius (units km to meters)
    sample_radius_m = sample_radius*1000 
    # input grid in swath form
    in_grid  = pr.geometry.SwathDefinition(lons=pr.utils.wrap_longitudes(in_lon), lats=in_lat)
    # output grid in defined projection
    out_grid = pr.geometry.SwathDefinition(lons=out_lon, lats=out_lat)
    # remap input_grid onto output_grid
    out_data = pr.kd_tree.resample_nearest(in_grid, in_data, out_grid, radius_of_influence=sample_radius_m, fill_value=np.nan)

    return out_data
#%%
def plot_thing1(lons, lats, atmos_vort, ice_vort, atmos_div, ice_div):

   #generate custom cmap1
    boundaries = np.arange(-1, 1, 0.05)
    cmap_custom = plt.cm.get_cmap('seismic',len(boundaries))
    colors = list(cmap_custom(np.arange(len(boundaries))))
    colors[0:int(-1+len(boundaries)/2)] = colors[1:int(len(boundaries)/2)]
    colors[int(1+len(boundaries)/2):int(len(boundaries))] = colors[int(len(boundaries)/2):int(-1+len(boundaries))]
    colors[int(-1+len(boundaries)/2)] = [1,1,1,1]
    colors[int(len(boundaries)/2)] = [1,1,1,1]
    cmap1 = mpl.colors.ListedColormap(colors, "k") #can replace index color eg: colors[0] = "white"
    cmap1.set_over(colors[-1]) # set over-color to last color of list 
    pc   = ccrs.PlateCarree()
    sps  = ccrs.SouthPolarStereo(central_longitude=0)
    e2n  = ccrs.LambertAzimuthalEqualArea(central_latitude=-90.0)
    orth = ccrs.Orthographic(central_longitude=0.0, central_latitude=-90.0)
    extent = [-180., 70., -60., -60.]

    
    fig = plt.figure(figsize=(20, 16))

    scale = 10000
    ax1 = plt.subplot(221, projection=sps)
    ax1.coastlines()
    ax1.gridlines()
    field = ax1.pcolormesh(lons, lats, atmos_vort*scale, transform=pc, cmap=cmap1, vmin = -1, vmax=1)
    ax1.set_extent(extent, crs=pc)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    order   = str(int(np.log10(scale)))
    cbar    = fig.colorbar(field, ax=ax1, location ='bottom', label='ERA5 Vorticity (x10$^{-'+order+'}$  s$^{-1}$)',pad=0.02)
    
    scale = 1000000
    ax2 = plt.subplot(222, projection=sps)
    ax2.coastlines()
    ax2.gridlines()
    field = ax2.pcolormesh(lons, lats, ice_vort*scale, transform=pc, cmap=cmap1, vmin = -1, vmax=1)
    ax2.set_extent(extent, crs=pc)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    order   = str(int(np.log10(scale)))
    cbar    = fig.colorbar(field, ax=ax2, location ='bottom', label='Ice Vorticity (x10$^{-'+order+'}$  s$^{-1}$)',pad=0.02)

    scale = 10000
    ax3 = plt.subplot(223, projection=sps)
    ax3.coastlines()
    ax3.gridlines()
    field = ax3.pcolormesh(lons, lats, atmos_div*scale, transform=pc, cmap=cmap1, vmin = -1, vmax=1)
    ax3.set_extent(extent, crs=pc)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    order   = str(int(np.log10(scale)))
    cbar    = fig.colorbar(field, ax=ax3, location ='bottom', label='ERA5 Divergence (x10$^{-'+order+'}$  s$^{-1}$)',pad=0.02)
    
    scale = 1000000
    ax4 = plt.subplot(224, projection=sps)
    ax4.coastlines()
    ax4.gridlines()
    field = ax4.pcolormesh(lons, lats, ice_div*scale, transform=pc, cmap=cmap1, vmin = -1, vmax=1)
    ax4.set_extent(extent, crs=pc)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.03, 0.7])
    order   = str(int(np.log10(scale)))
    cbar    = fig.colorbar(field, ax=ax4, location ='bottom', label='Ice Divergence (x10$^{-'+order+'}$  s$^{-1}$)',pad=0.02)


    plt.subplots_adjust(wspace=0.05, hspace=0.05)
#%%
def mask_era(atmos, ice):
    
    masked_era = (atmos*ice)/ice
    masked_ice = (ice*masked_era)/masked_era
    
    return masked_era, masked_ice
#%%
def compute_corr_row(dates, fields):

    date  = dates
  
    atmos    = fields[0]
    atmos_1d = list(atmos.flatten())
    atmos_1d = [x for x in atmos_1d if np.isnan(x) == False]

    ice    = fields[1]
    ice_1d = list(ice.flatten())
    ice_1d = [x for x in ice_1d if np.isnan(x) == False]

    if len(atmos_1d) == len(ice_1d):

        pearson          = np.corrcoef(atmos_1d, ice_1d, rowvar=True)[0][1]
        np.random.shuffle(ice_1d)
        pearson_shuffled = np.corrcoef(atmos_1d, ice_1d, rowvar=True)[0][1]
        row              = pd.DataFrame([[date, pearson, pearson_shuffled, len(atmos_1d)]], columns = ['date', 'pearson', 'pearson_shuffled', 'pixel_count'])

        return row

    elif len(atmos_1d) != len(ice_1d):
        print("ERROR: Different number of non-NaN pixels.")
        return 0

    else:
        print("UNKNOWN ERROR: Reached else-condition in compute_corr_row(...)")
        return 1
#%%
if __name__ == '__main__':

    import numpy as np
    import xarray as xr
    import metpy.calc as mpcalc
    from metpy.units import units
    import matplotlib as mpl
    import time
    import datetime
    import pyresample as pr
    from multiprocessing import Pool, current_process
    import pandas as pd
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt


    tic = time.time()

    #parameters
    hemisphere      = 'south'
    grid_proj       = 'ease2' #ease2 or polarstereo
    sea             = 'amund-bell'
    sample_radius   = 75
    cores           = 18
    sensor          = 'merged'
    flag_threshold  = 20
    version         = 'v002'
    
    #sectors
    sector = {
        'global'         :           np.arange(-180,    180+0.25, 0.25),
        'weddell'        :           np.arange( -70,    -14+0.25, 0.25),
        'king_hakon'     :           np.arange( -14,     71+0.25, 0.25),
        'east_antarctica':           np.arange(  71,    162+0.25, 0.25),
        'ross'           : np.append(np.arange( 162, 179.75+0.25, 0.25), np.arange(-180, -110.5, 0.25)),
        'amund-bell'     :           np.arange(-110,    -70+0.25, 0.25)
    }

    #directories
    dir_era5         = f'/home/waynedj/Data/era5/{hemisphere}/'
    dir_ice          = f'/home/waynedj/Data/seaicedrift/osisaf/24hr/south/{sensor}/'
    dir_output       = f'/home/waynedj/Projects/seaicedrift_correlation/intermediate_data/24hr/'

    #define generic input grid
    fname          = dir_era5 + 'ERA5_vorticity_1000hp_2020.nc'
    ds             = xr.open_dataset(fname)
    ds_sector      = ds.sel(longitude=sector[sea])
    in_lon         = ds_sector.longitude.values
    in_lat         = ds_sector.latitude.values
    in_lon, in_lat = np.meshgrid(in_lon, in_lat)

    #define generic output grid
    fname            = dir_ice + f'/2020/06/ice_drift_sh_ease2-750_cdr-v1p0_24h-202006301200.nc'
    ds               = xr.open_dataset(fname)
    out_lon          = ds.lon.values
    out_lat          = ds.lat.values

    #define period in units: years
    years = np.arange(1991, 2021)
    
    for year in years:
        print(f"PROCESSING YEAR: {year}")

        #define winter season
        date_0     = datetime.datetime(year, 4, 2,12,0,0) #2013, 4, 2,12,0,0
        date_1     = datetime.datetime(year,11,30,12,0,0) #2013,11,30,12,0,0
        delta_day  = datetime.timedelta(days=1)
        delta_hr   = datetime.timedelta(hours=1)
        date_list  = pd.date_range(date_0, date_1, freq='D')

        #fetch ERA5 Surface Fields, segment appropriate sector
        fname           = dir_era5 + f'ERA5_vorticity_1000hp_{year}.nc'
        ds              = xr.open_dataset(fname)
        ds_sector       = ds.sel(longitude=sector[sea])
        
        #compute 24hr mean field
        atmos_vort_list = []
        ice_vort_list   = []
        time0 = time.time()
        for date in date_list:
            atmos_vort_list.append(compute_mean_atmos_vort(date))
            ice_vort_list.append(compute_mean_ice_vort(date))
        atmos_vort_list = np.array(atmos_vort_list)
        ice_vort_list   = np.array(ice_vort_list)
        time1 = time.time()
        time_ave = time1-time0

        #POOL 1: remap era5_vort onto OSI455 grid
        time0 = time.time()
        p1 = Pool(cores)
        atmos_vort_list_remapped = p1.starmap(swath_remapping, zip(atmos_vort_list))
        p1.close()
        p1.join()
        time1 = time.time()
        time_p1 = time1-time0
        
        #POOL 3: mask out all era5_vort outside of ice field
        time0 = time.time()
        p3 = Pool(cores)
        masked_vort_fields = p3.starmap(mask_era, zip(atmos_vort_list_remapped, ice_vort_list))
        p3.close()
        p3.join()
        time1 = time.time()
        time_p3 = time1-time0
        
        #POOL 5: flatten both era5_vort and ice_vort fields, compute statistics, append as row on df
        time0 = time.time()
        p5 = Pool(cores)
        daily_vort_rows = p5.starmap(compute_corr_row, zip(date_list, masked_vort_fields))
        p5.close()
        p5.join()
        time1 = time.time()
        time_p5 = time1-time0

        #convert to user friendly dataframes, save to csv
        df_daily_corr_vort = pd.concat(daily_vort_rows, ignore_index=True)
        df_daily_corr_vort.to_csv(dir_output+f'vorticity/driftfieldcorrelation_24hr_{sensor}_{sea}_{year}_{version}.csv', index=False)
        
        print('PROCESSING TIME (s): MeanField00: {:.3f} | Pool01: {:.3f} | Pool03: {:.3f} | Pool05: {:.3f} |'.format(time_ave, time_p1, time_p3, time_p5))
        
    toc = time.time()
    print('> Script Completed in {:.3f} seconds'.format(toc-tic))