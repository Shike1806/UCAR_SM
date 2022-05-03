#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 23:24:35 2022

@author: shike
"""
import numpy as np
from ease_grid import EASE2_grid
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import time
import netCDF4 as nc
from global_land_mask import globe

def get_ease_grid():
    egrid = EASE2_grid(36000)
    assert egrid.shape == (406, 964)
    lat = egrid.latdim
    lon = egrid.londim
    
    lat = lat[(lat>25) & (lat<38)]
    lon = lon[(lon>-100) & (lon<-75)]
    lon, lat = np.meshgrid(lon,lat)
    mask = globe.is_land(lat, lon)

    lon_masked = lon[mask]
    lat_masked = lat[mask]
    return lon_masked, lat_masked

def create_grid(lon, lat):
    '''
    Input: array of lat and lon
    Ouput: grid subgrid as meshgrid
    '''
    dim1 = lon.shape[0]
    dim2 = 13
    dim3 = 12
    lat_all = np.empty([dim1, dim2])
    lon_all = np.empty([dim1, dim2])
    
    lat_center = np.empty([dim1, dim3])
    lon_center = np.empty([dim1, dim3])

    for i in range(lon.shape[0]):
        [lon_all[i,:], lat_all[i,:], lon_center[i,:], lat_center[i,:]] = create_subgrid(lon[i],lat[i])

    return lon_all, lat_all, lon_center, lat_center

def create_subgrid(lon, lat):
    d = 0.02705 # 3 km subgrid
   
    lon_sub = np.linspace(start=lon-6*d, stop=lon+6*d, num=13)
    lat_sub = np.linspace(start=lat-6*d, stop=lat+6*d, num=13)
    
    lon_center = np.linspace(start=lon-5.5*d, stop=lon+5.5*d, num=12)
    lat_center = np.linspace(start=lat-5.5*d, stop=lat+5.5*d, num=12)

    return lon_sub, lat_sub, lon_center, lat_center



if __name__ == '__main__':
    start_time = time.time()
    [lon_e2, lat_e2] = get_ease_grid()
    [lon_sub, lat_sub, lon_sub_c, lat_sub_c] = create_grid(lon_e2, lat_e2)
    
    
    # save to nc file
    net_fn = '/media/shike/Shared/data/brown_ocean/filtered_subgrid_MW.nc'
    
    # setup netCDF file
    subgrid = nc.Dataset(net_fn, "w", format="NETCDF4")
    subgrid.title = "filtered subgrid for SM estimation"
    
    # create dimensions
    dim1 = subgrid.createDimension("dim1", lat_e2.shape[0])
    dim2 = subgrid.createDimension("dim2", 13)
    dim3 = subgrid.createDimension("dim3", 12)

    # create variables
    lon_filtered = subgrid.createVariable("lon", np.float64, ("dim1"))
    lat_filtered = subgrid.createVariable("lat", np.float64, ("dim1"))
    
    lon_subgrid = subgrid.createVariable("lon_sub", np.float64, ("dim1", "dim2"))
    lat_subgrid = subgrid.createVariable("lat_sub", np.float64, ("dim1", "dim2"))
    
    lon_subgrid_center = subgrid.createVariable("lon_sub_c", np.float64, ("dim1", "dim3"))
    lat_subgrid_center = subgrid.createVariable("lat_sub_c", np.float64, ("dim1", "dim3"))
    
    
    lon_filtered[:] = lon_e2
    lat_filtered[:] = lat_e2
    
    lon_subgrid[:,:] = lon_sub
    lat_subgrid[:,:] = lat_sub
    
    lon_subgrid_center[:,:] = lon_sub_c
    lat_subgrid_center[:,:] = lat_sub_c
    subgrid.close()
    

    
    fig = plt.figure(figsize=(10, 8), dpi=80)
    m = Basemap(llcrnrlon=-100, llcrnrlat=25, urcrnrlon=-75,urcrnrlat=38, lat_0 = 40., lon_0 = -80) # show US
    # m = Basemap()
    m.drawcountries(color='black', linewidth=1) # countries
    m.drawstates(color='lightgrey', linewidth=1) # states
    m.drawcoastlines() # coastline
    m.scatter(x=lon_e2, y=lat_e2, s=0.4, c='b')

    
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    