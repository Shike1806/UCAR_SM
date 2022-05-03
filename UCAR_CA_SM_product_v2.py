#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:48:20 2022

@author: shike
"""

import time
import netCDF4 as nc
import os
import matplotlib.pyplot as plt
import socket
from warnings import filterwarnings
from UCAR_CA_SM_product_functions import *
import multiprocess as mp

if __name__ == '__main__':
    start_time = time.time()
    machine = socket.gethostname()
    # machine = "andoria"

    if machine == "shike-XPS-15-9560":
        # local machine
        dd = "/media/shike/Shared/data/cygnss/2021/"
        days = [152, 153, 154]
        # days = [152]
        day_num = len(days)
        cali_fn = "cali_data_10days.pkl"
        grid_fn = '/media/shike/Shared/data/brown_ocean/filtered_subgrid_MW.nc'
        cali_fn = "cali_data_10days.pkl"
        num_cpus = 4
        net_fn = 'test.nc'
    elif machine == "andoria.ecn.purdue.edu" or "cygnsslnx1":
        # server
        dd = '/scratch/andoria/a/shi443/brown_ocean/data/cygnss/2021/'
        days = os.listdir(dd)
        days.sort
        # days = [10:]
        day_num = len(days)
        cali_fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/calibration_AprJun2020.pkl'
        grid_fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/filtered_subgrid_MW.nc'
        net_fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/SM_product_Jun18.nc'
        num_cpus = 10
        


    p = mp.Pool(num_cpus)
    res = p.map(fun_wrapper, [(day, dd, cali_fn, grid_fn) for day in days])
    day = 169
    inp = [day, dd, cali_fn, grid_fn]
    
    (lon_grid, lat_grid, lon_sub, lat_sub) = read_subgrid(grid_fn)
    # setup netCDF file
    SM = nc.Dataset(net_fn, "w", format="NETCDF4")
    SM.title = "SM"
    
    # create dimensions
    days_dim = SM.createDimension("days", day_num)
    dim = SM.createDimension("dim", None)

    # create variables
    days_d = SM.createVariable("days", np.float64, ("days"))
    SM_d = SM.createVariable("SM", np.float64, ("days", "dim"))
    lon = SM.createVariable("lon", np.float64, ("days", "dim"))
    lat = SM.createVariable("lat", np.float64, ("days", "dim"))
    
    days_d[:] = days
    
    for i in range(len(res)):
        SM_d[i,:] = res[i][0]
        lon[i,:] = lon_grid[res[i][1]]
        lat[i,:] = lat_grid[res[i][1]]
        
    SM.close()
    
    print("--- %s seconds ---" % (time.time() - start_time))