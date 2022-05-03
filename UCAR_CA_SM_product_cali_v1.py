#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19, 2022
Modified on Tue Apr 19, 2022 

Author: Kevin Shi, Purdue University

This is a recreation of the official soil moisture product from the CYGNSS
mission based on doi:10.3390/rs12101558
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from warnings import filterwarnings
import time
from glob import glob
import os
from UCAR_CA_SM_product_cali_functions import *
import netCDF4 as nc
import h5py
import multiprocess as mp
import pickle
from scipy.optimize import curve_fit
import socket


if __name__ == '__main__':
    start_time = time.time()
    
    machine = socket.gethostname()
    # machine = "andoria"
     
    if machine == "shike-XPS-15-9560":
        # local machine
        grid_fn = '/media/shike/Shared/data/brown_ocean/filtered_subgrid_MW.nc'
        net_fn = '/media/shike/Shared/data/brown_ocean/calibration.nc'
        dir_smap = "/media/shike/Shared/data/smap/2021/*.h5"
        num_cpus = 4
        fn1 = '/home/shike/Documents/brown_ocean/UCAR_CA_SM/pre_cali_data_10days.pkl'
        fn2 = '/home/shike/Documents/brown_ocean/UCAR_CA_SM/cali_data_10days.pkl'
    elif machine == "andoria.ecn.purdue.edu" or "cygnsslnx1":
        # server
        # get calibration period info
        net_fn =  '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/calibration_AprJun.nc' # data
        grid_fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/filtered_subgrid_MW.nc'
        fn1 = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/pre_cali_data_AprJun2020.pkl'
        fn2 = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/calibration_AprJun2020.pkl'
        num_cpus = 10
    
    
    with open(fn1, "rb") as fp:
        [Preff_all, SM_all, cell_idx_all, subcell_idx_all] = pickle.load(fp)
        
    all_cell = list(set([x for x in cell_idx_all if cell_idx_all.count(x) > 1]))
    all_cell.sort()
    # print('start')
    p = mp.Pool(num_cpus)
    res = p.map(fun_wrapper, [(cell, cell_idx_all, subcell_idx_all, Preff_all, SM_all, i, len(all_cell)) for i, cell in enumerate(all_cell)])
    # inp = [cell, cell_idx_all, subcell_idx_all, Preff_all, SM_all, i, len(all_cell)]

    [cell_idx_new, subcell_idx_new, b_all, Pm_all, SMm_all, num_all] = unwrap(res)
    
    with open(fn2, "wb") as fp:
        pickle.dump([cell_idx_new, subcell_idx_new, b_all, Pm_all, SMm_all, num_all], fp)
    # [lon, lat, lon_sub, lat_sub] = get_loc(grid_fn, cell_idx_new, subcell_idx_new, b_all)
    
    # # setup netCDF file
    # cali = nc.Dataset(net_fn, "w", format="NETCDF4")
    # cali.title = "filtered land surface data for sm retrieval"
    
    # # create dimensions
    # cell_dim = cali.createDimension("cell", len(lon))
    # scell_dim = cali.createDimension("scell", len(lon_sub))
    
    
    # # create variables
    # lon_loc = cali.createVariable('lon', np.float64, ("cell"))
    # lat_loc = cali.createVariable('lat', np.float64, ("cell"))
    
    # lon_sc = cali.createVariable("lon_sc", np.float64, ("scell"))
    # lat_sc = cali.createVariable("lat_sc", np.float64, ("scell"))
    # b = cali.createVariable("b", np.float64, ("scell"))
    
    # lon_loc[:] = lon
    # lat_loc[:] = lat
    
    # lon_sc[:] = lon_sub
    # lat_sc[:] = lat_sub
    # b[:] = b_all
    
    # cali.close()
   
    
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    