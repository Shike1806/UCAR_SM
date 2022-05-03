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
from UCAR_CA_SM_product_pre_functions import *
import netCDF4 as nc
import h5py
import multiprocess as mp
import pickle
import sys
import socket

if __name__ == '__main__':
    start_time = time.time()

    machine = socket.gethostname()
    # machine = "andoria"
     
    if machine == "shike-XPS-15-9560":
        # local machine
        days = [152, 153, 154]
        net_fn = '/media/shike/Shared/data/brown_ocean/qf_data.nc'
        dir_smap = "/media/shike/Shared/data/smap/2021/*.h5"
        grid_fn = '/media/shike/Shared/data/brown_ocean/filtered_subgrid_NA.nc'
        num_cpus = 4
        fn = '/media/shike/Shared/data/brown_ocean/pre_cali_data_Apr2020.pkl'
    elif machine == "andoria.ecn.purdue.edu" or "cygnsslnx1":
        # server
        # get calibration period info
        net_fn =  '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/qf_data_AprJun2020.nc' # data
        dd = '/scratch/andoria/a/shi443/brown_ocean/data/cygnss/2020/'
        ds = nc.Dataset(net_fn)
        days = ds['days'][:]
        days = list(days[days.mask == False])
        days = [int(x) for x in days]
        fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/pre_cali_data_AprJun2020.pkl'
       
        # get smap data
        dir_smap = "/scratch/andoria/a/shi443/brown_ocean/data/smap/2020/*.h5"
        # get grid info
        grid_fn = '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/filtered_subgrid_MW.nc'
        num_cpus = 10
        
    # days_all = days
    # days = days_all[0:30]
    # smap data directory
    # days = days[0:2]
    
    fn_smap = sorted(glob(dir_smap))
    day_num = len(days)
    
    p = mp.Pool(num_cpus)
    res = p.map(fun_wrapper, [(day,days, net_fn, fn_smap, grid_fn) for day in range(day_num)])

    # inp = [day,days, net_fn, fn_smap, grid_fn]
    
    [Preff_all, SM_all, cell_idx_all, subcell_idx_all] = group_results(res)
            
    with open(fn, "wb") as fp:
        pickle.dump([Preff_all, SM_all, cell_idx_all, subcell_idx_all], fp)
    
    print("--- %s seconds ---" % (time.time() - start_time))
    
    