#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 29, 2022
Modified on Tue Mar 29, 2022 

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
from UCAR_CA_SM_product_qc_functions import *
import netCDF4 as nc
import time
import socket
import multiprocess as mp

if __name__ == '__main__':
    filterwarnings("ignore")
    start_time = time.time()
    machine = socket.gethostname()
    # machine = "andoria"

    if machine == "shike-XPS-15-9560":
        # local machine
        dd = "/media/shike/Shared/data/cygnss/2021/"
        days = [152, 153, 154]
        net_fn = '/media/shike/Shared/data/brown_ocean/qf_data_3day.nc'
        num_cpus = 4
    elif machine == "andoria.ecn.purdue.edu" or "cygnsslnx1":
        # server
        dd = '/scratch/andoria/a/shi443/brown_ocean/data/cygnss/2020/'
        days = os.listdir(dd)
        days = [int(x) for x in days]
        days.sort()
        num_cpus = 10
        net_fn =  '/scratch/andoria/a/shi443/brown_ocean/UCAR_CA_SM/qf_data_AprJun2020.nc'

    # os.remove(net_fn)
    
    day_num= len(days)
    
    p = mp.Pool(num_cpus)
    res = p.map(fun_wrapper, [(day, dd) for day in days])
    # inp = [day, dd]
    
    # setup netCDF file
    data_qf = nc.Dataset(net_fn, "w", format="NETCDF4")
    data_qf.title = "filtered land surface data for sm retrieval"
    
    # create dimensions
    days_dim = data_qf.createDimension("days", day_num)
    dim = data_qf.createDimension("dim", None)
    
    # create variables
    days = data_qf.createVariable('days', np.float64, ("dim"))
    
    Pr_all = data_qf.createVariable("Pr", np.float64, ("days", "dim"))
    Gr_all = data_qf.createVariable("Gr", np.float64, ("days", "dim"))
    Pt_all = data_qf.createVariable("Pt", np.float64, ("days", "dim"))
    Gt_all = data_qf.createVariable("Gt", np.float64, ("days", "dim"))
    Rr_all = data_qf.createVariable("Rr", np.float64, ("days", "dim"))
    Rt_all = data_qf.createVariable("Rt", np.float64, ("days", "dim"))
    
    lat_all = data_qf.createVariable("lat", np.float64, ("days", "dim"))
    lon_all = data_qf.createVariable("lon", np.float64, ("days", "dim")) 
    
    sp_inc_all = data_qf.createVariable("sp_inc", np.float64, ("days", "dim"))     
    
    save2nc(Pr_all, Gr_all, Pt_all, Gt_all, Rr_all, Rt_all, lat_all, lon_all,sp_inc_all, days, res)
        
    data_qf.close()
    print("--- %s seconds ---" %(time.time()-start_time))
