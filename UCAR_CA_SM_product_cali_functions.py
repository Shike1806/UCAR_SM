#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 12:16:39 2022

@author: shike
"""

import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
from warnings import filterwarnings
import time
from glob import glob
import os
import netCDF4 as nc
import h5py
import multiprocess as mp
import pickle
from scipy.optimize import curve_fit

def fun_wrapper(inp):
    cell = inp[0]
    cell_idx_all = inp[1]
    subcell_idx_all = inp[2]
    Preff_all = inp[3]
    SM_all = inp[4]
    
    cell_num = inp[5]
    len_all_cell = inp[6]
    
    cell_idx_new = []
    subcell_idx_new = []
    b_all = []
    Pm_all = []
    SMm_all = []
    num_all = []
    
    print('Processing Cell: '+str(cell_num+1)+'/'+str(len_all_cell))
    
    
    temp_idx = np.where(cell == np.array(cell_idx_all))[0] 
    
    if temp_idx.shape[0] > 1: # check how many observations of cell are available
        all_subcell_temp = [] 
        for j in temp_idx: # get all subcell idx from calibration period
            all_subcell_temp.extend(subcell_idx_all[j])
        
        # determine subcells with more than 2 observation
        all_subcell =list(set([x for x in all_subcell_temp if all_subcell_temp.count(x) > 2]))
        all_subcell.sort()
        
        temp_sc = []
        temp_b = []
        
        temp_Pm = []
        temp_SMm = []
        
        temp_num = []
        # iterate over available subcells with observations
        for sc in all_subcell:
            P_temp = []
            SM_temp = []
            
            for j in temp_idx: # iterate over all observations from different days
                temp = np.where(sc == np.array(subcell_idx_all[j]))[0] # get idx of P
                if temp.shape[0] != 0: # append idx if observation is in subcell
                    P_temp.append(Preff_all[j][temp[0]])
                    SM_temp.append(SM_all[j])
                
            if len(P_temp) > 2:
                [b, Pm, SMm] = calibrate(P_temp, SM_temp)
                temp_sc.append(sc)
                temp_b.append(b[0])
                temp_Pm.append(Pm)
                temp_SMm.append(SMm)
                temp_num.append(len(P_temp))
            
        # append all b values and subcell idx from one cell
        if len(temp_b) != 0:
            b_all.append(temp_b)
            subcell_idx_new.append(temp_sc)    
            cell_idx_new.append(cell)
            Pm_all.append(temp_Pm)
            SMm_all.append(temp_SMm)
            num_all.append(temp_num)
            
    
    return cell_idx_new, subcell_idx_new, b_all, Pm_all, SMm_all, num_all

def unwrap(res):
    cell_idx_new = []
    subcell_idx_new = []
    b_all = []
    Pm_all = []
    SMm_all = []
    num_all = []
    # b = []
    # Pm = []
    # SMm = []
    
    for i in range(len(res)):
        if len(res[i][0]) != 0:
            cell_idx_new.extend(res[i][0])
            subcell_idx_new.extend(res[i][1])
            b_all.extend(res[i][2])
            Pm_all.extend(res[i][3])
            SMm_all.extend(res[i][4])
            num_all.extend(res[i][5])
    
    # for i in range(len(b_all)):
    #     b.extend(b_all[i])
    #     Pm.extend(Pm_all[i])
    #     SMm.extend(SMm_all[i])
        
    
    return cell_idx_new, subcell_idx_new, b_all, Pm_all, SMm_all, num_all

def ft(x, b):
    return b*x

def calibrate(P,SM):
    P = np.array(P)
    SM = np.array(SM)
    
    Pm = P.mean()
    SMm = SM.mean()
    
    b, cov = curve_fit(f=ft, xdata=P-Pm, ydata=SM-SMm)
    
    return b, Pm, SMm

def get_loc(grid_fn, cell_idx_new, subcell_idx_new, b_all):
    grid_fn = '/media/shike/Shared/data/brown_ocean/filtered_subgrid_MW.nc'
    grid = nc.Dataset(grid_fn)
    
    lon_grid = grid['lon']
    lat_grid = grid['lat']
    lon_sub_grid = grid['lon_sub_c']
    lat_sub_grid = grid['lat_sub_c']
    lon = []
    lat = []
    
    lon_sub = []
    lat_sub = []
    b = []
    
    key_lat = []
    key_lon = []
    for i in range(12):
        key_lon.extend([i*1]*12) # lon
        key_lat.extend([x for x in range(12)])# lon  
        
    for i, cell in enumerate(cell_idx_new): # iterate over every cell
        lon.append(float(lon_grid[cell]))
        lat.append(float(lat_grid[cell]))
        for j in subcell_idx_new[i]: # iterate over every subcell
            
            lon_sub.append(float(lon_sub_grid[cell,key_lon[j]]))
            lat_sub.append(float(lat_sub_grid[cell,key_lat[j]]))
            b.append(float(lon_sub_grid[cell,key_lon[j]]))
        
    lon = np.array(lon)
    lat = np.array(lat)
    
    return lon, lat, lon_sub, lat_sub

