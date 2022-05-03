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
from datetime import datetime


def fun_wrapper(inp):
    day = int(inp[0])
    days = inp[1]
    fn_cygnss = inp[2]
    fn_smap = inp[3]
    grid_fn = inp[4]
    
    Preff_day = []
    SM_day = []
    cell_idx_day = []
    subcell_idx_day = []

    sm_days = check_smap(fn_smap)

    if days[day] in sm_days: # check if smap file available
        # print(days[day])
        (lon_grid, lat_grid, lon_sub, lat_sub) = read_subgrid(grid_fn) # read grid data
        ds_cygnss = nc.Dataset(fn_cygnss) # read cygnss qc
        
        day_idx = days.index(days[day])
        Preff = calc_eff_refl(ds_cygnss, day_idx)
        if Preff.shape[0] != 0:
            # read smap data (daily)
            sm_day = sm_days.index(days[day])
            
            
            ds_sm = filter_loc_smap(fn_smap[sm_day]) # filter North America location
            
            if ds_sm['flag']:
            
                lon_sm = np.around(ds_sm['lon_am'], decimals=3)
                lat_sm = np.around(ds_sm['lat_am'], decimals=3)
                
                for i in range(lon_sm.shape[0]): # iterate through every cell
                # for i in range(10): # iterate through every cell
                    print('Day ' +str(days[day])+', i: '+str(i+1)+'/'+str(lon_sm.shape[0]))
                    # check which sm location on ease2 grid
                    idx = np.where((lon_sm[i] == lon_grid) & (lat_sm[i] == lat_grid))[0] # index for ease2 grid
                    if idx.shape[0] != 0: # check if available
                        temp_P, temp_subcell_idx = filter_sub_loc(ds_cygnss, lon_sub[idx[0],:], lat_sub[idx[0],:], day, Preff)
                        
                        if len(temp_P) != 0:
                            # append to new overall
                            Preff_day.append(temp_P)
                            SM_day.append(ds_sm['sm_am'][i])
                            cell_idx_day.append(idx[0])
                            subcell_idx_day.append(temp_subcell_idx)
                            # print(idx)
    return Preff_day, SM_day, cell_idx_day, subcell_idx_day

def group_results(res):
    Preff_all = []
    SM_all = []
    cell_idx_all = []
    subcell_idx_all = []
    for i in range(len(res)):
        Preff_all.extend(res[i][0])
        SM_all.extend(res[i][1])
        cell_idx_all.extend(res[i][2])
        subcell_idx_all.extend(res[i][3])

    return Preff_all, SM_all, cell_idx_all, subcell_idx_all

def check_smap(fn_smap):
    days = []
    for f in fn_smap:
        metadata = f.split('/')[-1]
        date = metadata[13:21]
        date = date[0:4]+'/'+date[4:6]+'/'+date[6:8]
        days.append(int(datetime.strptime(date, '%Y/%m/%d').strftime('%j')))
    return days

def read_data(fn):
    '''
    Input: path and filename
    Output: data from file
        - sm am and pm
        - lon and lat
        - date        
    '''
    with h5py.File(fn, mode='r') as f:
        k = list(f.keys()) # variables
        # print(list(f[k[2]]))
        sm_am = f[k[1]]['soil_moisture'][:,:]
        lat_am = f[k[1]]['latitude'][:,:]
        lon_am = f[k[1]]['longitude'][:,:]
        
        sm_pm = f[k[2]]['soil_moisture_pm'][:,:]
        lat_pm = f[k[2]]['latitude_pm'][:,:]
        lon_pm = f[k[2]]['longitude_pm'][:,:]

        
    lat_am[lat_am == -9999] = np.nan
    lon_am[lon_am == -9999] = np.nan
    
    lat_pm[lat_pm == -9999] = np.nan
    lon_pm[lon_pm == -9999] = np.nan

    sm_am[sm_am == -9999] = np.nan
    sm_pm[sm_pm == -9999] = np.nan
    
    
    metadata = fn.split('/')[-1]
    date = metadata[13:21]
    date = date[0:4]+'/'+date[4:6]+'/'+date[6:8]
    ds = {
        'sm_am' : sm_am,
        'sm_pm' : sm_pm,
        'lon_am' : lon_am,
        'lat_am' : lat_am,
        'lon_pm' : lon_pm,
        'lat_pm' : lat_pm,
        'date' : date
        }
    return ds

def filter_loc_smap(fn_smap):
    '''
    Input: data, qf filtered indices, channel
    Output: indices filtered by location
    '''
    
    
    ds = read_data(fn_smap)
    lat_am = ds["lat_am"]
    lon_am = ds["lon_am"]
    
    lat_pm = ds["lat_pm"]
    lon_pm = ds["lon_pm"]
    
    lat_max = 38
    lat_min = 25
    lon_min = -100
    lon_max = -75
    
    # filter location of am data
    idx_loc_am = np.where((lat_am<lat_max) & (lat_am>lat_min) & (lon_am<lon_max) & (lon_am>lon_min))
    if idx_loc_am[0].shape[0] != 0:
    
        idx_lat_am = [np.min(idx_loc_am[0]),np.max(idx_loc_am[0])]
        idx_lon_am = [np.min(idx_loc_am[1]),np.max(idx_loc_am[1])]
        
        lat_am = ds["lat_am"][idx_lat_am[0]:idx_lat_am[1],idx_lon_am[0]:idx_lon_am[1]]
        lon_am = ds["lon_am"][idx_lat_am[0]:idx_lat_am[1],idx_lon_am[0]:idx_lon_am[1]]
        
        # # filter location of pm data
        # idx_loc_pm = np.where((lat_pm<lat_max) & (lat_pm>lat_min) & (lon_pm<lon_max) & (lon_pm>lon_min))
        # idx_lat_pm = [np.min(idx_loc_pm[0]),np.max(idx_loc_pm[0])]
        # idx_lon_pm = [np.min(idx_loc_pm[1]),np.max(idx_loc_pm[1])]
        
        # lat_pm = ds["lat_pm"][idx_lat_pm[0]:idx_lat_pm[1],idx_lon_pm[0]:idx_lon_pm[1]]
        # lon_pm = ds["lon_pm"][idx_lat_pm[0]:idx_lat_pm[1],idx_lon_pm[0]:idx_lon_pm[1]]
        
        sm_am = ds["sm_am"][idx_lat_am[0]:idx_lat_am[1],idx_lon_am[0]:idx_lon_am[1]]
        # sm_pm = ds["sm_pm"][idx_lat_pm[0]:idx_lat_pm[1],idx_lon_pm[0]:idx_lon_pm[1]]
        
        lat_am_all = []
        lon_am_all = []
        
        # lat_pm_all = []
        # lon_pm_all = []
        
        sm_am_all = []
        # sm_pm_all = []
        
        # convert from matrix to array
        for i in range(sm_am.shape[0]):
            idx_am = np.argwhere(sm_am[i,:]==sm_am[i,:])[:,0]
            sm_am_all.extend(sm_am[i,idx_am].tolist())
            lat_am_all.extend(lat_am[i,idx_am].tolist())
            lon_am_all.extend(lon_am[i,idx_am].tolist())
        
        # for i in range(sm_pm.shape[0]):
        #     idx_pm = np.argwhere(sm_pm[i,:]==sm_pm[i,:])[:,0]
        #     sm_pm_all.extend(sm_pm[i,idx_pm].tolist())
        #     lat_pm_all.extend(lat_pm[i,idx_pm].tolist())
        #     lon_pm_all.extend(lon_pm[i,idx_pm].tolist())
    
        
        ds_sm = {
            'sm_am' : np.array(sm_am_all),
            # 'sm_pm' : np.array(sm_pm_all),
            'lat_am' : np.around(np.array(lat_am_all), decimals=3),
            'lon_am' : np.around(np.array(lon_am_all), decimals=3),
            # 'lat_pm' : np.around(np.array(lat_pm_all), decimals=3),
            # 'lon_pm' : np.around(np.array(lon_pm_all), decimals=3),
            'date' : ds['date'],
            'flag': True
            }
    else:
        ds_sm={'flag':False}

    return ds_sm

def get_data(ds, var, i):
    '''
    Input: var is in string for desired variable
    '''
    temp = ds[var][:,:]
    temp = temp[i,:]
    var_out = temp[temp.mask == False]
    return var_out

def calc_eff_refl(ds, i):
    Pr = get_data(ds, 'Pr', i)
    Gr = get_data(ds, 'Gr', i)
    Pt = get_data(ds, 'Pt', i)
    Gt = get_data(ds, 'Gt', i)
    Rr = get_data(ds, 'Rr', i)
    Rt = get_data(ds, 'Rt', i)

    lamb = 0.19 # wavelength in m
    # Gamma eff in dB
    Geff = 10*np.log(Pr)-10*np.log(Pt)-10*np.log(Gt)-10*np.log(Gr)+20*np.log(Rt+Rr)-20*np.log(lamb)+20*np.log(4*np.pi)

    # Pref after incidence angle correction
    Peff = Geff
    return Peff

def read_subgrid(fn):
    ds = nc.Dataset(fn)
    lon_grid = np.around(ds['lon'], decimals=3)
    lat_grid = np.around(ds['lat'], decimals=3)
    lon_sub = ds['lon_sub']
    lat_sub = ds['lat_sub']
    return lon_grid, lat_grid, lon_sub, lat_sub

# def filter_loc(ds_cyg, ds_sm, lon_sub, lat_sub, Peff, day):
def filter_sub_loc(ds_cyg, lon_sub, lat_sub, day, Preff):
    '''
    Input:
        - smap location and soil moisture
        - cygnss location and effective reflectivity
        - subgrid location
    Output: observations idx in subgrid in list form. List has length of 144
    '''
    # # temp
    # lon_sub = lon_sub[idx[0],:]
    # lat_sub = lat_sub[idx[0],:]
    
    # get cygnss location
    temp_lon = ds_cyg['lon'][day,:]
    temp_lat = ds_cyg['lat'][day,:]
    
    lon_cyg = temp_lon[temp_lon.mask == False]
    lat_cyg = temp_lat[temp_lat.mask == False]
    
    # # get smap location
    # lon_sm = ds_sm['lon_am'][temp_sm]
    # lat_sm = ds_sm['lat_am'][temp_sm]
   
    # check grid cell observations
    lon_sub_max = lon_sub.max()
    lat_sub_max = lat_sub.max()
    
    lon_sub_min = lon_sub.min()
    lat_sub_min = lat_sub.min()

    temp_cyg = np.where((lat_cyg<lat_sub_max) & (lat_cyg>lat_sub_min) & (lon_cyg<lon_sub_max) & (lon_cyg>lon_sub_min))[0]
    # temp_cyg = np.where((lat_sm<lat_sub_max) & (lat_sm>lat_sub_min) & (lon_sm<lon_sub_max) & (lon_sm>lon_sub_min))[0]

    P = []
    subcell_idx = []
    counter = 0
    # check subgrid cell for cygnss observations
    
    if temp_cyg.shape[0] != 0:
        # print('found '+str(temp_cyg.shape[0]))
        for i in range(lon_sub.shape[0]-1): # iterate through every subcell
            lon_max = lon_sub[i+1]
            lon_min = lon_sub[i]
            for j in range(lat_sub.shape[0]-1):
                lat_max = lat_sub[j+1]
                lat_min = lat_sub[j]
        
                # find all observations in 
                idx_loc = np.where((lat_cyg[temp_cyg]<lat_max) & (lat_cyg[temp_cyg]>lat_min) & (lon_cyg[temp_cyg]<lon_max) & (lon_cyg[temp_cyg]>lon_min))[0]
                if idx_loc.shape[0] != 0:
                    P.append(np.mean(Preff[temp_cyg[idx_loc]]))
                    subcell_idx.append(counter)
                counter +=1
    return P, subcell_idx

    
    
    
    