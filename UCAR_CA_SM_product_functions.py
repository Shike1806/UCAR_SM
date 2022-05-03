#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 14:02:04 2022

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
import pickle

def fun_wrapper(inp):
    day = inp[0]
    dd = inp[1]
    cali_fn = inp[2]
    grid_fn = inp[3]
    
    Preff_all = []
    lon_cyg = []
    lat_cyg = []
    
    # print(str(day))
    d = dd+str(day)+"/*.nc"
    f = sorted(glob(d))
    # f = [f[0]]
    for file in f:
        print(file)
        data = read_data(file)
        # for day in days: # iterate through all days
        for ch in range(4):
            print("channel "+str(ch))
            idx_loc = filter_loc(data, ch) # filter location North America
            idx_qf = filter_qf(data['qf'][idx_loc,ch], idx_loc) # filter quality flags
            idx_qc = filter_qc(data, idx_qf, ch) # filter quality control
            
            
            if idx_qc.shape[0] != 0:
                ds = get_data(data, idx_qc, ch)
                Preff = calc_eff_refl(ds)
                
                Preff_all.extend(Preff.tolist())
                lon_cyg.extend(ds[7].tolist())
                lat_cyg.extend(ds[6].tolist())
                
    [SM_all, SM_cell_idx] = get_SM(Preff, lon_cyg, lat_cyg, cali_fn, grid_fn)
    
    return SM_all, SM_cell_idx

def filter_loc(data, ch):
    '''
    Input: data, qf filtered indices, channel
    Output: indices filtered by location
    '''
    
    # north america
    sp_lat = data["sp_lat"][:,ch]
    sp_lon = data["sp_lon"][:,ch]
    sp_lon[sp_lon>180] = sp_lon[sp_lon>180]-360
    lat_max = 40
    lat_min = 25
    lon_min = -100
    lon_max = -75
    
    idx_loc = np.where((sp_lat<lat_max) & (sp_lat>lat_min) & (sp_lon<lon_max) & (sp_lon>lon_min))[0]
    
    lat = sp_lat[idx_loc]
    lon = sp_lon[idx_loc]

    return idx_loc

def read_qf(qf):
    '''
    Input: raw quality flag of CYGNSS data for one channel as an array
    Output: quality flag as 32 bits
    '''

    temp = np.array([qf], dtype=np.uint32)
    temp = np.unpackbits(temp.view(dtype=np.uint8))
    temp = temp.reshape([qf.shape[0],32])
    qf_b = np.concatenate((np.fliplr(temp[:,0:8]), 
                            np.fliplr(temp[:,8:16]), 
                            np.fliplr(temp[:,16:24]), 
                            np.fliplr(temp[:,24:32])), axis=1)
    return qf_b

def read_data(fn):
    '''
    Input: path and filename
    Output: data from file:
        - DDM observation
        - quality flag
        - specular point lat and lon
        - receiver gain
        - transmitter power and gain
        - distance from specular point to receiver and transmitter
    '''
    ds = nc.Dataset(fn)
    data = {
        'DDM_obs' : ds['power_analog'],
        'qf' : ds['quality_flags'],
        "time": ds.variables["ddm_timestamp_utc"],
        'sp_inc' : ds['sp_inc_angle'], # specular point incidence angle
        'SNR' : ds['ddm_snr'], # SNR
        'sp_lat' : ds['sp_lat'],
        'sp_lon' : ds['sp_lon'],
        'Gr' : ds['sp_rx_gain'], # receiver gain in dBi
        'Pt' : ds['gps_tx_power_db_w'], # transmitter power dBW
        'Gt' : ds['gps_ant_gain_db_i'], # transmitter gain in dBi
        'Rt' : ds['tx_to_sp_range'], # transmitter to specular point in m
        'Rr' : ds['rx_to_sp_range'] # receiver to specular point in m
        }
    return data

def filter_qf(qf, idx_loc):
    '''
    Input: quality flags
    Output: indices of passed quality flags
        - Bit 2: ...
    '''
    qf_b = read_qf(qf) # qf of all channels

    idx_dict = {
        "idx_sb" : np.where(qf_b[:,1] == 0), # s-band power up flag Bit 2
        "idx_sa" : np.where(qf_b[:,2] == 0), # spacecraft attitude error flag Bit 3
        "idx_bb" : np.where(qf_b[:,4] == 0), # black body flag Bit 5
        "idx_tp" : np.where(qf_b[:,7] == 0), # DDM test pattern Bit 8
        "idx_lo" : np.where(qf_b[:,10] == 1), # land surface observation flag Bit 11
        "idx_ds" : np.where(qf_b[:,15] == 0), # direct signal Bit 16
        "idx_lc" : np.where(qf_b[:,16] == 0) # low confidence in gps eirp Bit 17
        }
    
    key_list = list(idx_dict.keys())
    temp = idx_dict[key_list[-1]]
    key_list.pop()
    
    
    for key in list(key_list):
        temp = np.intersect1d(temp, idx_dict[key])
    
    idx = idx_loc[temp]
    return idx

def filter_qc(data, idx, ch):
    '''
    Input:
    Output: filtered by quality control 
        - SNR > 2dB
        - Gr < 0dB
        - sp_inc > 65 deg
        - Pr in delay bin outside pixel 7-10
        - SNR <= Gr+14 dB
    '''
    SNR = data["SNR"][:,ch]
    Gr = data["Gr"][:,ch]
    sp_inc = data["sp_inc"][:,ch]
    

    temp = np.where(data["SNR"][:,ch]>2)
    idx = np.intersect1d(temp, idx)
    
    temp = np.where(data["Gr"][:,ch]>0)
    idx = np.intersect1d(temp, idx)
    
    temp = np.where(data["sp_inc"][:,ch]<65)
    idx = np.intersect1d(temp, idx)
    
    temp = np.where(sp_inc<65)
    idx = np.intersect1d(temp, idx)
    
    temp = np.where(data["SNR"][:,ch]<= data["Gr"][:,ch]+14)
    idx = np.intersect1d(temp, idx)
    
    delay_idx = np.empty(idx.shape[0])
    for j in range(idx.shape[0]):
        delay_idx[j] = np.floor(np.argmax(data['DDM_obs'][idx[j],ch])/11)
    
    temp = np.where((delay_idx > 7) & (delay_idx < 10))
    idx = np.intersect1d(idx[temp], idx)

    return idx
    


def ft(x,a,b,c):
    # input in deg
    x = x*np.pi/180
    return a-b**(x**c)

def calc_eff_refl(ds):
    Pr = ds[0]
    Gr = ds[1]
    Pt = ds[2]
    Gt = ds[3]
    Rr = ds[4]
    Rt = ds[5]
    
    sp_inc = ds[6]

    lamb = 0.19 # wavelength in m
    # Gamma eff in dB
    Geff = 10*np.log(Pr)-10*np.log(Pt)-10*np.log(Gt)-10*np.log(Gr)+20*np.log(Rt+Rr)-20*np.log(lamb)+20*np.log(4*np.pi)

    a = 2.005351730621715
    b = 1.1150009347864067
    c = 4.262909068248038
    
    f = ft(sp_inc,a,b,c)
    
    # Pref after incidence angle correction
    Peff = Geff/f

    return Peff

def get_data(data, idx, ch):
    '''
    Input:
    Output:
    '''
    Pr = np.empty(idx.shape[0])
    for i in range(idx.shape[0]):
        Pr[i] = np.amax(data['DDM_obs'][idx[i],ch])
    temp = np.argwhere(~np.isnan(Pr))[:,0] # discard all nan data
    Pr = Pr[temp]
    
    # lat/lon location
    lat = data['sp_lat'][idx[temp],ch]
    lon = data['sp_lon'][idx[temp],ch]
    lon[lon>180] = lon[lon>180]-360
    
    # get parameters
    Gr = data['Gr'][idx[temp],ch]
    Pt = data['Pt'][idx[temp],ch]
    Gt = data['Gt'][idx[temp],ch]
    Rr = data['Rr'][idx[temp],ch]
    Rt = data['Rt'][idx[temp],ch]
    sp_inc = data['sp_inc'][idx[temp],ch]

    return Pr, Gr, Pt, Gt, Rr, Rt, lat, lon, sp_inc

def read_subgrid(fn):
    ds = nc.Dataset(fn)
    lon_grid = np.around(ds['lon'], decimals=3)
    lat_grid = np.around(ds['lat'], decimals=3)
    lon_sub = ds['lon_sub']
    lat_sub = ds['lat_sub']
    return lon_grid, lat_grid, lon_sub, lat_sub

def get_SM(Preff_all, lon_cyg, lat_cyg, cali_fn, grid_fn):
    with open(cali_fn, "rb") as fp:
        [cell_idx, subcell_idx, b_all, Pm_all, SMm_all, num_all] = pickle.load(fp)
    (lon_grid, lat_grid, lon_sub, lat_sub) = read_subgrid(grid_fn) # read grid data
    
    
    SM_all = []
    SM_cell_idx = []
    c1 = 0
    for k, idx in enumerate(cell_idx): # iterate over all cells with calibration
        lon_sub_max = lon_sub[idx].max()
        lat_sub_max = lat_sub[idx].max()
        
        lon_sub_min = lon_sub[idx].min()
        lat_sub_min = lat_sub[idx].min()
        temp_SM = []
        c = 0
        for (lon,lat, Preff) in zip(lon_cyg, lat_cyg, Preff_all): # iterate over all cygnss oversvations
            
            # check if cygnss observation in cell
            if (lat<lat_sub_max) & (lat>lat_sub_min) & (lon<lon_sub_max) & (lon>lon_sub_min):
                c+=1
        print(c)
        if c != 0:
            c1+=1
                # find subcell
                temp_sc_lon_idx = np.where(lon_sub[idx]<lon)[0][-1]
                temp_sc_lat_idx = np.where(lat_sub[idx]<lat)[0][-1]
                
                temp_idx = temp_sc_lon_idx*12+temp_sc_lat_idx
                if temp_idx in subcell_idx[k]:
                    
                    b_idx = subcell_idx[k].index(temp_idx)
                    b = b_all[k][b_idx]
                    Pm = Pm_all[k][b_idx]
                    SMm = SMm_all[k][b_idx]
                    SM = b*(Preff-Pm)+SMm
                    temp_SM.append(SM)

        if len(temp_SM) != 0:
            SM = np.mean(np.array(temp_SM))
            print('found SM: '+str(SM))
            SM_all.append(SM)
            SM_cell_idx.append(idx)
    
    return SM_all, SM_cell_idx

















