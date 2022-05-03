#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 15 13:59:37 2022

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

def fun_wrapper(inp):
    day = inp[0]
    dd = inp[1]
    
    
    # print(str(day))
    d = dd+str(day)+"/*.nc"
    f = sorted(glob(d))
    Pr_day = []
    Gr_day = []
    Pt_day = []
    Gt_day = []
    Rr_day = []
    Rt_day = []
    lat_day = []
    lon_day = []
    sp_inc_day = []
        
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
                [Pr, Gr, Pt, Gt, Rr, Rt, lat, lon, sp_inc] = get_data(data, idx_qc, ch)

                Pr_day.extend(Pr.tolist())
                Gr_day.extend(Gr.tolist())
                Pt_day.extend(Pt.tolist())
                Gt_day.extend(Gt.tolist())
                Rr_day.extend(Rr.tolist())
                Rt_day.extend(Rt.tolist())
                lat_day.extend(lat.tolist())
                lon_day.extend(lon.tolist())
                sp_inc_day.extend(sp_inc.tolist())
    
    return Pr_day, Gr_day, Pt_day, Gt_day, Rr_day, Rt_day, lat_day, lon_day, sp_inc_day, day

def save2nc(Pr_all, Gr_all, Pt_all, Gt_all, Rr_all, Rt_all, lat_all, lon_all,sp_inc_all, days, res):
    
    temp_day = []
    for i in range(len(res)):
        print('Saving:' +str(res[i][9]))
        if len(res[i][0]):
            Pr_all[i,:] = res[i][0]
            Gr_all[i,:] = res[i][1]
            Pt_all[i,:] = res[i][2]
            Gt_all[i,:] = res[i][3]
            Rr_all[i,:] = res[i][4]
            Rt_all[i,:] = res[i][5]
            lat_all[i,:] = res[i][6]
            lon_all[i,:] = res[i][7]
            sp_inc_all[i,:] = res[i][8]
            temp_day.append(res[i][9])
    days[:] =  temp_day
    return

def filter_loc(data, ch):
    '''
    Input: data, qf filtered indices, channel
    Output: indices filtered by location
    '''
    
    # north america
    sp_lat = data["sp_lat"][:,ch]
    sp_lon = data["sp_lon"][:,ch]
    sp_lon[sp_lon>180] = sp_lon[sp_lon>180]-360
    lat_max = 38
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
