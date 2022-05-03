#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 18:32:20 2022

@author: shike
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def calc_refl(th,eps):
    Rvv = (eps*np.cos(th)-np.sqrt(eps-np.sin(th)*np.sin(th)))/(eps*np.cos(th)+np.sqrt(eps-np.sin(th)*np.sin(th)))       
    Rhh = (np.cos(th)-np.sqrt(eps-np.sin(th)*np.sin(th)))/(np.cos(th)+np.sqrt(eps-np.sin(th)*np.sin(th)))
    G = np.abs(1/2*(Rvv-Rhh))*np.abs(1/2*(Rvv-Rhh))
    return G

# def ft(x,a,b):
    # return a-np.exp(x**b)
    
def ft(x,a,b,c):
    return a-b**(x**c)

if __name__ == '__main__':
    
    n = 360
    th = np.linspace(0,90,n)*np.pi/180
    th_deg = th*180/np.pi
    eps_all = [3, 5, 20, 80]
    G_all = np.empty([len(eps_all), n])
    plt.figure(figsize=[10, 6], dpi=80)
    for i, eps in enumerate(eps_all):
        G = calc_refl(th,eps)
        G_nad = calc_refl(0,eps)
        plt.plot(th_deg, G, label='$\epsilon_s = $'+str(eps))
        # plt.plot(th_deg,G/G_nad, label='$\epsilon_s = $'+str(eps))
        G_all[i,:] = G/G_nad
        
    # uncorrected
    plt.legend()
    plt.grid()
    plt.xlabel(r'$Incidence angle \theta_i$', fontsize=18)
    plt.ylabel(r'$Surface reflectivity \Gamma_{rl}(\epsilon_s,\theta_i)$', fontsize=18)
    plt.show()
    
    G_mean = np.mean(G_all, axis=0)
    plt.figure(figsize=[10, 6], dpi=80)
    
    plt.plot(th_deg,G_all, '--')
    plt.legend()
   
    
    params, cov = curve_fit(f=ft, xdata=th, ydata=G_mean)
    a = params[0]
    b = params[1]
    c = params[2]
    G_fit = ft(th,a,b,c)
    
    # with open('inc_correc_params.txt', 'w') as f:
        # f.write(str(a)+','+str(b)+','+str(c))
        
    plt.figure(figsize=[10, 6], dpi=80)
    for i, eps in enumerate(eps_all):
        plt.plot(th_deg, G_all[i,:], label='$\epsilon_s = $'+str(eps))
    plt.plot(th_deg, G_fit,'k--', label=r'$f(\theta_i)$', linewidth=4)
    plt.legend()
    plt.grid()
    plt.xlabel(r'$Incidence angle \theta_i$', fontsize=18)
    plt.ylabel(r'$Surface reflectivity \frac{\Gamma_{rl}(\epsilon_s,\theta_i)}{\Gamma_{rl}(\epsilon_s,0)}$', fontsize=18)
    plt.show()
