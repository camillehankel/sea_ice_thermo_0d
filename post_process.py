#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 20:02:10 2021

@author: camillehankel
"""

import numpy as np
from math import *
import matplotlib.pyplot as plt


def get_monthly_avg(input_arr):
    months = [10,41,72,102,133,163,194,225,253,284,314,345,375]
    total_years = int(len(input_arr)/365)
    monthly_arr = np.zeros(total_years*12,)

    for i in range(total_years*12-1):
        month = i%12
        year = floor(i/12)
        start_day = months[month]+int(year*365)
        end_day = months[month+1]+int(year*365)
        month_avg = np.mean(input_arr[start_day:end_day])
        monthly_arr[i] = month_avg
    
    
    return monthly_arr


def plot_hysteresis(forward,back,variable,month):
    variable_names = ['Sea Ice Effective Thickness [m]','Sea Ice Fraction','Mixed Layer Temperature [C]','Sea Ice Temperature [C]']
    abbr = ['V','A','Tml','Ti']
    months = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
    co2_forward = forward[4,:]
    co2_back = back[4,:]
    plt.figure(1,figsize=(11,8))
    plt.rcParams.update({'font.size':16})
    plt.xlabel('CO2 Parameter')
    plt.ylabel(variable_names[variable])
    start = (month-7)%12
    plt.plot(co2_forward[start::12],forward[variable,start::12],label='Increasing CO2')
    plt.plot(co2_back[start::12],back[variable, start::12],label='Decreasing CO2')
    plt.legend()
    plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/hysteresis_'+abbr[variable]+'_'+months[month-1]+'.pdf')
    plt.clf()
    
    

def process_ML_exp(var,month):
    sweep_rates = np.arange(50,400,100)
    ML_depth = [50,100,200,300]
    ML_dict = np.load('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/saved_arrays/ML_exp.npy',allow_pickle=True)
    real_ML_dict = ML_dict.item()
    start = (month-7)%12
    
    variable_names = ['Sea Ice Effective Thickness [m]','Sea Ice Fraction','Mixed Layer Temperature [C]','Sea Ice Temperature [C]']
    abbr = ['V','A','Tml','Ti']
    months = ['Jan','Feb','Mar','Apr','May','June','July','Aug','Sept','Oct','Nov','Dec']
    for j,depth in enumerate(ML_depth):
        plt.figure(1,figsize=(11,8))
        plt.rcParams.update({'font.size':16})
        for rate in sweep_rates:
            monthly_series = get_monthly_avg(real_ML_dict[rate][j,var,:])
            monthly_co2 = get_monthly_avg(real_ML_dict[rate][j,-1,:])
            plt.plot(monthly_co2[start::12],monthly_series[start::12],label=str(rate)+'yr ramp')
        plt.legend()
        plt.xlabel('CO2 Parameter')
        plt.ylabel(variable_names[var])
        plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/ML_exp/'+str(depth)+'_depth_'+abbr[var]+'_'+months[month-1]+'.pdf')    
        plt.clf()