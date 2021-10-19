#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:47:12 2021

@author: camillehankel
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import ode
import matplotlib.pyplot as plt
from math import floor
from scipy import optimize 
import math
import sys
from model import *

yrs_to_eq = 30
sweep_rates = np.arange(100,450,50)
all_vars_dict = {}
rep_lat = 75
ML_depth = 50
start_co2 = 1
end_co2 = 1
N_mins = np.arange(1.4,6,.4)
N_amps = np.arange(1,6,.3)

min_vals = np.zeros((5,len(N_mins),len(N_amps)))
max_vals = np.zeros((5,len(N_mins),len(N_amps)))

for i,N_min in enumerate(N_mins):
    for j,N_amp in enumerate(N_amps):
        print(i,'\t',j)
        # if (N_amp > -2*N_min + 6 + .5) or (N_amp < -2*N_min + 6 - .5):
        #     min_vals[:,i,j] = [np.nan,np.nan,np.nan,np.nan]
        #     max_vals[:,i,j] = [np.nan,np.nan,np.nan,np.nan]
        #     print('SKIP')
        # else:
        V,A,Tml,Ti,co2 = run_batch(rep_lat,ML_depth,start_co2,end_co2,yrs_to_eq,0,N_min,N_amp,5,1,.02,-.011)
        h = V/A
        all_vars = np.array([V,A,h,Tml,Ti])
        min_vals[:,i,j] = np.amin(all_vars[:,-365:],axis=1)
        max_vals[:,i,j] = np.amax(all_vars[:,-365:],axis=1)
            
        
        
np.save('/n/home13/chankel/Arctic/Toy_models/Eisenman2007/Output/min_tuning_alb.npy',min_vals)
np.save('/n/home13/chankel/Arctic/Toy_models/Eisenman2007/Output/max_tuning_alb.npy',max_vals)
# for rate in sweep_rates:
#     V,A,Tml,Ti,co2 = run_batch(rep_lat,ML_depth,start_co2,end_co2,yrs_to_eq,rate)
#     all_vars = np.array([V,A,Tml,Ti,co2])
#     all_vars_dict[rate] = all_vars
    


# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_V = get_monthly_avg(all_vars_dict[rate][0])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[2::12],monthly_V[2::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('September Effective Ice Thickness [m]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/V_September.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_A = get_monthly_avg(all_vars_dict[rate][1])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[2::12],monthly_A[2::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('September Ice Area')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/A_September.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_Tml = get_monthly_avg(all_vars_dict[rate][2])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[2::12],monthly_Tml[2::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('September ML Temperature [C]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Tml_September.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_Ti = get_monthly_avg(all_vars_dict[rate][3])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[2::12],monthly_Ti[2::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('September Ice Temperature [C]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Ti_September.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_V = get_monthly_avg(all_vars_dict[rate][0])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[8::12],monthly_V[8::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('March Effective Ice Thickness [m]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/V_March.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_A = get_monthly_avg(all_vars_dict[rate][1])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[8::12],monthly_A[8::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('March Ice Area')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/A_March.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_Tml = get_monthly_avg(all_vars_dict[rate][2])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[8::12],monthly_Tml[8::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('March ML Temperature [C]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Tml_March.pdf')
# plt.clf()

# plt.figure(1,figsize=(11,8))
# for rate in sweep_rates:
#     monthly_Ti = get_monthly_avg(all_vars_dict[rate][3])
#     monthly_co2 = get_monthly_avg(all_vars_dict[rate][4])
#     plt.plot(monthly_co2[8::12],monthly_Ti[8::12],label=str(rate)+' yrs')
#     plt.legend()
# plt.xlabel('CO2 Parameter')
# plt.ylabel('March Ice Temperature [C]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Ti_March.pdf')
# plt.clf()