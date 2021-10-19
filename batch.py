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

yrs_to_eq = 15
sweep_rates = np.arange(50,400,100)
all_vars_dict = {}
rep_lat = 75
ML_depth = [50,100,200,300]
start_co2 = 1
end_co2 = 256
N_min = 1
N_amp = 4

fig = plt.figure(1,figsize = (11,8))

for j,rate in sweep_rates:
    all_vars = np.zeros(4,5,(yrs_to_eq+rate)*365)
    for i,Dml in enumerate(ML_depth):
        print(rate,Dml)
        V,A,Tml,Ti,co2 = run_batch(rep_lat,Dml,start_co2,end_co2,yrs_to_eq,rate,N_min,N_amp)
        all_vars[i,:] = np.array([V,A,Tml,Ti,co2])
    
    all_vars_dict[rate] = all_vars
    

# np.save('/n/home13/chankel/Arctic/Toy_models/Eisenman2007/Output/max_vals_tuning.npy',max_vals)
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