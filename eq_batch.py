#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 12:24:06 2021

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
from model_sbf import *

yrs_to_eq = 10
all_vars_dict = {}
rep_lat = 75
ML_depth = 100
start_co2 = 1
end_co2 = 16
N_min = 2.4#1.6 for nbf when k = 2 #2.4 for nbf when k = 2*2 # 2.4 for sbf
N_amp = .4#.8 for nbf when k = 2 #.6 for nbf when k = 2*2 # .4 for sbf

# co2_start = np.arange(1,6,.5)
# co2_middle = np.arange(5,6.1,.1)
# co2_end = np.arange(6.5,16.5,.5)
# co2_vals = np.concatenate([co2_start,co2_middle,co2_end])

# co2_start = np.arange(1,2.5,.5)
# co2_middle = np.arange(2.1,2.5,.1)
# co2_end = np.arange(2.5,16.5,.5)
# co2_vals_back = np.concatenate([co2_start,co2_middle,co2_end])

# co2_vals = np.arange(1,33,.75)
#co2_vals = np.arange(12.5,14.75,.25) #to fill in the area of fast transition
# co2_vals = [13, 13.25,13.5]# ones that didn't quite equilibrate

co2_vals = np.concatenate([np.arange(1,16,1),np.arange(16,21,.25),np.arange(21,25,1)])
co2_vals_back = np.flip(co2_vals)

init = [3.6,.85,.062,.0016] #(IC's for co2 = 1) [.065,.22685,.5156,.00096] #ICs for warmer 
all_vars_forward = np.zeros((len(co2_vals),5,(2*yrs_to_eq*365)))
j = 0
for co2 in co2_vals:
      print(co2)
      V,A,Tml,Ti,co2 = run_batch(rep_lat,100,co2,co2,yrs_to_eq,0,N_min,N_amp,init[0],init[1],init[2],init[3])
      Vars = np.array([V,A,Tml,Ti,co2])
      all_vars_forward [j,:] = Vars
      init = Vars[:,-1]
      print(init)
      j+=1
      
np.save('/n/home13/chankel/Arctic/Toy_models/Eisenman2007/Output/eq_runs_forward_sbf.npy',all_vars_forward)

init = [1e-5,1e-7,5,1e-6]
j = 0
all_vars_back = np.zeros((len(co2_vals),5,(2*yrs_to_eq*365)))

for co2 in co2_vals_back:
     V,A,Tml,Ti,co2 = run_batch(rep_lat,100,co2,co2,yrs_to_eq,0,N_min,N_amp,init[0],init[1],init[2],init[3])
     Vars = np.array([V,A,Tml,Ti,co2])
     all_vars_back[j,:] = Vars
     init = Vars[:,-1]
     j+=1
          
np.save('/n/home13/chankel/Arctic/Toy_models/Eisenman2007/Output/eq_runs_back_sbf.npy',all_vars_back)
