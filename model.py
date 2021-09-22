#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:59:13 2021

@author: camillehankel
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from math import floor
from scipy import optimize 
import math
import time

## CONSTANTS
DAY=86400
YEAR = DAY*365
c = 2e6 #ice heat capacity, J/m^3/K
L = 3e8 #ice latent heat of fusion J/m^3
cml = 4e6 #mixed layer heat capacity J/m^3/K
k = 2 #ice thermal conductivity W/m^2/K
a = 320 #LW rad W/m^2
b = 4.6 #LW rad W/m^2/K
S_0 = 1365 #solar constant

# T_midlat = 20#K
a_ice = .65 #ice albedo 
a_o= .2 #ocean albedo
a_mp = .55 #melt pond albedo
a_atm = .45 
gamma = 120 # ocean-ice heat exchange coeff # 120 W/m^2/K
Hml = 50 #m ML depth
Fentr = .5 #W/m^2 
Kd = 3.3 #Atmospheric heat transport constant W/m^2/K 
v0 = .10/YEAR #sea ice export /s
leads = .05 #minimum lead fraction
h0 = .5 #min thickness of new ice
# N = 4 #optical thickness of atm


def LW_imbalance(T,Ts,T_midlat,N):
   # Ts = A*Ti + (1-A)*Tml
    D = Kd*(T_midlat - Ts)
    lw_Sfc = 2*a/(2+N) - D/2 +b*(T-Ts+2*Ts/(2+N))
    return lw_Sfc
    

def ode_system(t,SW,T_mid,N,state_vars):
    [Tml,Ti,V,A] = state_vars
    Ts = A*Ti + (1-A)*Tml
    h = V/A
    Fsw = SW
    LWSW_flux = -LW_imbalance(Ti,Ts,T_mid,N)+(1-a_ice)*Fsw
    Fml = (1-A)*(-LW_imbalance(Tml,Ts,T_mid,N)+(1-a_o)*Fsw) -A*gamma*Tml + Fentr
    
    if (Tml <= 0) & (Fml < 0): #ice forming
        # print("if statement 1")
        F_ni = -Fml
        dTml_dt = 0
    else:
        F_ni = 0
        dTml_dt = Fml/(cml*Hml)
    
    if (Ti >= 0) & (LWSW_flux > 0):
        # print("if state¡ment 2")

        dTi_dt = 0
        dV_dt = (A*(LW_imbalance(Ti,Ts,T_mid,N) - (1-a_mp)*Fsw-gamma*Tml) -v0*L*V)/L
        
    else:
        dTi_dt = (LWSW_flux - k*Ti/h)*2/(c*h)
        dV_dt = (A*(-k*Ti/h-gamma*Tml) + F_ni - v0*L*V)/L
        
    if -dV_dt < 0: 
        # print("if stateme!nt 3")

        R = 0
    else:
        R = A/(2*V) * -dV_dt
        
    dA_dt =  F_ni/(L*h0) -R - v0*A
    # print(t)
    
    # print('vars',state_vars)
    
    # print('derivs',[dTml_dt,dTi_dt,dV_dt,dA_dt])
    
    return [dTml_dt,dTi_dt,dV_dt,dA_dt]

def solve_ODE_system(coupling,saves_per_day,SW,D,N,TML,TI,V0,A0):
    tspan=(0,DAY*coupling)
    dt = DAY/saves_per_day
    y_init = np.array([TML,TI,V0,A0])
    ts = np.arange(0,DAY*coupling+dt,dt)
    sol = solve_ivp(fun=lambda t,y: ode_system(t,SW,D,N,y) \
                ,vectorized=False,y0=y_init,t_span=tspan,t_eval=ts,rtol=1.e-6,method='LSODA')
    
    
    x = sol.y
    # print(x)
    return x


rep_lat = 80

LAT = rep_lat*2*np.pi/360
day = np.arange(91.25,456.25,1)
declination = 23.45*2*np.pi/360*np.sin(day*2*np.pi/365)
h_0 = np.arccos(-np.tan(LAT)*np.tan(declination))
Q_day = (1-a_atm)*S_0/np.pi*(h_0*np.sin(LAT)*np.sin(declination)+np.cos(LAT)*np.cos(declination)*np.sin(h_0))
Q_day[np.where(np.isnan(Q_day)&(declination<0))] = 0
Q_day[np.where(np.isnan(Q_day)&(declination>0))] = (1-a_atm)*S_0*np.sin(LAT)*np.sin(declination[np.where(np.isnan(Q_day)&(declination>0))])


cycle = np.arange(0,365,1)
D = 7*np.cos(2*np.pi*(cycle-20)/365) + 15 #seasonal cycle of midlatitude temp
N = 2.2*np.cos(2*np.pi*(cycle-45)/365) + 2 #seasonal cycle of optical thickness
N[N<2] = 2

total_yrs = 12
total_days = total_yrs*365
coupling_timestep = 1
saves_per_coupling = 2
total_saves = int(total_days/coupling_timestep*saves_per_coupling)
saves_per_day = saves_per_coupling/coupling_timestep

Ti_arr = np.zeros(total_saves,)
Tml_arr = 0*Ti_arr
V_arr = 0*Ti_arr
A_arr = 0*Ti_arr

#initial conditions
A = 1
V = 10
T_ml = 18
T_ice = -10

for day in np.arange(0,total_days,coupling_timestep):
        print(day)
        sw = Q_day[floor(day)%365] #set sw insolation for the day 
        d_day = D[floor(day)%365] #set midlatidude temp for the day
        n_day = N[floor(day)%365] #set optical depth for the day
        # print(sw)
        start_index = floor(day*saves_per_day)
        end_index = floor(start_index+saves_per_coupling)

        [temp_Tml_arr,temp_Ti_arr,temp_V_arr,temp_A_arr] = solve_ODE_system(coupling_timestep,saves_per_day,sw,d_day,n_day,T_ml,T_ice,V,A)
        Tml_arr[start_index:end_index] = temp_Tml_arr[1:]
        Ti_arr[start_index:end_index] = temp_Ti_arr[1:]
        V_arr[start_index:end_index] = temp_V_arr[1:]
        A_arr[start_index:end_index] = temp_A_arr[1:]
        
        A = temp_A_arr[-1]
        V = temp_V_arr[-1]
        T_ml = temp_Tml_arr[-1]
        T_ice = temp_Ti_arr[-1]


# t1 = time.perf_counter()
# [Tml_arr,Ti_arr,V_arr,A_arr] = solve_ODE_system(saves_per_day,total_days,Q_day,T_ml,T_ice,V,A)
# t2 = time.perf_counter()

plt.plot(Tml_arr)



        