#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:59:13 2021

@author: camillehankel
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import RK45
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


Fsw = 300 # SW seasonal cycle amplitude W/m^2 
T_midlat = 20#K
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
N = 4 #optical thickness of atm


def LW_imbalance(T,Ts):
   # Ts = A*Ti + (1-A)*Tml
    D = Kd*(T_midlat - Ts)
    lw_Sfc = 2*a/(2+N) - D/2 +b*(T-Ts+2*Ts/(2+N))
    return lw_Sfc
    

def ode_system(t,SW,state_vars):
    [Tml,Ti,V,A] = state_vars
    Ts = A*Ti + (1-A)*Tml
    h = V/A
    Fsw = SW[floor(t/DAY)%365]
    LWSW_flux = -LW_imbalance(Ti,Ts)+(1-a_ice)*Fsw
    Fml = (1-A)*(-LW_imbalance(Tml,Ts)+(1-a_o)*Fsw) -A*gamma*Tml + Fentr
    
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
        dV_dt = (A*(LW_imbalance(Ti,Ts) - (1-a_mp)*Fsw-gamma*Tml) -v0*L*V)/L
        
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

def solve_ODE_system(saves_per_day,num_days,SW,TML,TI,V0,A0):
    tspan=(0,DAY*num_days)
    dt = DAY/saves_per_day
    y_init = np.array([TML,TI,V0,A0])
    ts = np.arange(0,DAY*num_days+dt,dt)
    sol = solve_ivp(fun=lambda t,y: ode_system(t,SW,y) \
                ,vectorized=False,y0=y_init,t_span=tspan,t_eval=ts,rtol=1e-6,method='LSODA')
    
    x = sol.y
    
    return x


rep_lat = 80

LAT = rep_lat*2*np.pi/360
day = np.arange(91.25,456.25,1)
declination = 23.45*2*np.pi/360*np.sin(day*2*np.pi/365)
h_0 = np.arccos(-np.tan(LAT)*np.tan(declination))
Q_day = (1-a_atm)*S_0/np.pi*(h_0*np.sin(LAT)*np.sin(declination)+np.cos(LAT)*np.cos(declination)*np.sin(h_0))
Q_day[np.where(np.isnan(Q_day)&(declination<0))] = 0
Q_day[np.where(np.isnan(Q_day)&(declination>0))] = (1-a_atm)*S_0*np.sin(LAT)*np.sin(declination[np.where(np.isnan(Q_day)&(declination>0))])

total_yrs = 10
total_days = total_yrs*365
saves_per_day = 2
total_saves = int(total_days*saves_per_day)

Ti_arr = np.zeros(total_saves,)
Tml_arr = 0*Ti_arr
V_arr = 0*Ti_arr
A_arr = 0*Ti_arr

#initial conditions
A = 1
V = 10
T_ml = 18
T_ice = -10


t1 = time.perf_counter()
[Tml_arr,Ti_arr,V_arr,A_arr] = solve_ODE_system(saves_per_day,total_days,Q_day,T_ml,T_ice,V,A)
t2 = time.perf_counter()

plt.plot(Tml_arr)
plt.ylim(0,40)



        