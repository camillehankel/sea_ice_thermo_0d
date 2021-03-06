#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 21 11:59:13 2021

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
import time

## CONSTANTS
DAY=86400
YEAR = DAY*365
c = 2e6 #ice heat capacity, J/m^3/K
L = 3e8 #ice latent heat of fusion J/m^3
cml = 4e6 #mixed layer heat capacity J/m^3/K
k = 2*2 #ice thermal conductivity W/m^2/K (normally 2)
a = 320 #LW rad W/m^2
b = 4.6 #LW rad W/m^2/K
S_0 = 1365 #solar constant

# T_midlat = 20#K
<<<<<<< Updated upstream
a_ice = .65 #ice albedo 
a_o= .2 #ocean albedo
a_mp = .55 #melt pond albedo
a_atm = .45 
=======

#standard albedo params:
# a_ice = .75 #ice albedo 
# a_o= .1 #ocean albedo
# a_mp = .45 #melt pond albedo
# a_atm = .425 #atmoshperic albedo
                   
#suppressing bifurcation params:
a_ice = .45 #ice albedo 
a_o= .25 #ocean albedo
a_mp = .4 #melt pond albedo
a_atm = .425 #atmoshperic albedo

>>>>>>> Stashed changes
gamma = 120 # ocean-ice heat exchange coeff # 120 W/m^2/K
Hml = 100 #m ML depth
Fentr = .5 #W/m^2 
Kd = 3.3 #Atmospheric heat transport constant W/m^2/K 
v0 = .10/YEAR #sea ice export /s
leads = .05 #minimum lead fraction
h0 = .5 #min thickness of new ice
# N = 4 #optical thickness of atm
<<<<<<< Updated upstream

=======
# min_N = 1.25 #minimum co2 parameter throughout the year, tunable param
# amplitude_N = 4
D_amp = 7
delN = .192
>>>>>>> Stashed changes

def LW_imbalance(T,Ts,T_midlat,N):
   # Ts = A*Ti + (1-A)*Tml
    D = Kd*(T_midlat - Ts)
    lw_Sfc = 2*a/(2+N) - D/2 +b*(T-Ts+2*Ts/(2+N))
    return lw_Sfc
    

def ode_system(t,state_vars,SW,T_mid,N):
    [Tml,Ti,V,A] = state_vars
    # if V<0:
    #     V = 1e-6
    # if A<0:
    #     A = 1e-8
    Ts = A*Ti + (1-A)*Tml
    h = V/A
    Fsw = SW
    LWSW_flux = -LW_imbalance(Ti,Ts,T_mid,N)+(1-a_ice)*Fsw
    Fml = (1-A)*(-LW_imbalance(Tml,Ts,T_mid,N)+(1-a_o)*Fsw) -A*gamma*Tml + Fentr
    
    
    # if Fml < 0:
    #     if T_ml>10:
    #         F_ni = 0
    #         dTml_dt = Fml/(cml*Hml)
            
    #     elif T_ml<10:
    #         F_ni = -Fml
    #         dTml_dt = 0
    #     else:
    #         F_ni = -Fml/(1+math.exp(Tml*10))
    #         dTml_dt = Fml/(cml*Hml)/(1+math.exp(-Tml*10))
        
    # else:
    #     F_ni = 0
    #     dTml_dt = Fml/(cml*Hml)
           
    # no_sfc_melt_dTi = (LWSW_flux - k*Ti/h)*2/(c*h)
    # sfc_melt_dV = (A*(LW_imbalance(Ti,Ts,T_mid,N) - (1-a_mp)*Fsw-gamma*Tml) -v0*L*V)/L
    # no_sfc_melt_dV = (A*(-k*Ti/h-gamma*Tml) + F_ni - v0*L*V)/L
    # if (LWSW_flux < -100) or (Ti < -100):
    #     dTi_dt = no_sfc_melt_dTi
    #     dV_dt = no_sfc_melt_dV
    # else:
    #     dTi_dt = no_sfc_melt_dTi- no_sfc_melt_dTi/(1+(math.exp(-Ti*2)+math.exp(-LWSW_flux/2))/2)
    #     dV_dt = no_sfc_melt_dV - (no_sfc_melt_dV-sfc_melt_dV)/(1+(math.exp(-Ti*2)+math.exp(-LWSW_flux/2))/2)
    
    if (Tml <= 0) & (Fml < 0): #ice forming due to freezing mixed layer
        # print("if statment 1")
        F_ni = -Fml
        dTml_dt = 0
    else:
        # print("not if statement 1")
        F_ni = 0
        dTml_dt = Fml/(cml*Hml)
    
    
    if V < 0:
        dV_dt = (.01-V)/DAY
        dTi_dt = -Ti/DAY
    else:
        if (Ti >= 0) & (LWSW_flux>0): #ice surface melting
            # print("if statement 2")
            dTi_dt = 0
            dV_dt = (A*(LW_imbalance(Ti,Ts,T_mid,N) - (1-a_mp)*Fsw-gamma*Tml) -v0*L*V)/L
        
        else:
            # print('LWSW',LWSW_flux)
            # print('diffusion',k*Ti/h)
            # print("not if statementn 2")
            if V < .01:
                dTi_dt = -Ti/DAY
            else:
                dTi_dt = (LWSW_flux - k*Ti/h)*2/(c*h)
            
            dV_dt = (A*(-k*Ti/h-gamma*Tml) + F_ni - v0*L*V)/L
        
    # print(dV_dt)
    if dV_dt > 0: 
        # print("if statement 3")
        R = 0
    else:
        R = -dV_dt * A/(2*V) 
        # print("not if statement 3")
        
    if A < 0:
        dA_dt = (.001-A)/DAY
    else:
        dA_dt =  F_ni/(L*h0) -R - v0*A
    # print(Fml)
    
    # print('vars',state_vars)
    # print('derivs',[dTml_dt,dTi_dt,dV_dt,dA_dt])
    if np.sum(np.isnan(state_vars)):
        sys.exit()
    
    return [dTml_dt,dTi_dt,dV_dt,dA_dt]

def solve_ODE_system(coupling,saves_per_day,SW,D,N,TML,TI,V0,A0):
    tspan=(0,DAY*coupling)
    dt = DAY/saves_per_day
    y_init = np.array([TML,TI,V0,A0])
    ts = np.arange(0,DAY*coupling+dt,dt)
    
    # r = ode(ode_system).set_integrator('zvode',method='bdf',min_step=60)
    # y0 = y_init
    # t0 = 0
    # r.set_initial_value(y0,t0).set_f_params(SW,D,N)
    # i = 0
    # out = np.zeros((4,len(ts)))
    # while r.successful() and r.t<(DAY*coupling+dt):
    #     out[:,i] = r.integrate(r.t+dt)
    
    sol = solve_ivp(fun=lambda t,y: ode_system(t,y,SW,D,N) \
                ,vectorized=False,y0=y_init,t_span=tspan,t_eval=ts,atol=1.0e-2,rtol=1.0e-2,method='LSODA',min_step=1e-7*3600,max_step=.125*3600)
    # if sol.status != 0:
    #     print("convergence failed :( ")
    #     print(y_init)
    #     sol = solve_ivp(fun=lambda t,y: ode_system(t,y,SW,D,N) \
    #             ,vectorized=False,y0=y_init,t_span=tspan,t_eval=ts,atol=1.0e-6,rtol=1.0e-6,method='LSODA',min_step=.5e-5*3600,max_step=.5e-5*3600)
    x = sol.y
    # print(x)
    return x

def run_simulation(total_days,coupling_timestep,Q_day,N,D,saves_per_day,saves_per_coupling,T_ml,T_ice,V,A):
    total_saves = int(total_days/coupling_timestep*saves_per_coupling)

    Ti_arr = np.zeros(total_saves,)
    Tml_arr = 0*Ti_arr
    V_arr = 0*Ti_arr
    A_arr = 0*Ti_arr
    
    t1 = time.perf_counter()
    
    for day in np.arange(0,total_days,coupling_timestep):
        if day % 365 == 0:
            # print('year: ',day/365)
            t2 = time.perf_counter()
            # print('time elapsed :',t2-t1)
            t1 = t2
        sw = Q_day[floor(day)%365] #set sw insolation for the day 
        d_day = D[floor(day)] #set midlatidude temp for the day
        n_day = N[floor(day)] #set optical depth for the day
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
                
    return Ti_arr,Tml_arr,V_arr,A_arr
    

<<<<<<< Updated upstream


# will pass these params into function
rep_lat = 80
start_co2 = 2
final_co2 = 20 #factor times present co2
yrs_to_achieve = 200

LAT = rep_lat*2*np.pi/360
day = np.arange(91.25,456.25,1)
declination = 23.45*2*np.pi/360*np.sin(day*2*np.pi/365)
h_0 = np.arccos(-np.tan(LAT)*np.tan(declination))
Q_day = (1-a_atm)*S_0/np.pi*(h_0*np.sin(LAT)*np.sin(declination)+np.cos(LAT)*np.cos(declination)*np.sin(h_0))
Q_day[np.where(np.isnan(Q_day)&(declination<0))] = 0
Q_day[np.where(np.isnan(Q_day)&(declination>0))] = (1-a_atm)*S_0*np.sin(LAT)*np.sin(declination[np.where(np.isnan(Q_day)&(declination>0))])

coupling_timestep = 1
saves_per_coupling = 1
saves_per_day = saves_per_coupling/coupling_timestep


## First run the simulation with fixed co2
yrs_to_eq = 20
total_days = yrs_to_eq*365

cycles = np.arange(0,yrs_to_eq*365,1)
N_mean = 2 + .2*math.log(start_co2,2)
D_mean = 15 +3*math.log(start_co2,2)
N =  2.2*np.cos(2*np.pi*(cycles-45)/365) + N_mean
N[N<N_mean] = N_mean
D = 7*np.cos(2*np.pi*(cycles-20)/365) + D_mean

#initial conditions
A = .819
V = 1.42
T_ml = .1446
T_ice = .1205

Ti_arr_eq,Tml_arr_eq,V_arr_eq,A_arr_eq = run_simulation(total_days, coupling_timestep, Q_day, N, D, saves_per_day, saves_per_coupling, T_ml, T_ice, V, A)

if np.absolute(np.mean(V_arr_eq[-365*int(saves_per_day):])-np.mean(V_arr_eq[-2*365*int(saves_per_day):-365*int(saves_per_day)])) >.001:
    print("May not have achieved equilibrium")

#Last timestep used for new initial conditions
A_new = 1*A_arr_eq[-1]
V_new = 1*V_arr_eq[-1]
Tml_new = 1*Tml_arr_eq[-1]
Ti_new = 1*Ti_arr_eq[-1]

# Make an array of transient co2 forcing
total_days = yrs_to_achieve*365
cycles = np.arange(0,yrs_to_achieve*365,1)
starting_N = 2 + .2*math.log(start_co2,2)
ending_N = 2 + .2*math.log(final_co2,2)
N_mean = starting_N + cycles/(yrs_to_achieve*365-1)*(ending_N-starting_N)
N =  2.2*np.cos(2*np.pi*(cycles-45)/365) + N_mean
N[N<N_mean] = N_mean[N<N_mean]

starting_D = 15 + 3*math.log(start_co2,2)
ending_D = 15 + 3*math.log(final_co2,2)
D_mean = starting_D + cycles/(yrs_to_achieve*365-1)*(ending_D - starting_D)
D = 7*np.cos(2*np.pi*(cycles-20)/365) + D_mean #seasonal cycle of midlatitude temp

#run simulation with transient co2
Ti_arr,Tml_arr,V_arr,A_arr = run_simulation(total_days, coupling_timestep, Q_day, N, D, saves_per_day, saves_per_coupling, Tml_new, Ti_new, V_new, A_new)

total_V = np.concatenate([V_arr_eq,V_arr])
total_A = np.concatenate([A_arr_eq,A_arr])
total_Tml = np.concatenate([Tml_arr_eq,Tml_arr])
total_Ti = np.concatenate([Ti_arr_eq,Ti_arr])
days = range(len(total_V))

=======
def run_batch(rep_lat,ML_depth,start_co2,final_co2,yrs_to_eq,yrs_to_achieve,N_min,N_amp,V0,A0,Tml0,Ti0):
    # will pass these params into function
    global Hml 
    Hml = ML_depth
    
    LAT = rep_lat*2*np.pi/360
    day = np.arange(91.25,456.25,1)
    declination = 23.45*2*np.pi/360*np.sin(day*2*np.pi/365)
    h_0 = np.arccos(-np.tan(LAT)*np.tan(declination))
    Q_day = (1-a_atm)*S_0/np.pi*(h_0*np.sin(LAT)*np.sin(declination)+np.cos(LAT)*np.cos(declination)*np.sin(h_0))
    Q_day[np.where(np.isnan(Q_day)&(declination<0))] = 0
    Q_day[np.where(np.isnan(Q_day)&(declination>0))] = (1-a_atm)*S_0*np.sin(LAT)*np.sin(declination[np.where(np.isnan(Q_day)&(declination>0))])
    
    coupling_timestep = 1
    saves_per_coupling = 1
    saves_per_day = saves_per_coupling/coupling_timestep
    
    
    ## First run the simulation with fixed co2
    total_days = yrs_to_eq*365
    
    cycles = np.arange(0,yrs_to_eq*365,1)
    eq_N_mean = N_min + delN*math.log(start_co2,2)
    eq_D_mean = 15 +3*math.log(start_co2,2)
    N =  N_amp*np.cos(2*np.pi*(cycles-45)/365) + eq_N_mean
    N[N<eq_N_mean] =eq_N_mean
    D = D_amp*np.cos(2*np.pi*(cycles-20)/365) + eq_D_mean
    
    #initial conditions
    A = A0 #2.18e-13   #.9489
    V = V0 #8.42e-7 #3.337
    T_ml = Tml0 #6.4273 #.02257
    T_ice = Ti0 #2.1e-20 #.0004339
    
    Ti_arr_eq,Tml_arr_eq,V_arr_eq,A_arr_eq = run_simulation(total_days, coupling_timestep, Q_day, N, D, saves_per_day, saves_per_coupling, T_ml, T_ice, V, A)
    
    #Last timestep used for new initial conditions
    A_new = 1*A_arr_eq[-1]
    V_new = 1*V_arr_eq[-1]
    Tml_new = 1*Tml_arr_eq[-1]
    Ti_new = 1*Ti_arr_eq[-1]
    
    while np.absolute(np.mean(V_arr_eq[-365*int(saves_per_day):])-np.mean(V_arr_eq[-2*365*int(saves_per_day):-365*int(saves_per_day)])) >.001:
        print("May not have achieved equilibrium, run again")
        Ti_arr_eq,Tml_arr_eq,V_arr_eq,A_arr_eq = run_simulation(total_days, coupling_timestep, Q_day, N, D, saves_per_day, saves_per_coupling, Tml_new, Ti_new, V_new, A_new)
        A_new = 1*A_arr_eq[-1]
        V_new = 1*V_arr_eq[-1]
        Tml_new = 1*Tml_arr_eq[-1]
        Ti_new = 1*Ti_arr_eq[-1]
    
    # Make an array of transient co2 forcing
    total_days = yrs_to_achieve*365
    cycles = np.arange(0,yrs_to_achieve*365,1)
    starting_N = N_min +delN*math.log(start_co2,2)
    ending_N = N_min + delN*math.log(final_co2,2)
    N_mean = starting_N + cycles/(yrs_to_achieve*365-1)*(ending_N-starting_N)
    N =  N_amp*np.cos(2*np.pi*(cycles-45)/365) + N_mean
    N[N<N_mean] = N_mean[N<N_mean]
    
    starting_D = 15 + 3*math.log(start_co2,2)
    ending_D = 15 + 3*math.log(final_co2,2)
    D_mean = starting_D + cycles/(yrs_to_achieve*365-1)*(ending_D - starting_D)
    D = D_amp*np.cos(2*np.pi*(cycles-20)/365) + D_mean #seasonal cycle of midlatitude temp
    
    #run simulation with transient co2
    Ti_arr,Tml_arr,V_arr,A_arr = run_simulation(total_days, coupling_timestep, Q_day, N, D, saves_per_day, saves_per_coupling, Tml_new, Ti_new, V_new, A_new)
    
    #run with fixed co2 after co2 changes until eq
    if yrs_to_achieve !=0:   
        A_new = 1*A_arr[-1]
        V_new = 1*V_arr[-1]
        Tml_new = 1*Tml_arr[-1]
        Ti_new = 1*Ti_arr[-1]
    
    total_days = yrs_to_eq*365
    
    cycles = np.arange(0,yrs_to_eq*365,1)
    end_N_mean = N_min + delN*math.log(final_co2,2)
    end_D_mean = 15 +3*math.log(final_co2,2)
    N_end =  N_amp*np.cos(2*np.pi*(cycles-45)/365) + end_N_mean
    N_end[N_end<end_N_mean] =end_N_mean
    D_end = D_amp*np.cos(2*np.pi*(cycles-20)/365) + end_D_mean
    
    Ti_end,Tml_end,V_end,A_end = run_simulation(total_days, coupling_timestep, Q_day, N_end, D_end, saves_per_day, saves_per_coupling, Tml_new, Ti_new, V_new, A_new)
    
    total_V = np.concatenate([V_arr_eq,V_arr,V_end])
    total_A = np.concatenate([A_arr_eq,A_arr,A_end])
    total_Tml = np.concatenate([Tml_arr_eq,Tml_arr,Tml_end])
    total_Ti = np.concatenate([Ti_arr_eq,Ti_arr,Ti_end])
    
    total_co2 = np.concatenate([eq_N_mean*np.ones(V_arr_eq.shape),N_mean,end_N_mean*np.ones(V_arr_eq.shape)])
    
    
    return total_V,total_A,total_Tml,total_Ti,total_co2
    
    
>>>>>>> Stashed changes
# plt.plot(days,total_V);
# plt.xlabel('Days')
# plt.ylabel('Ice Volume [m/m$^2$ area]')
# plt.ylim(0,4)
# plt.xlim(5000,11000)

# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/V_problem_v5.pdf')
# plt.clf()

# plt.plot(days,total_A);
# plt.xlabel('Days')
# plt.ylabel('Ice Fraction')
# plt.xlim(5000,11000)
# plt.ylim(0,1)
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/A_problem_v5.pdf')
# plt.clf()

# plt.plot(days,total_Tml);
# plt.xlabel('Days')
# plt.ylabel('ML Temperature [C]')
# plt.xlim(5000,11000)
# plt.ylim(0,1)
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Tml_problem_v5.pdf')
# plt.clf()

# plt.plot(days,total_Ti);
# plt.xlabel('Days')
# plt.ylabel('Ice Temperature [C]')
# plt.xlim(5000,11000)
# plt.ylim(-30,1)
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/Ti_problem_v5.pdf')
# plt.clf()


# plt.plot(N_mean,V_arr)
# plt.xlabel('Mean CO2 Forcing')
# plt.ylabel('Ice Volume [m/m$^2$ area]')
# plt.savefig('/Users/camillehankel/Dropbox/Research/Output/Eisenmann_2007/VvsCO2_4xCO2_50yrs.pdf')

    
