using OrdinaryDiffEq
using Plots
using ClearSky: diurnalfluxfactor

##

const DAY = 86400
const YEAR = DAY*365
const c = 2e6
const L = 3e8   #ice latent heat of fusion J/m^3
const cml = 4e6 #mixed layer heat capacity J/m^3/K
const k = 2     #ice thermal conductivity W/m^2/K
const a = 320   #LW rad W/m^2
const b = 4.6   #LW rad W/m^2/K
const S0 = 1365 #solar constant W/m^2

const aice = 0.65 #ice albedo 
const ao = 0.2 #ocean albedo
const amp = 0.55 #melt pond albedo
const aatm = 0.45 
const γ = 120 # ocean-ice heat exchange coeff # 120 W/m^2/K
const Hml = 50 #m ML depth
const Fentr = 0.5 #W/m^2 
const Kd = 3.3 #Atmospheric heat transport constant W/m^2/K 
const v0 = 0.10/YEAR #sea ice export /s
const leads = 0.05 #minimum lead fraction
const h0 = 0.5 #min thickness of new ice

##

ramp(x) = x > zero(x) ? x : zero(x);

function LWimbalance(T, Tₛ, Tmidlat, N)
    D = Kd*(Tmidlat - Tₛ)
    2a/(2 + N) - D/2 + b*(T - Tₛ + 2Tₛ/(2 + N))
end

fD̄(t, co2₁, co2₂, days) = 15 + 3*log2(co2₁) + (t/days)*(3*log2(co2₂/co2₁))

fD(t, co2₁, co2₂, days) = fD̄(t, co2₁, co2₂, days) + 7*cos(2π*(t - 20)/365)

fN̄(t, co2₁, co2₂, days) = 2 + 0.2*log2(co2₁) + (t/days)*(2*log2(co2₂/co2₁))

fN(t, co2₁, co2₂, days) = ramp(2.2*cos(2π*(t - 45)/365)) + fN̄(t, co2₁, co2₂, days)

#this doesn't reproduce your Q_day calculation exactly, but it's close
function fSW(t)
    #convert seconds to days, with periodicity
    d = rem(t/DAY, 365.0)
    #convert to an angle for true anomaly (time of year, essentially)
    θ = 2π*(d/365.0)
    #calculate solar insolation
    latitude = 80*(π/180)
    obliquity = 23.5*(π/180) 
    (S0/π)*(1 - aatm)*diurnalfluxfactor(latitude, θ, obliquity)
end

function system!(du, u, param, t)

    #unpack prognostics
    Tml, Ti, V, A = u
    #unpack parameters
    co2₁, co2₂, days = param

    #quantities solely dependent on time
    SW = fSW(t)
    D = fD(t, co2₁, co2₂, days)
    N = fN(t, co2₁, co2₂, days)

    Ts = A*Ti + (1 - A)*Tml
    h = V/A
    Fsw = SW
    LWimb = LWimbalance(Ti, Ts, D, N)
    LWSW_flux = -LWimb + (1 - aice)*Fsw
    Fml = (1 - A)*(-LWimb + (1 - ao)*Fsw) - A*γ*Tml + Fentr #LOL FML

    if (Tml <= 0) & (Fml < 0)
        Fni = -Fml
        dTmldt = 0
    else
        Fni = 0
        dTmldt = Fml/(cml*Hml)
    end

    if (Ti >= 0) & (LWSW_flux > 0)
        dTidt = 0
        dVdt = (A*(LWimb - (1 - amp)*Fsw - γ*Tml) - v0*L*V)/L
    else
        dTidt = (LWSW_flux - k*Ti/h)*2/(c*h)
        dVdt = (A*(-k*Ti/h - γ*Tml) + Fni - v0*L*V)/L
    end

    if dVdt > 0
        R = 0
    else
        R = -dVdt*A/(2*V)
    end

    dAdt = Fni/(L*h0) - R - v0*A

    #set output
    @inbounds begin
        du[1] = dTmldt
        du[2] = dTidt
        du[3] = dVdt
        du[4] = dAdt
    end
end

##

#initial conditions
u₀ = [0.1446, 0.1205, 1.42, 0.819]
#time span to solve
tspan = (0.0, 10*YEAR)
#system parameters
param = (7.0, 7.8, tspan[2]/DAY)
#define the problem
prob = ODEProblem(system!, u₀, tspan, param)
#solve it
sol = solve(prob, RadauIIA3())