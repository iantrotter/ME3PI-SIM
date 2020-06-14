#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Macroeconomic-Epidemiological Emergency Policy Intervention Simulator (ME3PI-SIM)
=================================================================================

Pandemic in a model of economic growth, with an emergency policy parameter that
simulateneously reduces output and infection rate.

@author: Ian M. Trotter
"""

#########################################
# Importing packages and setting them up
#########################################
import sys
import math
import json
import numba
from numba import prange
import datetime
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from SimulateSIR import SimulateSIR
import scipy.optimize
import scipy.interpolate
from SimulateSIR import SimulateSIR
sns.set()
import argparse


#############################
# Parameters and Definitions
#############################

parser = argparse.ArgumentParser(description='Simulate a macroeconomic growth model under the impact of a pandemic.')

# Name of the configuration file
parser.add_argument('config_file', help='Name of the configuration file.', default='conf_backtest.json')
parser.add_argument('--p', help='Policy impact on GDP in percent.', type=float)
parser.add_argument('--sd', help='Intervention start date. (Format: YYYY-MM-DD)')
parser.add_argument('--ed', help='Intervention end date. (Format: YYYY-MM-DD)')
parser.add_argument('--out', help='Output file. (Pickle of a Pandas DataFrame.)')
args = parser.parse_args()

print "CONFIGURATION FILE:", args.config_file

# Open configuration file
with open(args.config_file) as json_data_file:
    config = json.load(json_data_file)

# Read the population growth parameters
a1, a2 = config['a1'], config['a2']

# Read the parameters for the economic model
delta = config['capital_depreciation_rate'] # Capital depreciation rate
alpha = config['output_elasticity_of_capital'] # Output elasticity of capital
g = config['productivity_growth_rate'] # Total factor productivity growth rate
discount_rate = config['discount_rate'] # Discount rate
beta = 1./(1+discount_rate) # Calculate the discount factor

# Read the parameters for the epidemiological model
u = config['admission_cost'] # Unit cost of hospital admission
h = config['hospitalisation_rate'] # Hospitalisation rate (of infected)
r = config['recovery_rate']  # Recovery rate
k1, k2 = config['k1'], config['k2'] # Mortality rate parameters
q1, q2 = config['q1'], config['q2'] # Infection rate parameters


# Read the scenario-specific values
N0 = config['N0'] # Population
I0 = config['I0'] # Initial number of infected individuals
R0 = config['R0'] # Initial number of recovered individuals
D0 = config['D0'] # Initial number of deceased individuals
b0 = config['b0'] # Initial/Uncontrolled infection rate
A0 = config['A0'] # Total Factor Productivity
K0 = config['K0'] # Capital stock

# These three arguments can be taken either from the command line or from the config file, but command
# line has precedence
if args.sd is None:
    intervention_start_date = datetime.datetime.strptime(config['intervention_start_date'], "%Y-%m-%d") # Intervention date
else:
    intervention_start_date = datetime.datetime.strptime(args.sd, "%Y-%m-%d") # Intervention date
if args.ed is None:
    intervention_end_date = datetime.datetime.strptime(config['intervention_end_date'], "%Y-%m-%d") # Intervention date
else:
    intervention_end_date = datetime.datetime.strptime(args.ed, "%Y-%m-%d") # Intervention date
if args.p is None:
    p = config['p'] # Policy parameter (proportion reduction in production)
else:
    p = args.p
if args.out is None:
    simulation_file = config['simulation_file']
else:
    simulation_file = args.out

# Details for discretisation of the capital stock
minK, maxK, resolutionK = config['minK'], config['maxK'], config['resolutionK']

# Starting and ending dates for the simulation
start_date = datetime.datetime.strptime(config['sim_start_date'], "%Y-%m-%d")
end_date = datetime.datetime.strptime(config['sim_end_date'], "%Y-%m-%d")

# Run-time parameters
eps = 1e-8
v_min = np.log(1e-16)
n_periods = config['n_periods'] # Number of periods for analysing

#############################
# Epidemionlogical Simulation
#############################
print "STARTING EPIDEMIOLOGICAL SIMULATION"
dates = pd.date_range(start_date, periods=n_periods)
df_NSIRD = SimulateSIR(dates, a1, a2, b0, r, intervention_start_date, intervention_end_date, p, q1, q2, k1, k2, N0, I0, R0, D0)

df_NSIRD.to_pickle("simulation_epidemiological_data.pkl")
df_NSIRD.to_excel("simulation_epidemiological_data.xls")

print "COMPLETED EPIDEMIOLOGICAL SIMULATION."


#####################
# Economic Simulation
#####################
print "STARTING ECONOMIC SIMULATION"
@numba.jit(nopython=True, parallel=True)
def createValueFunctions(n_periods, minK, maxK, resolutionK, alpha, beta, g, A0, b0, q1, q2, u, h, N, S, I, R, intervention_start_n, intervention_end_n, p):
    # Backwards induction to find the value functions
    KS = np.linspace(minK, maxK, resolutionK)
    VS = KS*0
    Kstep = (maxK-minK)/(resolutionK-1)
    Vfuncs = np.zeros((resolutionK,n_periods+1))
    
    for i in range(n_periods-1,-1,-1):
        # Calculate total factor productivity for this period
        At = A0*(1+g)**i
                
        # Get the work force for this period
        Nt = S[i]+R[i]

        # Reduce productivity and infection rate for the intervention
        if i >= intervention_start_n and i <= intervention_end_n:
            At = (1.-0.01*p)*At
            b = b0*(1.-0.01*q1*(p**q2))
        else:
            b = b0

        # Get the health costs for this period
        Ht = u*h*b*S[i]*I[i]
        
        # Iterate through the possible levels of capital, and calculate the
        # consumption corresponding to each level
        Vt = KS*0

        # This function does linear interpolation
        def VN(KN):
            kn_ind = (KN-minK)/Kstep
            kn_ind1 = int(math.floor(kn_ind))
            kn_ind2 = int(math.ceil(kn_ind))
            kw = kn_ind - kn_ind1
            return (1-kw)*VS[kn_ind1]+kw*VS[kn_ind2]
        
        for j in prange(resolutionK):
            # Get the capital stock level
            Kt = KS[j]
    
            # Calculate the production
            Yt = At*(Kt**(alpha))*(Nt**(1-alpha))

            minC = 1e6
            #maxC = (1-delta)*Kt + Yt
            maxC = Yt - Ht

            if maxC < 0:
                Vt[j] = v_min
                continue

            # Solve the optimisation problem to find consumption
            def f(Ct):
                KN = (1-delta)*Kt + Yt - Ct - Ht
                if KN < minK: return -1*v_min
                if KN > maxK: KN = maxK
                interpV = VN(KN)
                return -1*(Nt*math.log(Ct/Nt) + beta*interpV)

            # Golden Searh Algorithm
            invphi = (math.sqrt(5) - 1)/2 # 1/phi
            invphi2 = (3-math.sqrt(5))/2 # 1/phi^2
            tol = 1e-8
            (a,b) = (min(minC, maxC), max(minC, maxC))
            h = b - a
            if h <= tol:
                Vt[j] = -1*f(0.5*(a+b))
            else:
                n = int(math.ceil(math.log(tol/h)/math.log(invphi)))
                c = a+invphi2*h
                d = a+invphi*h
                yc = f(c)
                yd = f(d)
                for z in range(n-1):
                    if yc < yd:
                        b = d
                        d = c
                        yd = yc
                        h = invphi*h
                        c = a+invphi2*h
                        yc = f(c)
                    else:
                        a = c
                        c = d
                        yc = yd
                        h = invphi*h
                        d = a+invphi*h
                        yd = f(d)
                if yc < yd: Vt[j] = -1*yc
                else: Vt[j] = -1*yd
            # Done Golden Search Algorithm
            
        Vfuncs[:,i] = Vt
        VS = Vt

    # Reverse the order of the CTG functions
    return Vfuncs

intervention_start_n = np.argmax(pd.date_range(start_date, periods=n_periods)>=intervention_start_date)
intervention_end_n = np.argmax(pd.date_range(start_date, periods=n_periods)>=intervention_end_date)

value_functions = createValueFunctions(n_periods, minK, maxK, resolutionK, alpha, beta, g, A0, b0, q1, q2, u, h, np.array(df_NSIRD['N']), np.array(df_NSIRD['S']), np.array(df_NSIRD['I']), np.array(df_NSIRD['R']), intervention_start_n, intervention_end_n, p)
#sns.heatmap(value_functions)
#plt.show()

KS = np.linspace(minK, maxK, resolutionK)
plt.plot(KS, value_functions[:,0])

# Forward simulation of the specified scenario
data = {'A': [], 'K': [], 'Y': [], 'C': [], 'N': [], 'I': [], 'R': [], 'D': [], 'H': []}
Kt = K0
dates = pd.date_range(start_date, end_date)
Kstep = (maxK-minK)/(resolutionK-1)
for i in range(len(dates)):
    # Get the current date
    date = dates[i]
  
    # Calculate total factor productivity for this period
    At = A0*(1+g)**i

    # Calculate the "work force" for this period: S+R
    Nt = df_NSIRD['S'].iloc[i] + df_NSIRD['R'].iloc[i]

    # Reduce productivity and infection rate for the intervention
    if i >= intervention_start_n and i <= intervention_end_n:
        At = (1.-0.01*p)*At
        b = b0*(1.-0.01*q1*(p**q2))
    else:
        b = b0

    # Get the health costs for this period
    Ht = u*h*b*df_NSIRD['S'].iloc[i]*df_NSIRD['I'].iloc[i]

    # Calculate the production
    Yt = At*(Kt**(alpha))*(Nt**(1-alpha))

    # Get the Cost-to-go
    VS = value_functions[:,i+1]

    # This function does linear interpolation
    def VN(KN):
        kn_ind = (KN-minK)/Kstep
        kn_ind1 = int(math.floor(kn_ind))
        kn_ind2 = int(math.ceil(kn_ind))
        kw = kn_ind - kn_ind1
        return (1-kw)*VS[kn_ind1]+kw*VS[kn_ind2]

    minC = 1e6
    #maxC = (1-delta)*Kt + Yt
    maxC = Yt - Ht

    if maxC < 0:
        print "Spending exceeds production on", dates[i]
        print u, h, b, df_NSIRD['S'].iloc[i], df_NSIRD['I'].iloc[i], Ht, maxC
        sys.exit()
  
    # Solve the optimisation problem to find consumption
    def f(Ct):
        KN = (1-delta)*Kt + Yt - Ct - Ht
        if KN < minK: return -1*v_min
        if KN > maxK: KN = maxK
        interpV = VN(KN)
        return -1*(Nt*math.log(Ct/Nt) + beta*interpV)

    # Golden Searh Algorithm
    invphi = (math.sqrt(5) - 1)/2 # 1/phi
    invphi2 = (3-math.sqrt(5))/2 # 1/phi^2
    tol = 1e-8
    (a,b) = (min(minC, maxC), max(minC, maxC))
    h = b - a
    if h<= tol:
        Ct = 0.5*(a+b)
    else:
        n = int(math.ceil(math.log(tol/h)/math.log(invphi)))
        c = a+invphi2*h
        d = a+invphi*h
        yc = f(c)
        yd = f(d)
        for z in range(n-1):
            if yc < yd:
                b = d
                d = c
                yd = yc
                h = invphi*h
                c = a+invphi2*h
                yc = f(c)
            else:
                a = c
                c = d
                yc = yd
                h = invphi*h
                d = a+invphi*h
                yd = f(d)
        if yc < yd: Ct = c
        else: Ct = d      
    # Done Golden Search Algorithm

    # Save the data
    data['A'].append(At)
    data['K'].append(Kt)
    data['Y'].append(Yt)
    data['C'].append(Ct)
    data['N'].append(Nt)
    data['I'].append(df_NSIRD['I'].iloc[i])
    data['R'].append(df_NSIRD['R'].iloc[i])
    data['D'].append(df_NSIRD['D'].iloc[i])
    data['H'].append(Ht)

    # Update next-period capital for next iteration
    KN = (1-delta)*Kt + Yt - Ct - Ht
    Kt = KN

print "COMPLETED ECONOMICS SIMULATION."

df = pd.DataFrame(data, index=dates)
print "Writing simulation results to file", simulation_file
df.to_pickle(simulation_file)
df.to_excel(config['simulation_file'].replace("pkl", "xls"))
