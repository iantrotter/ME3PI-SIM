#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@author: Ian Michael Trotter
"""
import pandas as pd

def SimulateSIR(dates, a1, a2, b0, r, intervention_start_date, intervention_end_date, p, q1, q2, k1, k2, N0, I0, R0, D0):
  S0 = N0 - I0 - R0 # Initial number of susceptible individuals

  # Initialise the dataset
  data = {'N': [N0], 'S': [S0], 'I': [I0], 'R': [R0], 'D': [D0]}
  Nm, Sm, Im, Rm, Dm = N0, S0, I0, R0, D0

  # Iterate through the dates
  for i in range(len(dates)):
    # Set appropriate infection rate
    if dates[i] >= intervention_start_date and dates[i] <= intervention_end_date:
      b = b0*(1.-0.01*q1*(p**q2))
      # Set the mortality rate
      m = k1*(b**k2)
    else:
      b = b0      

    # Set the mortality rate
    m = k1*(b**k2)

    births = a1*Nm + a2*Nm*Nm - Nm # Net births
    infections = b*Sm*Im # New infections
    recoveries = r*Im # New recoveries
    deaths = m*Im # New deaths due to pandemic

    N = Nm + births - deaths
    S = Sm + births - infections
    I = Im + infections - recoveries - deaths
    R = Rm + recoveries
    D = Dm + deaths

    data['N'].append(N)
    data['S'].append(S)
    data['I'].append(I)
    data['R'].append(R)
    data['D'].append(D)
        
    # Update the variables for the next iteration
    Nm, Sm, Im, Rm, Dm = N, S, I, R, D

  return pd.DataFrame(data)
    
"""
def SimulateSIR(n_periods, a1, a2, b, r, m, N0, I0, R0, D0):
  S0 = N0 - I0 # Initial number of susceptible individuals

  # Initialise before iterating
  data = {'N': [N0], 'S': [S0], 'I': [I0], 'R': [R0], 'D': [D0]}
  Nm, Sm, Im, Rm, Dm = N0, S0, I0, R0, D0
  # Iterate through the dates
  for i in range(n_periods-1):
    births = a1*Nm + a2*Nm*Nm - Nm # Net births
    infections = b*Sm*Im # New infections
    recoveries = r*Im # New recoveries
    deaths = m*Im # New deaths due to pandemic

    N = Nm + births - deaths
    S = Sm + births - infections
    I = Im + infections - recoveries - deaths
    R = Rm + recoveries
    D = Dm + deaths

    data['N'].append(N)
    data['S'].append(S)
    data['I'].append(I)
    data['R'].append(R)
    data['D'].append(D)
        
    # Update the variables for the next iteration
    Nm, Sm, Im, Rm, Dm = N, S, I, R, D
    
  return pd.DataFrame(data)
"""
