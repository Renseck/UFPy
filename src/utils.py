# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:20:12 2024

@author: rens_
"""
from colorsys import hsv_to_rgb, rgb_to_hsv
from sklearn.metrics import mean_squared_error, r2_score
from tqdm import tqdm

import numpy as np


def order_of_magnitude(number):
    """Calculates order of magnitude of number"""
    return int(np.floor(np.log10(number)))

def q2RH(q, p, T):
    """Calculates relative humidity from specific humidity"""
    T0 = 273.15
    return 0.263*p*q / (np.exp((17.67*(T - T0)) / (T - 29.65)))

def RH2q(RH, p, T):
    """Calculates specific humidity from relative humidity"""
    T0 = 273.15
    return RH * (np.exp((17.67*(T - T0)) / (T - 29.65))) / (0.263*p)

def complementary_color(r, g, b):
   """returns RGB components of complementary color"""
   hsv = rgb_to_hsv(r, g, b)
   return hsv_to_rgb((hsv[0] + 0.5) % 1, hsv[1], hsv[2])

def print_underlined(text):
    print('\033[4m' + text + '\033[0m')
    
def find_optimum(model_dist, measurement_dist):
    """Calculates time of optimum rmse, r2 and corr between two distributions"""
    corrs = []
    rmses = []
    rsquareds = []
    for row in tqdm(range(len(model_dist)-1), desc = "Calculating optimum time..."):
        rmse = mean_squared_error(model_dist[measurement_dist.columns].iloc[row], measurement_dist.iloc[0])
        r2 = r2_score(measurement_dist.iloc[0], model_dist[measurement_dist.columns].iloc[row])
        corr = np.corrcoef(model_dist[measurement_dist.columns].iloc[row], measurement_dist.iloc[0])[0,1]
        rmses.append(rmse)
        rsquareds.append(r2)
        corrs.append(corr)

    print(f"\nMin RMSE of {np.min(rmses):.2e} at time = {np.argmin(rmses)}s")
    print(f"Max R2 of {np.max(rsquareds):.2f} at time = {np.argmax(rsquareds)}s")
    print(f"Max correlation of {np.max(corrs):.2f} at time = {np.argmax(corrs)}s")
    
def lognormal(x, sigma, center_x, scale = 1):
    """Returns **normalized** lognormal curve"""
    prefactor = 1/(np.log(sigma)*np.sqrt(2*np.pi))
    variable = np.exp((-np.log(x*(1/center_x))**2)/(2*np.log(sigma)**2))
    return scale*prefactor*variable