# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:20:12 2024

@author: rens_
"""
from colorsys import hsv_to_rgb, rgb_to_hsv
from sklearn.metrics import mean_squared_error, r2_score
from tqdm import tqdm
from json import loads

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
    """print the input argument with a line under it"""
    print('\033[4m' + text + '\033[0m')
    
def lognormal(x, sigma, center_x, scale = 1):
    """Returns **normalized** lognormal curve"""
    prefactor = 1/(np.log(sigma)*np.sqrt(2*np.pi))
    variable = np.exp((-np.log(x*(1/center_x))**2)/(2*np.log(sigma)**2))
    return scale*prefactor*variable

def parse_metadata(metadata):
    """
    Parse metadata string to dictionary for easier variable reading

    Parameters
    ----------
    metadata : STR
        String of metadata as outputted by read_model_metadata.

    Returns
    -------
    parsed : DICT
        Dictionary of metadata.

    """
    parsed = {}
    for line in metadata.split("\n"):
        if "=" in line:
            name, value = line.split("=")
            name = name.strip()
            try:
                parsed[name] = float(value)
                if name not in ["pt", "pqm1", "pap", "dispersion_rate", "particle_flux"]:
                    parsed[name] = int(value)
            except:
                if name == "particle_flux" :
                    # This is a disgusting hack to make a string of a list into an actual list
                    parsed[name] = loads(value)
                else:
                    parsed[name] = value.strip()
        else:
            break
            
    return parsed