# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 10:20:12 2024

@author: rens_
"""
import numpy as np
from colorsys import rgb_to_hsv, hsv_to_rgb

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