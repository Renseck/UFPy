# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from PychamAnalysis import PychamSimulation
from SalsaAnalysis import SalsaSimulation

pycham1 = PychamSimulation("Example_Run_Output_2")
salsa = SalsaSimulation()

pycham1_bins = pycham1.get_data("size_bin_bounds")
pycham1_nc = pycham1.get_data("particle_number_concentration_dry").set_axis(pycham1_bins.iloc[-1][1:], axis = 1)

salsa_bins = salsa.get_data("radii")
salsa_nc = salsa.get_data("output")

pm01 = pycham1_nc[pycham1_nc.columns[pycham1_nc.columns<=0.1]]
pm25 = pycham1_nc[pycham1_nc.columns[pycham1_nc.columns<=2.5]]
pm01_total = pm01.sum(axis = 1)
pm25_total = pm25.sum(axis = 1)


plt.figure(figsize = (10,6))
plt.imshow(pycham1_nc.values.T, cmap = "viridis", aspect = "auto", interpolation = "none")
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()


plt.figure(figsize = (10,6))
plt.plot(pm01_total)
plt.plot(pm25_total, "--")
plt.show()