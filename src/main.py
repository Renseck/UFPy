# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import PychamAnalysis as pyan
import SalsaAnalysis as saan

pycham0 = pyan.PychamSimulation("Example_Run_Output_2")
pycham1 = pyan.PychamSimulation("Slow_Chemical_Reactions")
salsa = saan.SalsaSimulation()

pycham0_bins = pycham0.get_data("size_bin_bounds")
pycham0_nc = pycham0.get_data("particle_number_concentration_dry").set_axis(
    pycham0_bins.iloc[-1][1:], axis=1
)

salsa_bins = salsa.get_data("radii")
salsa_nc = salsa.get_data("output")

pycham0.plume_plot()

fig, ax = plt.subplots(figsize  = (10,6))
pycham0.pnc_plot(ax = ax, label = "Base")
pycham1.pnc_plot(ax = ax, label = "Slow")
plt.show()

pm01_bound = 0.1
pm25_bound = 2.5
pm01 = pycham0_nc[pycham0_nc.columns[pycham0_nc.columns < pm01_bound/2]]
pm25 = pycham0_nc[pycham0_nc.columns[pycham0_nc.columns < pm25_bound/2]]
pm01_total = pm01.sum(axis=1)
pm25_total = pm25.sum(axis=1)

plt.figure(figsize=(10, 6))
plt.plot(pm01_total.index/60, pm01_total, label = f"Dp < {pm01_bound} $\mu m$")
plt.plot(pm25_total.index/60, pm25_total, label = f"Dp < {pm25_bound} $\mu m$")
plt.xlabel("Time (hr)")
plt.ylabel("Concentration (# cm$^{-3}$)")
plt.legend()
plt.show()


salsa_nc[np.linspace(1, len(salsa_nc.columns)-1, len(salsa_nc.columns)-1, dtype = int)].plot()