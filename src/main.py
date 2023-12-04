# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

GITHUB_PATH = os.path.join(os.curdir, "..\\..") # Two levels up, back to the general GitHub directory
PYCHAM_PATH = os.path.join(GITHUB_PATH, "PyCHAM\\PyCHAM\\output\\ex_chem_scheme")

class Simulation:
    def __init__(self, simname):
        self.simname = simname
        filedir = os.listdir(os.path.join(PYCHAM_PATH, self.simname))
        self.output_files = {file: os.path.join(os.path.join(PYCHAM_PATH, self.simname), file) for file in filedir if not file.endswith("npy")}
        self.data = {}
        self.read_output_files_to_dataframes()
        
    def read_output_files_to_dataframes(self):
        for file_name, file_path in self.output_files.items():
            if file_name not in ["inputs", "model_and_component_constants", "simulation_self.pickle"]:
                print(file_name)
                df = pd.read_csv(file_path, delimiter = ",", comment = "#")
                self.data[file_name] = df
            
    def get_data(self, file_name):
        return self.data.get(file_name)
    
    def plot(self, file_name):
        self.get_data(file_name).plot()
    
    def __str__(self):
        return f"{self.simname}\n{list(self.data.keys())}"

sim1 = Simulation("Example_Run_Output_2")

sim1_bins = sim1.get_data("size_bin_bounds")
sim1_nc = sim1.get_data("particle_number_concentration_dry").set_axis(sim1_bins.iloc[-1][1:], axis = 1)

pm01 = sim1_nc[sim1_nc.columns[sim1_nc.columns<=0.1]]
pm25 = sim1_nc[sim1_nc.columns[sim1_nc.columns<=2.5]]
pm01_total = pm01.sum(axis = 1)
pm25_total = pm25.sum(axis = 1)


plt.figure(figsize = (10,6))
plt.imshow(sim1_nc.values.T, cmap = "viridis", aspect = "auto", interpolation = "none")
plt.gca().invert_yaxis()
plt.colorbar()
plt.show()


plt.figure(figsize = (10,6))
plt.plot(pm01_total)
plt.plot(pm25_total, "--")
plt.show()