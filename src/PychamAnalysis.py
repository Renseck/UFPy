# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 09:07:21 2023

@author: rens_
"""
import os
import pandas as pd 
import matplotlib.pyplot as plt

GITHUB_PATH = os.path.join(os.curdir, "..\\..") # Two levels up, back to the general GitHub directory
PYCHAM_PATH = os.path.join(GITHUB_PATH, "PyCHAM\\PyCHAM\\output\\ex_chem_scheme")

class PychamSimulation:
    def __init__(self, simname):
        self.simname = simname
        filedir = os.listdir(os.path.join(PYCHAM_PATH, self.simname))
        self.output_files = {file: os.path.join(os.path.join(PYCHAM_PATH, self.simname), file) for file in filedir if not file.endswith("npy")}
        self.data = {}
        self.read_output_files_to_dataframes()
        self.bins = self.get_data("size_bin_bounds")
        self.nc = self.get_data("particle_number_concentration_dry").set_axis(
            self.bins.iloc[-1][1:], axis=1
        )
        
    def read_output_files_to_dataframes(self):
        for file_name, file_path in self.output_files.items():
            if file_name not in ["inputs", "model_and_component_constants", "simulation_self.pickle"]:
                # print(file_name)
                df = pd.read_csv(file_path, delimiter = ",", comment = "#")
                self.data[file_name] = df
            
    def get_data(self, file_name):
        return self.data.get(file_name)
    
    def plot(self, file_name):
        self.get_data(file_name).plot()
        
    def plume_plot(self):
        plt.figure(figsize=(10, 6))
        plt.imshow(self.nc.values.T, cmap="viridis", aspect="auto", interpolation="none")
        plt.gca().invert_yaxis()
        plt.colorbar()
        plt.show()
        
    def pnc_plot(self, ax = None, label = None):
        pm25 = self.nc[self.nc.columns[self.nc.columns <= 2.5]]
        pm25_total = pm25.sum(axis = 1)
        
        if ax == None:
            plt.figure(figsize=(10, 6))
            #plt.plot(pm01_total)
            plt.plot(pm25_total.index/60, pm25_total)
            plt.xlabel("Time (hr)")
            plt.ylabel("Concentration (# cm$^{-3}$)")
            plt.show()
        else: 
            ax.plot(pm25_total.index/60, pm25_total, label = label)
            ax.set_xlabel("Time (hr)")
            ax.set_ylabel("Concentration (# cm$^{-3}$)")
            ax.legend()
    
    def __str__(self):
        return f"{self.simname}\n{list(self.data.keys())}"