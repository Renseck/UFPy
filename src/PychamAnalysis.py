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
    
if __name__ == "__main__":
    pycham0 = PychamSimulation("Example_Run_Output_2")
    pycham1 = PychamSimulation("Slow_Chemical_Reactions")
    
    pycham0_bins = pycham0.get_data("size_bin_bounds")
    pycham0_nc = pycham0.get_data("particle_number_concentration_dry").set_axis(
        pycham0_bins.iloc[-1][1:], axis=1
    )
    
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