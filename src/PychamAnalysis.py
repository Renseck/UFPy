# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 09:07:21 2023

@author: rens_
"""
import os
import pandas as pd 

GITHUB_PATH = os.path.join(os.curdir, "..\\..") # Two levels up, back to the general GitHub directory
PYCHAM_PATH = os.path.join(GITHUB_PATH, "PyCHAM\\PyCHAM\\output\\ex_chem_scheme")

class PychamSimulation:
    def __init__(self, simname):
        self.simname = simname
        filedir = os.listdir(os.path.join(PYCHAM_PATH, self.simname))
        self.output_files = {file: os.path.join(os.path.join(PYCHAM_PATH, self.simname), file) for file in filedir if not file.endswith("npy")}
        self.data = {}
        self.read_output_files_to_dataframes()
        
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
    
    def __str__(self):
        return f"{self.simname}\n{list(self.data.keys())}"