# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 09:08:27 2023

@author: rens_
"""
import os
import pandas as pd 

GITHUB_PATH = os.path.join(os.curdir, "..\\..") # Two levels up, back to the general GitHub directory
SALSA_PATH = os.path.join(GITHUB_PATH, "SALSA-standalone")

class SalsaSimulation:
    def __init__(self):
        self.simname = "SalsaSim"
        self.output_files = {"output": os.path.join(SALSA_PATH, "output.dat"),
                             "radii": os.path.join(SALSA_PATH, "radii2.dat")}
        self.data = {}
        self.read_output_files_to_dataframes()
        self.output = self.get_data("output")
        self.radii = self.get_data("radii")
        
    def read_output_files_to_dataframes(self):
        for file_name, file_path in self.output_files.items():
            # print(file_name)
            df = pd.read_fwf(file_path, header = None)
            self.data[file_name] = df
            
    def get_data(self, file_name):
        return self.data.get(file_name)
    
    def plot(self, file_name):
        self.get_data(file_name).plot()
    
    def __str__(self):
        return f"{self.simname}\n{list(self.data.keys())}"