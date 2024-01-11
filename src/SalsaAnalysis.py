# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 09:08:27 2023

@author: rens_
"""

# =============================================================================
# This file functions to analyse the output of the *FULL* SALSA2.0 Standalone model, which includes 
# atmospheric chemistry and cloud formation as well. Doumentation is extremely limited, but the repo can
# be found on https://github.com/UCLALES-SALSA/SALSA-standalone. Investigation of the source code (driver.f90) 
# seems to indicate that number concentration are not currently written to output.
# =============================================================================

import os
import pandas as pd 

GITHUB_PATH = os.path.join(os.curdir, "..\\..") # Two levels up, back to the general GitHub directory
SALSA_PATH = os.path.join(GITHUB_PATH, "SALSA-standalone")

class SalsaSimulation:
    def __init__(self):
        self.simname = "SalsaSim"
        self.files = {"output": os.path.join(SALSA_PATH, "output.dat"),
                      "radii": os.path.join(SALSA_PATH, "radii2.dat"),
                      "input": os.path.join(SALSA_PATH, "input.dat")}
        self.data = {}
        self.read_output_files_to_dataframes()
        self.input = self.get_data("input")
        self.output = self.get_data("output")
        self.radii = self.get_data("radii")
        
    def read_output_files_to_dataframes(self):
        for file_name, file_path in self.files.items():
            # print(file_name)
            df = pd.read_fwf(file_path, header = None)
            self.data[file_name] = df
            
    def get_data(self, file_name):
        return self.data.get(file_name)
    
    def plot(self, file_name):
        self.get_data(file_name).plot()
    
    def __str__(self):
        return f"{self.simname}\n{list(self.data.keys())}"