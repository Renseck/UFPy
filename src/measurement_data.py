# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 08:20:04 2024

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

DATA_FOLDER = "../data"
MEASUREMENTS_FOLDER = os.path.join(DATA_FOLDER, "Measurements")
RESULTS_FOLDER = "../results"

cabauw_df = pd.read_csv(os.path.join(MEASUREMENTS_FOLDER, "Cabauw_UFP_Hourly_Avgs.csv"), skiprows = 6, header = 1)
erzeij_df = pd.read_csv(os.path.join(MEASUREMENTS_FOLDER, "Erzeijstraat_UFP_Hourly_Avgs.csv"), skiprows = 6, header = 1)
grift_df = pd.read_csv(os.path.join(MEASUREMENTS_FOLDER, "Griftpark_UFP_Hourly_Avgs.csv"), skiprows = 6, header = 1)