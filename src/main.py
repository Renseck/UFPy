# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from sklearn.metrics import r2_score, mean_squared_error
import matplotlib.pyplot as plt
import measurement_data as md
import HAM_box as ham
import HAM_plot as hp
from utils import *

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")


if __name__ == "__main__":
    # Start by reading measurement data and filtering/resampling
    smps_df = md.read_measurement_data("SMPS")
    weather_df = md.read_measurement_data("Davis")
    
    nwsw_dist = md.get_nwsw_dist(smps_df, weather_df)
    
    # Read model data
    experiment_name = "test"
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    bin_boundaries = hp.define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    
    num, metadata = ham.read_model_data(experiment_name)
    metadata = ham.parse_metadata(metadata)
    
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100
    
    corrs = []
    rmses = []
    rsquareds = []
    for row in tqdm(range(len(num)-1), desc = "Calculating optimum time..."):
        rmse = mean_squared_error(num[nwsw_dist.columns].iloc[row], nwsw_dist.iloc[0])
        r2 = r2_score(nwsw_dist.iloc[0], num[nwsw_dist.columns].iloc[row])
        corr = np.corrcoef(num[nwsw_dist.columns].iloc[row], nwsw_dist.iloc[0])[0,1]
        rmses.append(rmse)
        rsquareds.append(r2)
        corrs.append(corr)

    print(f"\nMin RMSE of {np.min(rmses):.2f} at time = {np.argmin(rmses)}s")
    print(f"Max R2 of {np.max(rsquareds):.2f} at time = {np.argmax(rsquareds)}s")
    print(f"Max correlation of {np.max(corrs):.2f} at time = {np.argmax(corrs)}s")

    fig, axes = hp.plot_size_dist(rdry, nwsw_dist, rows = [0], populations = ["a"], ymin = 1, linestyle = "dashed")
    hp.plot_size_dist(rdry, num[nwsw_dist.columns], rows = [0, 20, 40, 2720], populations = ["a"], ymin = 1e6, ymax = 1e12, fig = fig, axes = axes)


