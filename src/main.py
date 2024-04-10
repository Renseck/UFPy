# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score
from tqdm import tqdm

import HAM_box as ham
import HAM_plot as hp
import measurement_data as md
import utils

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")


if __name__ == "__main__":
    # Start by reading measurement data and filtering/resampling
    smps_df = md.read_measurement_data("SMPS")
    smps_df = md.filter_outliers(smps_df)
    
    smps_bins = smps_df.filter(regex="^[.0-9]+$").columns
    smps_bins_float = [7.] + [float(binbound) for binbound in smps_bins]
    
    weather_df = md.read_measurement_data("Davis")
    weather_df["Spec Hum"] = utils.RH2q(weather_df["In Hum"], weather_df["Bar"]*100, weather_df["Temp Out"] + 273.15)
    
    # Merge weather and SMPS data together
    merged_df = pd.merge_asof(smps_df, weather_df, on = "datetime", direction = "nearest")
    merged_df = merged_df.drop(["time24"], axis = 1)
    merged_df = merged_df.dropna()
    
    salsa_bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[0:8] * 1e9)
    salsa_bins_float = [float(binbound) for binbound in salsa_bins]
    
    nwsw_dist_orig, environmental_data = md.get_nwsw_dist(smps_df, weather_df)
    nwsw_dist = md.translate_particles(smps_bins_float, nwsw_dist_orig.values, salsa_bins_float)
    
    # Read model data
    experiment_name = "test"
    # environmental_data = [298, 0.0058535, 101325]
    ham.write_environmental_data(environmental_data)
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    bin_boundaries = hp.define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    
    num, metadata = ham.read_model_data(experiment_name)
    num5 = pd.read_csv(os.path.join(os.path.join(MODEL_L0_FOLDER, experiment_name), "num_5.dat"), sep = r"\s+")
    metadata = ham.parse_metadata(metadata)
    
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100
    
    dist_to_check = num
    corrs = []
    rmses = []
    rsquareds = []
    for row in tqdm(range(len(dist_to_check)-1), desc = "Calculating optimum time..."):
        rmse = mean_squared_error(dist_to_check[nwsw_dist.columns].iloc[row], nwsw_dist.iloc[0])
        r2 = r2_score(nwsw_dist.iloc[0], dist_to_check[nwsw_dist.columns].iloc[row])
        corr = np.corrcoef(dist_to_check[nwsw_dist.columns].iloc[row], nwsw_dist.iloc[0])[0,1]
        rmses.append(rmse)
        rsquareds.append(r2)
        corrs.append(corr)

    print(f"\nMin RMSE of {np.min(rmses):.2e} at time = {np.argmin(rmses)}s")
    print(f"Max R2 of {np.max(rsquareds):.2f} at time = {np.argmax(rsquareds)}s")
    print(f"Max correlation of {np.max(corrs):.2f} at time = {np.argmax(corrs)}s")

    fig, axes = hp.plot_size_dist(rdry, nwsw_dist, rows = [0], populations = ["a"], ymin = 1, linestyle = "dashed")
    hp.plot_size_dist(rdry, dist_to_check[nwsw_dist.columns], rows = [0, 20, 70, np.argmin(rmses)], populations = ["a"], ymin = 1e6, ymax = 1e12, fig = fig, axes = axes)
    

