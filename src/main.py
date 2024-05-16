# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

@author: rens_
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    davis_201 = md.read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    
    davis_641 = md.read_measurement_data("Davis641", show_comments = False)
    davis_641 = davis_641.rename(columns = {"Waarden": "Bar",
					 "Waarden.1": "Rain minute",
					 "Waarden.2": "Rain day",
					 "Waarden.3": "Rain pulse",
					 "Waarden.4": "In Hum",
					 "Waarden.5": "Out Hum", 
					 "Waarden.6": "Temp In", 
					 "Waarden.7": "Temp Out", 
					 "Waarden.8": "winddir_deg", 
					 "Waarden.9": "Wind Speed"})
    davis_641["Spec Hum"] = utils.RH2q(davis_641["Out Hum"], davis_641["Bar"]*100, davis_641["Temp Out"] + 273.15)
    davis_641["datetime"] = pd.to_datetime(davis_641["Begintijd"], format = "%d-%m-%Y %H:%M")
    
    rivm_201 = md.read_measurement_data("SMPS", show_comments = False)
    rivm_201 = md.smps_filter(rivm_201)
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    rivm_641 = md.smps_filter(rivm_641)
    
    smps_bins = rivm_201.filter(regex="^[.0-9]+$").columns
    smps_bins_float = [7.] + [float(binbound) for binbound in smps_bins]
    
    salsa_bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[0:8] * 1e9)
    salsa_bins_float = [float(binbound) for binbound in salsa_bins]
    
# =============================================================================
#     # Merge weather and SMPS data together
#     merged_df = pd.merge_asof(rivm_201, davis_201, on = "datetime", direction = "nearest")
#     merged_df = merged_df.drop(["time24"], axis = 1)
#     merged_df = merged_df.dropna()
# =============================================================================
    
    highway_dist_orig, environmental_data = md.get_directional_dist(rivm_641, davis_641, min_angle = 202.5, max_angle = 337.5)
    background_dist_orig, _ = md.get_directional_dist(rivm_201, davis_201, min_angle = 45, max_angle = 135)
    highway_dist = md.translate_particles(smps_bins_float, highway_dist_orig.values, salsa_bins_float)
    background_dist = md.translate_particles(smps_bins_float, background_dist_orig.values, salsa_bins_float)
    
    # Read model data
    experiment_name = "test"
    # environmental_data = [291, 0.0096788, 100901]
    ham.write_environmental_data([environmental_data]*6)
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    bin_boundaries = hp.define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    
    num, metadata = ham.read_model_data(experiment_name)
    num5 = pd.read_csv(os.path.join(os.path.join(MODEL_L0_FOLDER, experiment_name), "num_5.dat"), sep = r"\s+")
    metadata = ham.parse_metadata(metadata)
    
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100
    
    dist_to_check = num
    
    fig, axes = hp.plot_size_dist(rdry.iloc[1:], highway_dist, rows=[0], ymin=1, xmin = -20e-9, xmax = 400e-9,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "-.", label = "NL10641 hw", color = "black")
    
    fig, axes = hp.plot_size_dist(rdry.iloc[1:], background_dist, rows=[0], ymin=1, xmin = -20e-9, xmax = 400e-9,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      fig = fig, axes = axes, linestyle = "-.", label = "N201 bg", color = "green")
    
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1], ymin=1, xmin = -20e-9, xmax = 400e-9,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      fig = fig, axes = axes, linestyle = "dashed", label = "Model")
    
    hp.plot_size_dist(rdry, num5*1e-6, rows=[500], ymin=1, ymax = 2e3, xmin = -20e-9, xmax = 400e-9,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      fig = fig, axes = axes, linestyle = "-", label = "Model")
    

