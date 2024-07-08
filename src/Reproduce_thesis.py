# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:37:39 2024

The entire point of this file is simply to recreate all data and figures shown in the master's thesis that was written
with this project. All the functions in it aim to generate one set of results each, which will all be collected into 
one 'main' function, which is then capable of executing all of them in one go. The running time of this will be significant,
and it really only serves a purpose in the interest of open science. I'm not going to bother commenting on each and 
every function either. Just run the file, it'll take an age and reproduce everything seen in the thesis.

@author: rens_
"""
import os
from itertools import product

import cmocean as co
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

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
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")
smps_bins_float = [3.0, 11.5, 15.4, 20.5, 27.4, 36.5, 48.7, 64.9, 86.6, 115.5, 154.0, 205.4, 273.8, 365.2]

def show_windroses():
    davis_201 = md.read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    # plot_windrose(davis_201, "Wind Speed")
    md.plot_windrose(davis_201, "Wind Speed", title = "Wind speed (N201)", min_angle = 202.5, max_angle = 270)
    
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
    # plot_windrose(davis_641, "Wind Speed")
    md.plot_windrose(davis_641, "Wind Speed", title = "Wind speed (NL10641)", min_angle = 202.5, max_angle = 337.5)

def show_directional_dists():    
    # Read data, add spechum and make time usable
    davis_201 = md.read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    
    rivm_201 = md.read_measurement_data("SMPS", show_comments = False)
    rivm_201 = md.smps_filter(rivm_201)
    
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
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    rivm_641 = md.smps_filter(rivm_641)
    
    # Calculate directional distributions
    highway_dist_201, _ = md.get_directional_dist(rivm_201, davis_201, min_angle = 202.5, max_angle = 270)
    background_dist_201, _ = md.get_directional_dist(rivm_201, davis_201, min_angle = 90, max_angle = 180)
    
    highway_dist_641, _ = md.get_directional_dist(rivm_641, davis_641, min_angle = 202.5, max_angle = 337.5)
    background_dist_641, _ = md.get_directional_dist(rivm_641, davis_641, min_angle = 45, max_angle = 135)

    # Show the various mean distribution of particles
    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    md.plot_data(highway_dist_201, ax = axd["main"], label = "W - SSW", title="Mean size distribution (N201)")
    md.plot_data(background_dist_201, ax = axd["main"], label = "E-S")
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/N201_Dist.jpg"), bbox_inches = "tight", dpi = 150)
    plt.show()
    
    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    md.plot_data(highway_dist_641, ax = axd["main"], label = "NNW - SSW", title="Mean size distribution (NL10641)")
    md.plot_data(background_dist_641, ax = axd["main"], label = "NE-SE")
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/N641_Dist.jpg"), bbox_inches = "tight", dpi = 150)
    plt.show()
    
def show_resampled_dists():
    salsa_bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[0:8] * 1e9)
    salsa_bins_float = [float(binbound) for binbound in salsa_bins]
    salsa_upper_boundaries = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[1:8] * 1e9)
    
    # Read data, add spechum and make time usable
    davis_201 = md.read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    
    rivm_201 = md.read_measurement_data("SMPS", show_comments = False)
    rivm_201 = md.smps_filter(rivm_201)
    
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
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    rivm_641 = md.smps_filter(rivm_641)
    
    # Calculate directional distributions
    highway_dist_201, _ = md.get_directional_dist(rivm_201, davis_201, min_angle = 202.5, max_angle = 270)
    background_dist_201, _ = md.get_directional_dist(rivm_201, davis_201, min_angle = 90, max_angle = 180)
    
    highway_dist_641, _ = md.get_directional_dist(rivm_641, davis_641, min_angle = 202.5, max_angle = 337.5)
    background_dist_641, _ = md.get_directional_dist(rivm_641, davis_641, min_angle = 45, max_angle = 135)
    
    highway_201_resampled = md.translate_particles(smps_bins_float[1:], highway_dist_201, salsa_bins_float[1:])
    highway_641_resampled = md.translate_particles(smps_bins_float[1:], highway_dist_641, salsa_bins_float[1:])
        
    # Show the resampling of the nw_sw distribution
    fig = plt.figure(layout="tight", figsize=(16, 6))
    fig.suptitle("Resampled distributions", x = 0.52, fontsize = 15)
    fig.text(0.49, -0.02, "Diameter (nm)", fontsize = 15)

    axd = fig.subplot_mosaic([["left", "right"]], sharey = True, sharex = True)
    axd["right"].plot(salsa_upper_boundaries, highway_201_resampled.values[0], label = "SALSA2.0")
    axd["right"].plot(highway_dist_201.index.astype(float), highway_dist_201.values, label = "SMPS")
    
    axd["right"].set_title("N201", fontsize = 14)
    axd["right"].tick_params(axis = "both", labelsize = 12)
    axd["right"].legend(fontsize = 14)

    axd["left"].plot(salsa_upper_boundaries, highway_641_resampled.values[0], label = "SALSA2.0")
    axd["left"].plot(highway_dist_641.index.astype(float), highway_dist_641.values, label = "SMPS")
    highway_641_resampled.values[0][1] = 1760
    axd["left"].plot(salsa_upper_boundaries, highway_641_resampled.values[0], label = "Corrected", color = "green")
    axd["left"].set_ylabel("# particles cm$^{-3}$", fontsize = 15)
    axd["left"].set_title("NL10641", fontsize = 14)
    axd["left"].tick_params(axis = "both", labelsize = 12)
    axd["left"].legend(fontsize = 14)
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/resampling_subfig.png"), bbox_inches = "tight", dpi = 150)
    plt.show()

def show_correlations():
    # Read data, add spechum and make time usable
    davis_201 = md.read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    
    rivm_201 = md.read_measurement_data("SMPS", show_comments = False)
    rivm_201 = md.smps_filter(rivm_201)
    
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
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    rivm_641 = md.smps_filter(rivm_641)
    
    merged_201_df = pd.merge_asof(rivm_201, davis_201, on = "datetime", direction = "nearest")
    merged_201_df = merged_201_df.drop(["time24"], axis = 1)
    merged_201_df = merged_201_df.dropna()
    
    merged_641_df = pd.merge_asof(rivm_641, davis_641, on = "datetime", direction = "nearest")
    merged_641_df = merged_641_df.dropna()
    
    # Plot (immediate) correlations
    column_order = ['11.5', '15.4', '20.5', '27.4', '36.5', '48.7', '64.9',
                    '86.6', '115.5','154.0', '205.4', '273.8', '365.2',
                    'winddir_deg', 'Wind Speed', 'Out Hum', 'Temp Out','Spec Hum']
    
    plt.figure(figsize = (10,6))
    plt.title("Correlations N201")
    corr_201 = merged_201_df.filter(regex = "(?i)^[.0-9]+$|Wind|Temp|Hum|Status").\
           drop(["Wind Dir", "Wind Tx", "Wind Samp", "Wind Chill", "Wind Run",
                 "In  Temp", "Hi Temp", "Low Temp", "In Hum"], axis = 1).query("Status == 'No errors'").\
           drop(["Status"], axis = 1).reindex(columns = column_order).rename(columns = {"winddir_deg": "Wind Direction"}).corr()
           
    cmap = sb.diverging_palette(-1, 1, as_cmap = True)
    sb.heatmap(corr_201, cmap = co.cm.balance_r, mask = np.triu(corr_201), center = 0)
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/measurment_correlations_N201.jpg"), dpi = 150)

    # Plot (immediate) correlations
    plt.figure(figsize = (10,6))
    plt.title("Correlations NL10641")
    corr_641 = merged_641_df.filter(regex = "(?i)^[.0-9]+$|Wind|Temp|Hum|Status").\
           drop([col for col in merged_641_df.columns if "Status " in col], axis = 1).query("Status == 'No errors'").\
           drop(["Status", "In Hum", "Temp In"], axis = 1).reindex(columns = column_order).rename(columns = {"winddir_deg": "Wind Direction"}).corr()
           
           
    cmap = sb.diverging_palette(-1, 1, as_cmap = True)
    sb.heatmap(corr_641, cmap = co.cm.balance_r, mask = np.triu(corr_641), center = 0)
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/measurment_correlations_N641.jpg"), dpi = 150)

def diesel_run():
    experiment_name = "Diesel_emissions"
    
    environmental_data = [293, 0.0096610, 99890]
    ham.write_environmental_data([environmental_data]*6)
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.68
    scale = 17e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    diesel_lognormal = utils.lognormal(x, sigma, center_x = 90, scale = scale)
    diesel_flux = np.interp(salsa_bins, x, diesel_lognormal)
    
    particle_flux = diesel_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.012
    ham.write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)
    
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    num, metadata = ham.read_model_data(experiment_name, gridcell = 1)
    num5, _ = ham.read_model_data(experiment_name, gridcell = 5)
    metadata = ham.parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    fig, axes = hp.plot_size_dist(rdry, num, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_close = [1669, 1981, 816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                      fig = fig, axes = axes, label = "Model")
    
    # axes.plot(smps_bins_float[1:], smps_close, label = "60m", linestyle = "-.")
    axes.plot(smps_bins_float[1:], smps_far, label = "Measurement (320m)", linestyle = "-.")
    axes.legend()
    fig_name = "size_distribution.png"
    savepath = os.path.join(MODEL_PLOT_FOLDER, experiment_name)
    full_savepath = os.path.join(savepath, fig_name)
    plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
    plt.show()
    
    hp.stacked_timeseries_plot(num, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 1)")
    fig, axes = hp.stacked_timeseries_plot(num5, populations = ["a"], ymin = 1, exp_name = experiment_name,
                                           title = "Size distribution evolution (cell 5)",
                                           highlights = [600, 1500], highlight_colors = ["green", "red"])
    
def diesel_run_secondaries():
    experiment_name = "Diesel_emissions_with_secondary"
    
    environmental_data = [293, 0.0096610, 99890]
    ham.write_environmental_data([environmental_data]*6)
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.68
    scale = 17e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    diesel_lognormal = utils.lognormal(x, sigma, center_x = 90, scale = scale)
    diesel_flux = np.interp(salsa_bins, x, diesel_lognormal)
    
    secondary_lognormal = utils.lognormal(x, sigma, center_x = 20, scale = scale/2)
    additional_flux = np.interp(salsa_bins, x, secondary_lognormal)
    
    particle_flux = diesel_flux + additional_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.012
    ham.write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)
    
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    num, metadata = ham.read_model_data(experiment_name, gridcell = 1)
    num5, _ = ham.read_model_data(experiment_name, gridcell = 5)
    metadata = ham.parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    fig, axes = hp.plot_size_dist(rdry, num, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_close = [1669, 1981, 816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                      fig = fig, axes = axes, label = "Model")
    
    # axes.plot(smps_bins_float[1:], smps_close, label = "60m", linestyle = "-.")
    axes.plot(smps_bins_float[1:], smps_far, label = "Measurement (320m)", linestyle = "-.")
    axes.legend()
    fig_name = "size_distribution.png"
    savepath = os.path.join(MODEL_PLOT_FOLDER, experiment_name)
    full_savepath = os.path.join(savepath, fig_name)
    plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
    plt.show()
    
    hp.stacked_timeseries_plot(num, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 1)")
    fig, axes = hp.stacked_timeseries_plot(num5, populations = ["a"], ymin = 1, exp_name = experiment_name,
                                           title = "Size distribution evolution (cell 5)",
                                           highlights = [600, 1500], highlight_colors = ["green", "red"])
    
    
def gasoline_run():
    experiment_name = "gasoline_emissions"
    
    environmental_data = [293, 0.0096610, 99890]
    ham.write_environmental_data([environmental_data]*6)
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.7
    scale = 8e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    main_lognormal = utils.lognormal(x, sigma, center_x = 20, scale = scale)
    main_flux = np.interp(salsa_bins, x, main_lognormal)

    particle_flux = main_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.006
    ham.write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)
    
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    num, metadata = ham.read_model_data(experiment_name, gridcell = 1)
    num5, _ = ham.read_model_data(experiment_name, gridcell = 5)
    metadata = ham.parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    fig, axes = hp.plot_size_dist(rdry, num, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                      fig = fig, axes = axes, label = "Model")
    
    # axes.plot(smps_bins_float[1:], smps_close, label = "60m", linestyle = "-.")
    axes.plot(smps_bins_float[1:], smps_far, label = "Measurement (320m)", linestyle = "-.")
    axes.legend()
    fig_name = "size_distribution.png"
    savepath = os.path.join(MODEL_PLOT_FOLDER, experiment_name)
    full_savepath = os.path.join(savepath, fig_name)
    plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
    plt.show()
    
    hp.stacked_timeseries_plot(num, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 1)")
    fig, axes = hp.stacked_timeseries_plot(num5, populations = ["a"], ymin = 1, exp_name = experiment_name,
                                           title = "Size distribution evolution (cell 5)",
                                           highlights = [600, 1500], highlight_colors = ["green", "red"])
    
def gasoline_run_secondaries():
    experiment_name = "gasoline_emissions_with_background"
    
    environmental_data = [293, 0.0096610, 99890]
    ham.write_environmental_data([environmental_data]*6)
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.7
    scale = 8e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    main_lognormal = utils.lognormal(x, sigma, center_x = 20, scale = scale)
    main_flux = np.interp(salsa_bins, x, main_lognormal)
    
    secondary_flux = np.array([0., 0., 0., 1355., 1271., 177., 0., .0, 0.])

    particle_flux = main_flux + 3000*secondary_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.006
    ham.write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)
    
    ham.run_model(experiment_name=experiment_name, recompile=True)
    
    num, metadata = ham.read_model_data(experiment_name, gridcell = 1)
    num5, _ = ham.read_model_data(experiment_name, gridcell = 5)
    metadata = ham.parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    fig, axes = hp.plot_size_dist(rdry, num, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                      fig = fig, axes = axes, label = "Model")
    
    # axes.plot(smps_bins_float[1:], smps_close, label = "60m", linestyle = "-.")
    axes.plot(smps_bins_float[1:], smps_far, label = "Measurement (320m)", linestyle = "-.")
    axes.legend()
    fig_name = "size_distribution.png"
    savepath = os.path.join(MODEL_PLOT_FOLDER, experiment_name)
    full_savepath = os.path.join(savepath, fig_name)
    plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
    plt.show()
    
    hp.stacked_timeseries_plot(num, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 1)")
    fig, axes = hp.stacked_timeseries_plot(num5, populations = ["a", "b"], ymin = 1, exp_name = experiment_name,
                                           title = "Size distribution evolution (cell 5)",
                                           highlights = [600, 1500], highlight_colors = ["green", "red"])

def temp_humi_variation():
    temps = np.linspace(273, 298, num = 15)
    humis = np.linspace(0, 0.006, num = 15)
    pap = 101325
    
    environmental_values = [[temp, humi, pap] for temp, humi in product(temps, humis)]
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.68
    scale = 17e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    diesel_lognormal = utils.lognormal(x, sigma, center_x = 90, scale = scale)
    diesel_flux = np.interp(salsa_bins, x, diesel_lognormal)
    
    secondary_lognormal = utils.lognormal(x, sigma, center_x = 20, scale = scale/2)
    additional_flux = np.interp(salsa_bins, x, secondary_lognormal)
    
    particle_flux = diesel_flux + additional_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.012

    experiment_name = f"Diesel_Flux_TempHumidityVariation{pap}"
    
    ham.run_variation(experiment_name, environmental_values, [particle_flux], [dispersion_rate])
    
    numdict, metadict = ham.read_model_data(experiment_name, gricell = 5)
    hp.plot_variation_surface(numdict, metadict, experiment_name, "1a1", 1500)
    
def temp_pressure_variation():
    temps = np.linspace(273, 298, num = 15)
    humi = 0.0033
    paps = np.linspace(97000, 105000, num = 15)
    
    environmental_values = [[temp, humi, pap] for temp, pap in product(temps, paps)]
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.68
    scale = 17e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    diesel_lognormal = utils.lognormal(x, sigma, center_x = 90, scale = scale)
    diesel_flux = np.interp(salsa_bins, x, diesel_lognormal)
    
    secondary_lognormal = utils.lognormal(x, sigma, center_x = 20, scale = scale/2)
    additional_flux = np.interp(salsa_bins, x, secondary_lognormal)
    
    particle_flux = diesel_flux + additional_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.012

    experiment_name = f"Diesel_Flux_TempPressureVariation{humi}"
    
    ham.run_variation(experiment_name, environmental_values, [particle_flux], [dispersion_rate])
    
    numdict, metadict = ham.read_model_data(experiment_name, gricell = 5)
    hp.plot_variation_surface(numdict, metadict, experiment_name, "1a1", 1500)

def pressure_vs_particle_count():
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
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    # rivm_641 = filter_outliers(rivm_641)
    rivm_641 = md.smps_filter(rivm_641)
    
    merged_641_df = pd.merge_asof(rivm_641, davis_641, on = "datetime", direction = "nearest")
    merged_641_df = merged_641_df.dropna()
    
    directional_df = merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]
    
    coeffs = np.polyfit(directional_df["Spec Hum"], directional_df["11.5"], 1)

    plt.figure(figsize = (10,6))
    plt.scatter(directional_df["Spec Hum"], directional_df["11.5"], facecolor = 'none', edgecolor = 'blue',
                label = "Measurement")
    plt.plot(directional_df["Spec Hum"], directional_df["Spec Hum"]*coeffs[0] + coeffs[1],
             color = "orange", label = f"Regression \n(slope = {coeffs[0]:.2f})")
    plt.ylim(0, 1e4)

    plt.title("Relation between ambient pressure and particle counts ($\leq 11.5$ nm)")
    plt.xlabel("q ($\\frac{kg}{kg}$)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.legend()
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/pressure_vs_11.5_NL641.jpg"), dpi = 150)

def spechum_vs_particle_count():
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
    
    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    # rivm_641 = filter_outliers(rivm_641)
    rivm_641 = md.smps_filter(rivm_641)
    
    merged_641_df = pd.merge_asof(rivm_641, davis_641, on = "datetime", direction = "nearest")
    merged_641_df = merged_641_df.dropna()
    
    coeffs = np.polyfit(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Spec Hum"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["11.5"], 1)

    plt.figure(figsize = (10,6))
    plt.scatter(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Spec Hum"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["11.5"], facecolor = 'none', edgecolor = 'blue', label = "Measurement")
    plt.plot(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Spec Hum"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Spec Hum"]*coeffs[0] + coeffs[1], color = "orange", label = f"Regression \n(slope = {coeffs[0]:.2f})")
    plt.ylim(0, 1e4)

    plt.title("Relation between specific humidity and particle counts ($\leq 11.5$ nm)")
    plt.xlabel("q ($\\frac{kg}{kg}$)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.legend()
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/spechum_vs_11.5_NL641.jpg"), dpi = 150)

def temp_vs_particle_count():
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

    rivm_641 = md.read_measurement_data("RIVM", show_comments = False)
    # rivm_641 = filter_outliers(rivm_641)
    rivm_641 = md.smps_filter(rivm_641)

    merged_641_df = pd.merge_asof(rivm_641, davis_641, on = "datetime", direction = "nearest")
    merged_641_df = merged_641_df.dropna()

    coeffs = np.polyfit(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Temp Out"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["11.5"], 1)

    plt.figure(figsize = (10,6))
    plt.scatter(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Temp Out"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["11.5"], facecolor = 'none', edgecolor = 'blue', label = "Measurement")
    plt.plot(merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Temp Out"], merged_641_df[(merged_641_df["winddir_deg"] >= 202.5) &(merged_641_df["winddir_deg"] <= 337.5)]["Temp Out"]*coeffs[0] + coeffs[1], color = "orange", label = f"Regression \n(slope = {coeffs[0]:.2f})")
    plt.ylim(0, 1e4)

    plt.title("Relation between temperature and particle counts ($\leq 11.5$ nm)")
    plt.xlabel("T ($C^{\\circ}$)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.legend()
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/temp_vs_11.5_NL641.jpg"), dpi = 150)

def dispersion_variation_2():
    exp_name = "Diesel_flux_dispersion_variation_2"
    environmentals = [298, 0.005, 101325]
    
    particle_flux = [3.60323517e-03, 6.28993054e+01, 4.82238719e+04, 1.62828747e+06,
           2.97977956e+06, 1.16208335e+06, 9.65807270e+04, 1.71058143e+03,
           1.01516992e+02]
           
    dispersions = np.linspace(0, 0.02, num = 10)
    
    ham.run_variation(exp_name, environmentals, [particle_flux], dispersions)
    numdict, metadict = ham.read_model_data(exp_name, gridcell = 5)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")
    
    hp.plot_dispersion_variation(numdict, metadict, rdry, time = 1500)
    
def dispersion_variation_3():
    exp_name = "Diesel_flux_dispersion_variation_3"
    environmentals = [298, 0.005, 101325]
    
    particle_flux = [3.60323517e-03, 6.28993054e+01, 4.82238719e+04, 1.62828747e+06,
           2.97977956e+06, 1.16208335e+06, 9.65807270e+04, 1.71058143e+03,
           1.01516992e+02]
           
    dispersions = np.linspace(0, 0.005, num = 20)
    
    ham.run_variation(exp_name, environmentals, [particle_flux], dispersions)
    numdict, metadict = ham.read_model_data(exp_name, gridcell = 5)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")
    
    hp.plot_dispersion_variation(numdict, metadict, rdry, time = 1500)

def main():
    """
    Be incredibly careful with running this function, because depending on your device,
    it may take quite a long time and approximately 6.5 gigabytes of storage before it's completed.
    """
    utils.print_underlined("[WARNING]")
    print("You're about to start the function that reproduces all data and figures shown in the thesis.\nThis will take a considerable amount of time (approx. 12 minutes on an AMD Ryzen 5 5600X 6-Core Processor 3.70 GHz with plenty of RAM), and aroud \033[4m6.5 gigabytes\033[0m of storage. Are you sure you want to continue?\n")
    
    run = input("Continue? yes/no\n")
    
    if run in ["yes", "y", "ye", "yep", "sure", "go"]:
        print("Good luck!")
        # Methods
        
        utils.print_underlined("Figure 7")
        hp.show_normalised_lognormal_flux(sigma = 1.68, center_x = 90, scale = 1,
                                          title = "Normalized lognormal diesel emissions", label = "", file_addition = "normalized")
        
        utils.print_underlined("Figure 9")
        show_windroses()
        
        utils.print_underlined("Figure 10")
        show_directional_dists()
        
        utils.print_underlined("Figure 11")
        show_resampled_dists()
        
        utils.print_underlined("Figure 12")
        show_correlations()
        
        utils.print_underlined("Figure 14")
        md.show_haarrijn_data()
        
        utils.print_underlined("Figure 15")
        md.show_haarrijn_init()
        
        # Results
        
        utils.print_underlined("Figure 16")
        hp.show_model_lognormal_flux(sigma = 1.68, center_x = 90, scale = 17, label = "", file_addition = "diesel")
        
        utils.print_underlined("Figures 17 and 18")
        diesel_run()
        
        utils.print_underlined("Figure 19")
        diesel_run_secondaries()
        
        utils.print_underlined("Figure 20")
        hp.show_model_lognormal_flux(sigma = 1.7, center_x = 20, scale = 17, label = "gasoline", file_addition = "gasoline")
        
        utils.print_underlined("Figures 21 and 22")
        gasoline_run()
        
        utils.print_underlined("Figure 23")
        hp.show_model_flux(np.array([2.56258763e+06, 2.56258763e+06, 6.00964113e+06, 5.41928197e+06,
               3.88603985e+06, 5.31839479e+05, 2.05616901e+00, 1.07327004e-03,
               9.47644140e-06])*1e-6)
        
        utils.print_underlined("Figures 24 and 25")
        gasoline_run_secondaries()
        
        utils.print_underlined("Figure 26")
        temp_humi_variation()
        
        utils.print_underlined("Figure 27")
        spechum_vs_particle_count()
        
        utils.print_underlined("Figure 28")
        temp_vs_particle_count()()
        
        utils.print_underlined("Figure 29")
        temp_pressure_variation()
        
        utils.print_underlined("Figure 30")
        pressure_vs_particle_count()
        
        utils.print_underlined("Figure 31, first image")
        dispersion_variation_2()
        
        utils.print_underlined("Figure 32, second image")
        dispersion_variation_3()
        
        # Appendix
        
        utils.print_underlined("Figure 33, second image")
        rivm_201 = md.read_measurement_data("SMPS", show_comments = False)
        rivm_201 = md.smps_filter(rivm_201)
        md.show_bin_difference(rivm_201)
        
    else:
        print("Alright, stopping the process.")
        
if __name__ == "__main__":
    main()