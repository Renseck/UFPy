# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 13:37:39 2024

The entire point of this file is simply to recreate all data and figures shown in the master's thesis that was written
with this project. All the functions in it aim to generate one set of results each, which will all be collected into 
one 'main' file, which is then capable of executing all of them in one go. The running time of this will be significant,
and it really only serves a purpose in the interest of open science. 

@author: rens_
"""
import os

import cmocean as co
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb

import HAM_plot as hp
import HAM_box as ham
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

def show_methodology():
    md.show_bin_difference(rivm_201)
    show_model_flux(sigma = 1.68, center_x = 90, scale = 1)
    show_haarrijn_data()
    show_haarrijn_init()
    
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
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_close = [1669, 1981, 816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5*1e-6, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
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
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_close = [1669, 1981, 816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5*1e-6, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
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
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5*1e-6, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
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
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    
    smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
    
    fig, axes = hp.plot_size_dist(rdry, num5*1e-6, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
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

    plt.title("Relation between ambient pressure and particle counts (<= 11.5nm)")
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

    plt.title("Relation between specific humidity and particle counts (<= 11.5nm)")
    plt.xlabel("q ($\\frac{kg}{kg}$)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.legend()
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/spechum_vs_11.5_NL641.jpg"), dpi = 150)

def dispersion_variation_1():
    exp_name = "Diesel_flux_dispersion_variation"
    environmentals = [298, 0.005, 101325]
    
    particle_flux = [3.60323517e-03, 6.28993054e+01, 4.82238719e+04, 1.62828747e+06,
           2.97977956e+06, 1.16208335e+06, 9.65807270e+04, 1.71058143e+03,
           1.01516992e+02]
           
    dispersions = np.linspace(0, 0.014, num = 15)
    
    ham.run_variation(exp_name, environmentals, [particle_flux], dispersions)
    numdict, metadict = ham.read_model_data(exp_name, gridcell = 5)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")
    
    hp.plot_dispersion_variation(numdict, metadict, rdry, time = 1500)

def main():
    # Methods
    hp.show_normalised_lognormal_flux(sigma = 1.68, center_x = 90, scale = 1, label = "diesel", file_addition = "normalized")
    hp.show_model_lognormal_flux(sigma = 1.68, center_x = 90, scale = 17, label = "gasoline", file_addition = "diesel")
    hp.show_model_lognormal_flux(sigma = 1.7, center_x = 20, scale = 17, label = "gasoline", file_addition = "gasoline")
    
    # Results
    diesel_run()
    diesel_run_secondaries()
    gasoline_run()
    