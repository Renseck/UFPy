# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 09:12:34 2023

Showcases the workflow.

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

# Some useful data paths and stuff
HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")
smps_bins_float = [3.0, 11.5, 15.4, 20.5, 27.4, 36.5, 48.7, 64.9, 86.6, 115.5, 154.0, 205.4, 273.8, 365.2]


#%%
# =============================================================================
# Input of the model
# =============================================================================

# Give the simulation a name
experiment_name = "Test_1"

# Declare the environmental variables
temp = 290 # Kelvin
humi = 0.005 # kg water per kg air, specific humidity
pressure = 101000 # Pascal

environmental_variables = [temp, humi, pressure] # Order is important

# Particle influx and dispersion

# As it stands, particle_flux is a list with length 9, and dispersion is just a float < 1
# Let's showcase a lognormal input, which needs to be interpolated to salsa's bins
x = np.linspace(3, 1000, 10000) 
lognormal_flux = utils.lognormal(x, sigma = 1.68, center_x = 90, scale = 17e6) # 1e6 for cm^-3 -> m^-3

# Interpolation
bin_boundaries = hp.define_bin_boundaries()
salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)

particle_flux = np.interp(salsa_bins, x, lognormal_flux)
dispersion_rate = 0.012

#%%
# =============================================================================
# Write the input to the data file, and run the model
# =============================================================================

# These need to be written to every horizontal grid-cell in the system
# The best way of doing a homogenous system, is just by copying the list n-times, for n grid cells
# Right now, it's setup for 6, so
nested_environmental_variables = [environmental_variables]*6
ham.write_environmental_data(nested_environmental_variables)

ham.write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)

# recompile can be set to False if you're positive you've not changed the f90 source code of the model
ham.run_model(experiment_name=experiment_name, recompile=False)

#%%
# =============================================================================
# Read the model data
# =============================================================================

# The read_model_data() function can discern between grid cell
# Right now, the model outputs data for cells 1 and 5
num, metadata = ham.read_model_data(experiment_name, gridcell = 1)
num5, _ = ham.read_model_data(experiment_name, gridcell = 5)

# Some plotting functions require the rdry file, which contains the immediate bin boundaries.
# I've not implemented a function for it, because it's a one-liner
rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

# rdry has radii which are off by 2 orders of magnitudes, because SALSA works
# with cm, "for some reason". Divide everything by 100 to make it SI compliant.
rdry = rdry/100

# Parse metadata into more legible form
metadata = utils.parse_metadata(metadata)

# %%
# =============================================================================
# Plotting the data in question
# =============================================================================

# Num files are always in m^-3, but these plot functions scale them down to cm^-3 automatically
fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
                  exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                  label = "Model")

fig, axes = hp.stacked_timeseries_plot(num5, populations = ["a"], ymin = 1, exp_name = experiment_name,
                                       title = "Size distribution evolution (cell 5)",
                                       highlights = [600, 1500], highlight_colors = ["green", "red"])