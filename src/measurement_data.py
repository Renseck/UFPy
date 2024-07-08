# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 08:20:04 2024

@author: rens_
"""
import os

import cmocean as co
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sb
from windrose import WindroseAxes

import utils
from HAM_plot import define_bin_boundaries, show_model_lognormal_flux
from utils import complementary_color, order_of_magnitude

DATA_FOLDER = "../data"
MEASUREMENTS_FOLDER = os.path.join(DATA_FOLDER, "Measurements")
RESULTS_FOLDER = "../results"

smps_bins = [11.5, 15.4, 20.5, 27.4, 36.5, 48.7, 64.9, 86.6, 115.5, 154.0, 205.4, 273.8, 365.2]
smps_60m = [1669, 1981,	816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
smps_320m = [418, 622,	429, 658, 990, 1292, 1455, 1391, 1026, 481, 13, 0, 0]
smps_200m = [474, 634,	431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]
smps_220m = [566, 851,	526, 673, 975, 1287, 1484, 1445, 1073, 498, 7, 0, 0]

def plot_data(dataframe, ax=None, title="", label="",
              xlabel="Diameter (nm)", ylabel="# particles cm${^-3}$",
              linestyle = "solid", alpha = 1):
    """
    Generates a plot of whatever you throw into it, but makes it easy to put things into the same axes.

    Parameters
    ----------
    dataframe : DATAFRAME
        DataFrame containing data you want to plot.
    ax : AXES, optional
        Axes you want to plot in. Default is None, which creates a new Axes.
    title : STRING, optional
        Title of figure. The default is "", which maintains title if already present.
    label : STRING, optional
        Label of plot. The default is "".
    xlabel : STRING, optional
        X-label of axes. The default is "Bin boundary [nm]".
    ylabel : STRING, optional
        Y-label of axes. The default is "[#particles cm${^-3}$]".

    Returns
    -------
    ax : AXES
        Axes that was plotted in.

    """
    if not ax:  # If no axis is given, make one
        fig, ax = plt.subplots(figsize=(10, 6))

    if not ax.get_title():  # If there is no title yet, set one.
        ax.set_title(title, fontsize = 15)
    else:
        if title:
            ax.set_title(title, fontsize = 15)  # If title is given, override current

    ax.set_xlabel(xlabel, fontsize = 15)
    ax.set_ylabel(ylabel, fontsize = 15)

    dataframe.plot(ax=ax, label=label, linestyle = linestyle, alpha = alpha)

    if label:
        ax.legend(fontsize = 14)

    ax.tick_params(axis = "both", labelsize = 12)    

    return ax


def read_measurement_data(data_name, show_comments=True):
    """
    Function that helps to read in measurement data. Just give (part of) the filename, and it'll figure out the rest

    Parameters
    ----------
    data_name : STR
        Name (or part of the name) of the measurement file.
    show_comments : BOOL, optional
        Print out at any metadata/comments found above the data. The default is True.

    Returns
    -------
    dataframe : DATAFRAME
        Pandas DataFrame containing the measuremnt data.

    """
    measurement_files = os.listdir(MEASUREMENTS_FOLDER)
    file = None
    comments = []

    for filename in measurement_files:
        if data_name.lower() in filename.lower():
            file = filename
            break

    if file is not None:
        # Dynamically determine comments and split them off - maybe return them for safekeeping metadata?
        # RIVM files have comments between double quotes instead of the usual #, so we're checking for both.
        file_path = os.path.join(MEASUREMENTS_FOLDER, file)

        try:
            skip_rows = 0
            with open(file_path, 'r', encoding = "utf-8") as f:
                for line in f:
                    if line.startswith("#") or line.startswith("\ufeff#"):
                        comment = line.strip().strip('"')
                        comments.append(comment)
                        skip_rows += 1
                    else:
                        break
                    
        except UnicodeDecodeError:
            skip_rows = 0
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith("#"):
                        comment = line.strip().strip('"')
                        comments.append(comment)
                        skip_rows += 1
                    else:
                        break
            
        try:
            dataframe = pd.read_csv(file_path, skiprows=skip_rows)
            comments = [comment for comment in comments if comment]  # Clean these up a little bit
            if len(dataframe.columns) == 1:
                # Let's just assume that this means the reader picked the wrong separator - switch to ;
                dataframe = pd.read_csv(file_path, skiprows = skip_rows, sep = ";")
            
        except UnicodeDecodeError:
            dataframe = pd.read_csv(file_path, skiprows=skip_rows, encoding = "ISO-8859-1", sep = ";")
            comments = [comment for comment in comments if comment]  # Clean these up a little bit

    else:
        print(f"No matching file found for {data_name}")
        dataframe = None

    if show_comments and comments != []:
        print(". ".join([comment.strip(".") for comment in comments]) + ".")
    
    if "datetime" in dataframe.columns:
        dataframe["datetime"] = pd.to_datetime(dataframe["datetime"])
    elif "Datetime Raw" in dataframe.columns:
        dataframe["datetime"] = pd.to_datetime(dataframe["Datetime Corr"])
        
    return dataframe


def plot_windrose(df, feature, title = "Figure Title", max_order = None, min_angle = 0, max_angle = 0):
    """
    Plots a windrose of the entered feature in the dataframe.

    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe containing info to be windrosed.
    feature : String
        Feature to be windrosed.
    title : String, optional
        Figure title. The default is "Figure Title".
    max_order : Float, optional33
        Order of magnitude to scale the windrose colours to. The default is None.

    Returns
    -------
    None.

    """
    ax = WindroseAxes.from_ax()
    
    feature_stats = df[feature].describe()
    feat_pct75 = feature_stats.iloc[6]
    feat_max = feature_stats.iloc[7]
    
    if order_of_magnitude(feat_pct75) == order_of_magnitude(feat_max):
        bins = np.linspace(0, np.floor(feat_max) - 2, 8)
        
    else:
        bins = np.linspace(0, 10**order_of_magnitude(feat_pct75) + 1, 8)
    
    ax.bar(df["winddir_deg"], df[feature], normed = True, bins =  bins, opening = 0.8,
           edgecolor = "white", cmap = co.cm.dense, alpha = 1)
    
    ax.set_legend(title = "m/s" if feature == "Wind Speed" else "# cm$^{-3}$", bbox_to_anchor=(-0.2, 0))
    
    # For highlighting certain bars
    # Convert direction range to radians
    direction_range = (min_angle, max_angle)
    direction_range_rad = np.deg2rad(direction_range)

    # Get the bin numbers corresponding to the direction range
    bin_range = (np.floor(direction_range[0] / 22.5), np.ceil(direction_range[1] / 22.5))

    for i, bar in enumerate(ax.patches):
        if not len(bins)*bin_range[0] <= i <= len(bins)*bin_range[1] + 8:
            bar.set_alpha(0.3)
    
    if title != "Figure Title":
        ax.set_title(title)
    else:
        ax.set_title(feature)
    
def stacked_timeseries_plot(df):
    """
    Generates vertically stacked timeseries plot, with each layer showing one bin of the particle distribution.

    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe containing the data to be shown (smps_df).

    Returns
    -------
    None.

    """
    channels = df.filter(regex="^[.0-9]+$|datetime")
    channels.set_index("datetime", inplace = True)
    
    fig, ax = plt.subplots(dpi = 80)
    colormap = plt.cm.rainbow_r
    
    colors = colormap(np.linspace(1, 0, len(channels.columns)))

    for ind, col in enumerate(channels.columns):
        ax.set_yscale("log")

        if ind == 0:
            ax.plot(channels.index, channels[col], color = colors[ind])
            ax.fill_between(channels.index, channels[col], 0, color = colors[ind], label = f"7 - {col} nm")

        elif ind < len(channels.columns) - 1:
            sum_so_far = channels[channels.columns[range(ind)]].sum(axis = 1)
            ax.plot(channels.index, channels[col] + sum_so_far, color = colors[ind])
            ax.fill_between(channels.index, channels[col] + sum_so_far, sum_so_far, color = colors[ind], label = f"{col} - {channels.columns[ind+1]} nm")

    ax.legend(title = "Bins", bbox_to_anchor = (1.0, 1.02))
    ax.set_ylim(bottom=1, top=1e6)
    ax.set_title("Timeseries of measurement")
    ax.set_xlabel("Time")
    ax.set_ylabel("# particles cm$^{-3}$")
    plt.xticks(rotation = 45)

    # num_ticks = 6
    # step = len(channels) // num_ticks
    # plt.xticks(df.index[::step], df["Datetime Corr"].dt.date.iloc[::step], rotation = 45)
    plt.show()

def show_bin_difference(smps_dataframe):
    """
    Generates schematic images of the difference in bin distribution and how to "resample" them.

    Parameters
    ----------
    smps_dataframe : Pandas DataFrame
        Contains the SMPS measurement series.

    Returns
    -------
    None.

    """
    # This part shows the difference pre-resampling
    ###############################################
    bin_boundaries = define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    
    salsa_boundaries = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[:8] * 1e9)
    smps_boundaries = np.concatenate([[7], smps_dataframe.filter(regex = "^[.0-9]+$").columns.astype(float).values])
    
    plt.figure(figsize = (10,6))
    plt.vlines(salsa_boundaries, 0, 1, label = "SALSA2.0")
    plt.vlines(smps_boundaries, 1.25, 2.25, color = "orange", label = "SMPS")
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.xlabel("Bin boundaries (nm)")
    plt.xscale('log')
    plt.legend(bbox_to_anchor = (1.17,1.018))
    plt.title("Schematic difference in bins SMPS / SALSA2.0")
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/bin_difference.jpg"), dpi = 75)
    plt.show()
    ###############################################
    
    plt.figure(figsize = (10,6))
    plt.vlines(salsa_boundaries[3:6], 0, 1, label = "SALSA2.0")
    plt.vlines(smps_boundaries[8:10], 1.25, 2.25, color = "orange", label = "SMPS")
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.xlabel("Bin boundaries (nm)")
    plt.xscale('log')
    plt.legend()
    plt.title("Bin definitions SMPS / SALSA2.0")
    plt.text(salsa_boundaries[3] - 1, -0.08, "$x_1$", fontsize = 14)
    plt.text(salsa_boundaries[4] - 1, -0.08, "$x_2$", fontsize = 14)
    plt.text(salsa_boundaries[5] - 1, -0.08, "$x_3$", fontsize = 14)
    plt.text(smps_boundaries[8] - 1, 1.25-0.08, "$y_1$", fontsize = 14)
    plt.text(smps_boundaries[9] - 1, 1.25-0.08, "$y_2$", fontsize = 14)
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/bin_definitions.jpg"), dpi = 75)
    plt.show()
    
    # This part shows schematically how to resample them
    ###############################################
    fig, ax = plt.subplots(figsize = (10,6))

    smps_ymin = 0
    smps_ymax = 1

    alpha = 1
    colormap = plt.cm.Blues
    colors = colormap(np.linspace(1, 0, len(salsa_boundaries[:-1])), alpha = alpha)
    colors_full = colormap(np.linspace(1, 0, len(salsa_boundaries[:-1])))

    line_color = complementary_color(*np.mean(colors_full, axis = 0)[:-1])
    vlines = ax.vlines(smps_boundaries, smps_ymin, smps_ymax, color = line_color, label = "SMPS", linewidth = 3)

    patches = []

    for ind, salsa_lower, salsa_upper, salsa_bin_name in zip(range(len(salsa_boundaries[:-1])), salsa_boundaries[:-1],
                                                             salsa_boundaries[1:], bin_names[:6]):
        # Fill between the boundaries
        ax.axvspan(xmin = salsa_lower, xmax = salsa_upper, ymin = 0.045, ymax = 0.955, color = colors[ind])
        patches.append(mpatches.Patch(color = colors[ind], label = salsa_bin_name))

    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlabel("Bin boundaries (nm)")
    ax.set_xscale('log')
    plt.title("Remapping of SMPS bins to SALSA2.0 bins")
    plt.legend(handles = [vlines] + patches, bbox_to_anchor = (1.13, 1.018))
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/bin_remapping.jpg"), dpi = 75)
    plt.show()
    ###############################################

def translate_particles(original_bins, original_counts, new_bins):
    """
    Translated particles from one system of bins to another

    Parameters
    ----------
    original_bins : List
        System of bins the data was given in.
    original_counts : List
        Number of particles per bin.
    new_bins : List
        System of bins to translate to.

    Returns
    -------
    output : List
        Number of particles in the new bin system.

    """
    columns = ['1a1', '1a2', '1a3', '2a1', '2a2', '2a3', '2a4', '2a5', '2a6', '2a7',
           '2b1', '2b2', '2b3', '2b4', '2b5', '2b6', '2b7']
    translated_counts =  np.interp(new_bins, original_bins, original_counts)
    output = pd.DataFrame(data = [translated_counts], columns = columns[:len(new_bins)])
    
    return output

def filter_outliers(df):
    """
    Filters outliers based on mean +- 2 sigma

    Parameters
    ----------
    df : DataFrame
        Needs to contain particle counts as given by the SMPS device.

    Returns
    -------
    df : DataFrame
        All outliers have been removed.

    """
    channels = df.filter(regex="^[.0-9]+$")
    
    for col in channels:
        df.loc[df[col] > np.mean(df[col]) + 2*np.std(df[col]), col] = np.nan
        df.loc[df[col] < np.mean(df[col]) - 2*np.std(df[col]), col] = np.nan
        
    return df

def smps_filter(df):
    """Quick filter based on the 3 common rules for SMPS data"""
    return df[(df["Status"] == "No errors") &
              (df["Total Conc"] >= 1000) &
              ~((df["datetime"].dt.day == 12) & (df["datetime"].dt.month == 7))]

def get_directional_dist(smps_dataframe, weather_dataframe, min_angle = 202.5, max_angle = 337.5):
    """
    Just a quick function that wraps things together, to make importing the relevant distribution easier
    
    Parameters
    ----------
    smps_dataframe : Pandas DataFrame
        Pandas DataFrame of the SMPS measurement series.
    weather_dataframe : Pandas DataFrame
        Pandas DataFrame of the DAVIS weather measurement series.
    
    Returns
    -------
    directional_dist : Pandas Series
        Series containing the mean directional distribution in #/m^-3.
    environmental_vars: List
        List containing the (ordered and) filtered environemntal variables.

    """    
    smps_bins = smps_dataframe.filter(regex="^[.0-9]+$").columns
    
    merged_df = pd.merge_asof(smps_dataframe, weather_dataframe, on = "datetime", direction = "nearest")
    
    try:
        merged_df = merged_df.drop(["time24"], axis = 1)
    except:
        pass
    
    merged_df = merged_df.dropna()
        
    directional_dist = merged_df[(merged_df["winddir_deg"] <= max_angle) &
                           (merged_df["winddir_deg"] >= min_angle) &
                           (merged_df["Status"] == "No errors")][smps_bins].mean()
    
    env = merged_df[(merged_df["winddir_deg"] <= max_angle) &
                    (merged_df["winddir_deg"] >= min_angle) &
                    (merged_df["Status"] == "No errors")][['Temp Out', 'Spec Hum','Bar']].mean()
    
    temp = env["Temp Out"] + 273.15
    hum = env["Spec Hum"]
    press = env["Bar"] * 100
    
    environmental_vars = [temp, hum, press]
    
    return directional_dist, environmental_vars
    
def show_haarrijn_data():
    """Quick visualisation of the measurement data at Haarrijn"""
    plt.figure(figsize = (10,6))
    
    plt.plot(smps_bins, smps_60m, label = "60m", color = "red")
    plt.scatter(smps_bins, smps_60m, color = "red")
    
    plt.plot(smps_bins, smps_200m, label = "200m", color = "green")
    plt.scatter(smps_bins, smps_200m, color = "green")
    
    plt.plot(smps_bins, smps_220m, label = "220m", color = "magenta")
    plt.scatter(smps_bins, smps_220m, color = "magenta")
    
    plt.plot(smps_bins, smps_320m, label = "320m", color = "blue")
    plt.scatter(smps_bins, smps_320m, color = "blue")
    
    
    plt.title("Particle distributions at varying distance")
    plt.xlabel("Diameter (nm)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.legend(title = "Distance")
    
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/Haarrijn_dists.jpg"), dpi = 150)
    
def show_haarrijn_init():
    """Quick visualisation of the 60m measurement data at Haarrijn, translated to SALSA2.0 bins"""
    salsa_bin_boundaries = define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[0:8] * 1e9)
    
    salsa_init = np.interp(salsa_bins, smps_bins, smps_60m)
    bar_widths = np.append(np.diff(salsa_bins), 0)
    bar_positions = salsa_bins + bar_widths / 2

    plt.figure(figsize = (10,6))
    plt.title("Translating measurements to model bins")

    plt.plot(smps_bins, smps_60m, label = "Measurement 60m")
    plt.bar(bar_positions, salsa_init, width = bar_widths, label = "SALSA2.0", color = "orange", alpha = 0.6, edgecolor = 'black')

    plt.xlabel("Diameter (nm)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.xlim((0, 400))

    plt.legend()
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/Haarrijn_init.jpg"), dpi = 150)


def show_methodology():
    show_bin_difference(rivm_201)
    show_model_lognormal_flux(sigma = 1.68, center_x = 90, scale = 1)
    show_haarrijn_data()
    show_haarrijn_init()


if __name__ == "__main__":
    davis_201 = read_measurement_data("Davis201", show_comments = False)
    davis_201["Spec Hum"] = utils.RH2q(davis_201["Out Hum"], davis_201["Bar"]*100, davis_201["Temp Out"] + 273.15)
    plot_windrose(davis_201, "Wind Speed", title = "Wind speed (N201)", min_angle = 202.5, max_angle = 270)
    
    davis_641 = read_measurement_data("Davis641", show_comments = False)
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
    plot_windrose(davis_641, "Wind Speed", title = "Wind speed (NL10641)", min_angle = 202.5, max_angle = 337.5)
    
    rivm_201 = read_measurement_data("SMPS", show_comments = False)
    # rivm_201 = filter_outliers(rivm_201)
    rivm_201 = smps_filter(rivm_201)
    stacked_timeseries_plot(rivm_201[rivm_201["Status"] == "No errors"])
    
    rivm_641 = read_measurement_data("RIVM", show_comments = False)
    # rivm_641 = filter_outliers(rivm_641)
    rivm_641 = smps_filter(rivm_641)
    stacked_timeseries_plot(rivm_641[rivm_641["Status"] == "No errors"])
    

    show_methodology()
    
    salsa_bin_boundaries = define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[0:8] * 1e9)
    salsa_bins_float = [float(binbound) for binbound in salsa_bins]
    salsa_upper_boundaries = np.unique(np.concatenate(list(salsa_bin_boundaries.values()), 0)[1:8] * 1e9)
    smps_bins = rivm_201.filter(regex="^[.0-9]+$").columns
    smps_bins_float = [3.] + [float(binbound) for binbound in smps_bins]
    
    merged_201_df = pd.merge_asof(rivm_201, davis_201, on = "datetime", direction = "nearest")
    merged_201_df = merged_201_df.drop(["time24"], axis = 1)
    merged_201_df = merged_201_df.dropna()
    
    merged_641_df = pd.merge_asof(rivm_641, davis_641, on = "datetime", direction = "nearest")
    merged_641_df = merged_641_df.dropna()
    
    # Calculate directional distributions
    highway_dist_201, _ = get_directional_dist(rivm_201, davis_201, min_angle = 202.5, max_angle = 270)
    background_dist_201, _ = get_directional_dist(rivm_201, davis_201, min_angle = 90, max_angle = 180)
    
    highway_dist_641, _ = get_directional_dist(rivm_641, davis_641, min_angle = 202.5, max_angle = 337.5)
    background_dist_641, _ = get_directional_dist(rivm_641, davis_641, min_angle = 45, max_angle = 135)
    
    highway_201_resampled = translate_particles(smps_bins_float[1:], highway_dist_201, salsa_bins_float[1:])
    highway_641_resampled = translate_particles(smps_bins_float[1:], highway_dist_641, salsa_bins_float[1:])
    

    # Show the various mean distribution of particles
    
    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    plot_data(highway_dist_201, ax = axd["main"], label = "W - SSW", title="Mean size distribution (N201)")
    plot_data(background_dist_201, ax = axd["main"], label = "E-S")
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/N201_Dist.jpg"), bbox_inches = "tight", dpi = 150)
    plt.show()
    
    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    plot_data(highway_dist_641, ax = axd["main"], label = "NNW - SSW", title="Mean size distribution (NL10641)")
    plot_data(background_dist_641, ax = axd["main"], label = "NE-SE")
    plt.savefig(os.path.join(RESULTS_FOLDER, "Measurement Figures/N641_Dist.jpg"), bbox_inches = "tight", dpi = 150)
    plt.show()
    
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
