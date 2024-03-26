# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 08:20:04 2024

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.animation import FuncAnimation

import seaborn as sb
import cmocean as co
from windrose import WindroseAxes

from utils import order_of_magnitude, complementary
from HAM_plot import define_bin_boundaries

DATA_FOLDER = "../data"
MEASUREMENTS_FOLDER = os.path.join(DATA_FOLDER, "Measurements")
RESULTS_FOLDER = "../results"


def plot_data(dataframe, ax=None, title="", label="",
              xlabel="Bin boundary [nm]", ylabel="[#particles cm${^-3}$]",
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
        ax.set_title(title)
    else:
        if title:
            ax.set_title(title)  # If title is given, override current

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    dataframe.plot(ax=ax, label=label, linestyle = linestyle, alpha = alpha)

    if label:
        ax.legend()

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

        skip_rows = 0

        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('"') or line.startswith("#"):
                    comment = line.strip().strip('"')
                    comments.append(comment)
                    skip_rows += 1
                else:
                    break

        dataframe = pd.read_csv(file_path, skiprows=skip_rows)
        comments = [comment for comment in comments if comment]  # Clean these up a little bit

    else:
        print(f"No matching file found for {data_name}")
        dataframe = None

    if show_comments and comments != []:
        print(". ".join([comment.strip(".") for comment in comments]) + ".")

    return dataframe


def plot_windrose(df, feature, title = "Figure Title", max_order = None):
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
    max_order : Float, optional
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
        bins = np.linspace(0, np.floor(feat_max) - 1, 8)
        
    else:
        bins = np.linspace(0, 10**order_of_magnitude(feat_pct75) + 1, 8)
    
    ax.bar(df["winddir_deg"], df[feature], normed = True, bins =  bins, opening = 0.8, edgecolor = "white", cmap = co.cm.dense)
    ax.set_title(title)
    ax.set_legend(title = "m/s")

def create_animated_plot(df, window_size=100, fps=24, save_path='animation.gif'):
    fig, ax = plt.subplots()

    vmin = df.min().min()
    vmax = df.max().max()

    def update(frame):
        ax.clear()
        start_index = frame
        end_index = start_index + window_size
        data_slice = df.iloc[start_index:end_index, :]
        im = ax.imshow(data_slice.values.T, cmap='plasma', aspect='auto',
                       vmin=vmin, vmax=vmax, interpolation="none")
        plt.gca().invert_yaxis()
        ax.set_title(f'Frame: {frame}')

        ax.set_yticks(range(len(data_slice.columns)))
        ax.set_yticklabels(data_slice.columns)
        ax.set_ylabel("Bin boundary [nm]")

        return im

    cbar = plt.colorbar(update(0), ax=ax)
    cbar.set_label("[#particles cm${^-3}$]")

    num_frames = len(df) - window_size + 1
    animation = FuncAnimation(fig, update, frames=num_frames, repeat=False)

    # Save the animation
    animation.save(os.path.join(RESULTS_FOLDER, save_path), fps=fps, writer='pillow')

def stacked_timeseries_plot(df):
    """
    Generates vertically stacked timeseries plot, with each layer showing one bin of the particle distribution.

    Parameters
    ----------
    df : Pandas DataFrame
        Dataframe containing the data to be shown (merged_df).

    Returns
    -------
    None.

    """
    fig, ax = plt.subplots()
    colormap = plt.cm.rainbow_r
    
    colors = colormap(np.linspace(1, 0, len(df.columns)))

    for ind, col in enumerate(df.columns):
        ax.set_yscale("log")

        if ind == 0:
            ax.plot(df.index, df[col], color = colors[ind])
            ax.fill_between(df.index, df[col], 0, color = colors[ind], label = f"7 - {col} nm")

        elif ind < len(df.columns) - 1:
            sum_so_far = df[df.columns[range(ind)]].sum(axis = 1)
            ax.plot(df.index, df[col] + sum_so_far, color = colors[ind])
            ax.fill_between(df.index, df[col] + sum_so_far, sum_so_far, color = colors[ind], label = f"{col} - {df.columns[ind+1]} nm")

    ax.legend(title = "Bins", bbox_to_anchor = (1.0, 1.02))
    ax.set_ylim(bottom=1, top=1e5)
    ax.set_title("Timeseries of measurement")
    ax.set_xlabel("Time")
    ax.set_ylabel("# particles cm$^{-3}$")

    num_ticks = 6
    step = len(df) // num_ticks
    plt.xticks(df.index[::step], df["Datetime Corr"].dt.date.iloc[::step], rotation = 45)
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
    
    salsa_boundaries = np.concatenate(list(bin_boundaries.values()), 0)[:8] * 1e9
    smps_boundaries = smps_dataframe.filter(regex = "^[.0-9]+$").columns.astype(float).values
    
    plt.figure(figsize = (10,6))
    plt.vlines(salsa_boundaries, 0, 1, label = "SALSA2.0")
    plt.vlines(smps_boundaries, 1.25, 2.25, color = "orange", label = "SMPS")
    plt.gca().axes.get_yaxis().set_visible(False)
    plt.xlabel("Bin boundaries (nm)")
    plt.xscale('log')
    plt.legend(bbox_to_anchor = (1.17,1.018))
    plt.title("Schematic difference in bins SMPS / SALSA2.0")
    plt.show()
    ###############################################
    
    # This part shows how to resample them
    ###############################################
    fig, ax = plt.subplots(figsize = (10,6))

    smps_ymin = 0
    smps_ymax = 1

    alpha = 1
    colormap = plt.cm.Blues
    colors = colormap(np.linspace(1, 0, len(salsa_boundaries[:-1])), alpha = alpha)
    colors_full = colormap(np.linspace(1, 0, len(salsa_boundaries[:-1])))

    line_color = complementary(*np.mean(colors_full, axis = 0)[:-1])
    vlines = ax.vlines(smps_boundaries, smps_ymin, smps_ymax, color = line_color, label = "SMPS", linewidth = 3)

    patches = []

    for ind, salsa_lower, salsa_upper, salsa_bin_name in zip(range(len(salsa_boundaries[:-1])), salsa_boundaries[:-1],
                                                             salsa_boundaries[1:], bin_names[:6]):
        
        ax.axvspan(xmin = salsa_lower, xmax = salsa_upper, ymin = 0.045, ymax = 0.955, color = colors[ind])
        patches.append(mpatches.Patch(color = colors[ind], label = salsa_bin_name))

    ax.axes.get_yaxis().set_visible(False)
    ax.set_xlabel("Bin boundaries (nm)")
    ax.set_xscale('log')
    plt.title("Remapping of SMPS bins to SALSA2.0 bins")
    plt.legend(handles = [vlines] + patches, bbox_to_anchor = (1.13, 1.018))
    plt.show()
    ###############################################

def resample_to_salsa(series):
    """
    Resamples the SMPS dataseries (mean) into the bins of SALSA2.0.

    Parameters
    ----------
    series : Pandas Series
        Pandas Series of the mean NW-SW distribution.

    Returns
    -------
    resampled : Pandas Series
        Pandas Series of the resampled mean NW-SW distribution..

    """
    bin_boundaries = define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    
    salsa_boundaries = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[:8]  * 1e9)
    smps_boundaries = np.append(np.array([7.0]), series.index.astype(float).values)
    new_data = {}

    for salsa_lower, salsa_upper, salsa_bin_names in zip(salsa_boundaries[:-1], salsa_boundaries[1:], bin_names[:6]):
        for smps_lower, smps_upper, smps_num in zip(smps_boundaries[:-1], smps_boundaries[1:], series[:-1]):

            # Fully contained within salsa bounds
            if (smps_lower > salsa_lower) and (smps_upper < salsa_upper):
                new_data[salsa_upper] = new_data.get(salsa_upper, 0) + smps_num

            # Partially contained above
            elif (smps_lower > salsa_lower) and (smps_upper > salsa_upper) and not (smps_lower > salsa_upper):
                weight = (salsa_upper - smps_lower) / (smps_upper - smps_lower)
                new_data[salsa_upper] = new_data.get(salsa_upper, 0) + smps_num*weight

            # Partially contained below
            elif (smps_lower < salsa_lower) and (smps_upper < salsa_upper) and not (smps_upper < salsa_lower):
                weight = (smps_upper - salsa_lower) / (smps_upper - smps_lower)
                new_data[salsa_upper] = new_data.get(salsa_upper, 0) + smps_num*weight

            else:
                new_data[salsa_upper] = new_data.get(salsa_upper, 0)
                
    resampled = pd.Series(new_data)
    return resampled

if __name__ == "__main__":
    weather_df = read_measurement_data("Davis")
    weather_df["datetime"] = pd.to_datetime(weather_df["datetime"])
    smps_df = read_measurement_data("SMPS")
    smps_df["datetime"] = pd.to_datetime(smps_df["Datetime Raw"])
    
    show_bin_difference(smps_df)

    channels = smps_df.filter(regex="^[.0-9]+$")
    merged_df = pd.merge_asof(smps_df, weather_df, on = "datetime", direction = "nearest")
    
    merged_df["winddir_deg"] = (merged_df["winddir_deg"] + 180 ) % 360 
    nw_sw_dist = merged_df[(merged_df["winddir_deg"] <= 315) & (merged_df["winddir_deg"] >= 225)][channels.columns].mean()
    nw_sw_resampled = resample_to_salsa(nw_sw_dist)

    # Show the mean distribution of particles
    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    plot_data(nw_sw_dist, ax = axd["main"], label = "SW - NW")
    plot_data(channels.mean(), ax=axd["main"], label="All", title="Mean size distribution", linestyle = "dashed", alpha = 0.6)
    plot_data(channels[smps_df["Status"] != "No errors"].mean(), ax=axd["main"], label="No errors", linestyle = "dashed", alpha = 0.6)
    plot_data(channels[smps_df["Status"] == "No errors"].mean(), ax=axd["main"], label="Only error", linestyle = "dashed", alpha = 0.6)
    plt.show()
    
    # Show the resampling of the nw_sw distribution
    plt.plot(nw_sw_resampled, label = "SALSA2.0")
    plt.plot(nw_sw_dist.index.astype(float), nw_sw_dist.values, label = "SMPS")
    plt.legend()
    plt.xlabel("Upper bin boundary (nm)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.title("Resampled distribution")
    plt.show()
    
# =============================================================================
#     plt.figure(figsize = (10,6))
#     plt.title("Correlations")
#     corr = merged_df.filter(regex = "(?i)^[.0-9]+$|Wind|Temp|Hum|Rain").drop(["Wind Dir", "Wind Tx"], axis = 1).corr()
#     cmap = sb.diverging_palette(5, 250, as_cmap = True)
#     sb.heatmap(corr, cmap = "Blues", mask = np.triu(corr), center = 0)
# =============================================================================

    # create_animated_plot(channels.iloc[0:500], fps = 50, window_size=100, save_path='SMPS_animation.gif')
