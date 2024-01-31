# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 08:20:04 2024

@author: rens_
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
from matplotlib.animation import FuncAnimation

DATA_FOLDER = "../data"
MEASUREMENTS_FOLDER = os.path.join(DATA_FOLDER, "Measurements")
RESULTS_FOLDER = "../results"


def plot_data(dataframe, ax=None, title="", label="", xlabel="Bin boundary [nm]", ylabel="[#particles cm${^-3}$]"):
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

    dataframe.plot(ax=ax, label=label)

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


if __name__ == "__main__":
    weather_df = read_measurement_data("Davis")
    weather_df["datetime"] = pd.to_datetime(weather_df["datetime"])
    smps_df = read_measurement_data("SMPS")
    smps_df["datetime"] = pd.to_datetime(smps_df["Datetime Raw"])

    channels = smps_df.filter(regex="^[.0-9]+$")
    merged_df = pd.merge_asof(smps_df, weather_df, on = "datetime", direction = "nearest")

    fig = plt.figure(layout="tight", figsize=(10, 6))
    axd = fig.subplot_mosaic([["main"]])
    plot_data(channels.mean(), ax=axd["main"], label="All", title="Mean size distribution")
    plot_data(channels[smps_df["Status"] != "No errors"].mean(), ax=axd["main"], label="No errors")
    plot_data(channels[smps_df["Status"] == "No errors"].mean(), ax=axd["main"], label="Only error")
    plt.show()

    plt.figure(figsize = (10,6))
    plt.title("Correlations")
    corr = merged_df.filter(regex = "(?i)^[.0-9]+$|Wind|Temp|Hum|Rain").drop(["Wind Dir", "Wind Tx"], axis = 1).corr()
    cmap = sb.diverging_palette(5, 250, as_cmap = True)
    sb.heatmap(corr, cmap = "Blues", mask = np.triu(corr), center = 0)

    # create_animated_plot(channels.iloc[0:500], fps = 50, window_size=100, save_path='SMPS_animation.gif')
