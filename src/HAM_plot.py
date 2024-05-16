# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:28:14 2024

@author: rens_
"""
import os
import re

import cmocean as co
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mc
from scipy.interpolate import griddata

from utils import order_of_magnitude

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")

def plot_size_dist(
    rdry, num, rows=[0], populations=['a', 'b'],
    xmin=None, xmax=None,
    ymin=None, ymax=None,
    fig=None, axes=None,
    exp_name="", title="", name_addition="",
    **kwargs
):
    """
    Plots size distributions at various timesteps.

    Parameters
    ----------
    rdry : DataFrame
        Contains the bin boundaries.
    num : DataFrame
        Contains the numbers of particles per bin.
    rows : List, optional
        Time steps to be read and plotted. The default is [0].
    populations : List, optional
        Which populations of particles to show. The default is ['a', 'b'].
    xmin : Float, optional
        Left x-axis limit. The default is None.
    xmax : Float, optional
        Right x-axis limit. The default is None.
    ymin : Float, optional
        Bottom y-axis limit. The default is None.
    ymax : Float, optional
        Top y-axis limit. The default is None.
    exp_name : String, optional
        Name of the experiment. The default is "". Leave empty to forego saving the image.
    title : String, optional
        Title of the image. The default is "".
    name_addition: String. The default is ""
        Suffix for filename.
    Returns
    -------
    None.

    """

    # make sure that row_nr is a list-like object
    try:
        iter(rows)
    except Exception:
        rows = [rows]

    nrows = 1
    ncols = len(populations)

    if fig == None and axes == None:
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(9*ncols, 6*nrows),
            sharex=True,
            sharey=True,
            dpi = 75
        )
        
        plt.tight_layout()
        fig.text(-0.01, 0.5, "# particles cm$^{-3}$", va="center", rotation="vertical")
        fig.text(0.47, 0.04, "Diameter (m)", ha="center")

    if ncols > 1:
        label = kwargs.get("label", "")
        if label != "":
            kwargs.pop("label")
            
        for pop, ax in zip(populations, axes.ravel()):
            bins = [col for col in rdry.columns if pop in col]
            bins = [binName for binName in bins if binName in num.columns]
    
            for n in rows:
    
                r_row = rdry.iloc[n][bins]
                N_row = num.iloc[n][bins]

                ax.plot(r_row, N_row, label=f"{label + ' '}{n} s" if label == "" else f"{label}", **kwargs)
                
            ax.set_title(f"Population {pop}")
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)
            if pop == "b":
                ax.legend(title="Time", bbox_to_anchor=(1.22, 1.02))
                
            if title != "":
                fig.suptitle(title, y = 1, x = 0.47)
                
    else:
        label = kwargs.get("label", "placeholder")
        if label != "placeholder":
            kwargs.pop("label")
            
        for pop in populations:
            ax = axes
            bins = [col for col in rdry.columns if pop in col]
            bins = [binName for binName in bins if binName in num.columns]

            for n in rows:
    
                r_row = rdry.iloc[n][bins]
                N_row = num.iloc[n][bins]
    
                ax.plot(r_row, N_row, label=f"{label + ' '}{n} s" if label == "Model" else f"{label}", **kwargs)
                
            ax.set_title(f"Population {pop}")
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)

            ax.legend(title="Data", bbox_to_anchor=(1.01, 1.02))
                
            if title != "":
                fig.suptitle(title, y = 1, x = 0.52)

    if exp_name != "":
        figure_name = f"size_distribution_{name_addition}.png" if name_addition != "" else "size_distribution.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
        
    return fig, axes

def define_bin_boundaries(populations=['1a', '2a', '2b']):
    """
    Defines bin boundaries for the various particle populations.

    Parameters
    ----------
    populations : List, optional
        List of particle populations. The default is ['1a', '2a', '2b'].

    Returns
    -------
    dict
        Bin boundaries per population.

    """
    return {
        pop: (
            np.logspace(np.log10(3e-9), np.log10(50e-9), 4) if pop[0] == '1' else
            np.concatenate([
                np.logspace(np.log10(50e-9), np.log10(700e-9), 5, endpoint=True)[:-1],
                np.logspace(np.log10(700e-9), np.log10(10000e-9), 4, endpoint=True),
            ])
        )
        for pop in populations
    }


def plot_size_dist_evolution(
        rdry, num, populations=['a', 'b'],
        xmin=None, xmax=None,
        ymin=None, ymax=None,
        vmin=None, vmax=None,
        exp_name="", title="", name_addition=""):
    """


    Parameters
    ----------
    rdry : DataFrame
        Contains the bin boundaries.
    num : DataFrame
        Contains the numbers of particles per bin.
    rows : List, optional
        Time steps to be read and plotted. The default is [0].
    populations : List, optional
        Which populations of particles to show. The default is ['a', 'b'].
    xmin : Float, optional
        Left x-axis limit. The default is None.
    xmax : Float, optional
        Right x-axis limit. The default is None.
    ymin : Float, optional
        Bottom y-axis limit. The default is None.
    ymax : Float, optional
        Top y-axis likmit. The default is None.
    vmin : Float, optional
        Bottom limit of the colorbar. The default is None.
    vmax : FLoat, optional
        Top limit of the colorbar. The default is None.
    exp_name : String, optional
        Name of the experiment. The default is "". Leave empty to forego saving the image.
    title : String, optional
        Title of the image. The default is "".
    name_addition: String. The default is ""
        Suffix for filename.

    Returns
    -------
    None.

    """

    bin_boundaries = define_bin_boundaries()

    nrows = len(populations)
    ncols = 1

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(12*ncols, 4*nrows),
        sharex=True,
        sharey=True,
    )

    # computing common color bar bounds:
    vmin = num.values.min() if vmin is None else vmin
    vmax = num.values.max() if vmax is None else vmax

    if nrows > 1:
        for pop, ax in zip(populations, axes.ravel()):
            bins = [col for col in rdry.columns if pop in col]
            bins = [binName for binName in bins if binName in num.columns]
    
            # combining bin boundaries (omitting double values)
            rbounds = [bounds for k, bounds in sorted(bin_boundaries.items()) if pop in k]
            rbounds = np.concatenate([rb[:-1] for rb in rbounds[:-1]]+[rbounds[-1]])
    
            # time bounds
            tbounds = np.arange(num.shape[0]+1)
    
            # generating meshgrid
            t, r = np.meshgrid(tbounds, rbounds)
    
            norm = mc.LogNorm(vmin=vmin, vmax=vmax)
            cls = ax.pcolormesh(t, r, num[bins].T, norm=norm, cmap=co.cm.dense)
    
            ax.set_title(f"Population {pop}")
            ax.set_yscale('log')
            ax.set_ylabel("Diameter (m)")
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)
    else:
        for pop in populations:
            ax = axes
            bins = [col for col in rdry.columns if pop in col]
            bins = [binName for binName in bins if binName in num.columns]
    
            # combining bin boundaries (omitting double values)
            rbounds = [bounds for k, bounds in sorted(bin_boundaries.items()) if pop in k]
            rbounds = np.concatenate([rb[:-1] for rb in rbounds[:-1]]+[rbounds[-1]])
    
            # time bounds
            tbounds = np.arange(num.shape[0]+1)
    
            # generating meshgrid
            t, r = np.meshgrid(tbounds, rbounds)
    
            norm = mc.LogNorm(vmin=vmin, vmax=vmax)
            cls = ax.pcolormesh(t, r, num[bins].T, norm=norm, cmap=co.cm.dense)
    
            ax.set_title(f"Population {pop}")
            ax.set_yscale('log')
            ax.set_ylabel("Diameter (m)")
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)

    if title != "":
        fig.suptitle(title, y = 1, x = 0.435)

    if nrows > 1:
        cbar = fig.colorbar(cls, ax=axes.ravel().tolist())
    else:
        cbar = fig.colorbar(cls)
        
    cbar.set_label("# particles m$^{-3}$")
    ax.set_xlabel("Time (s)")

    if exp_name != "":
        figure_name = f"size_distribution_LES_box_{name_addition}.png" if name_addition != "" else "size_distribution_LES_box.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
    # plt.close()
    

def stacked_timeseries_plot(
        num, populations = ['a', 'b'],
        xmin = None, xmax = None,
        ymin = None, ymax = None,
        exp_name = "", title = "", name_addition = "",
        colormap = None):
    """
    Generates vertically stacked timeseries, each layer being a bin of the model.

    Parameters
    ----------
    num : DataFrame
        Contains the numbers of particles per bin.
    populations : List, optional
        Which populations of particles to show. The default is ['a', 'b'].
    xmin : Float, optional
        Left x-axis limit. The default is None.
    xmax : Float, optional
        Right x-axis limit. The default is None.
    ymin : Float, optional
        Bottom y-axis limit. The default is None.
    ymax : Float, optional
        Top y-axis limit. The default is None.
    exp_name : String, optional
        Name of the experiment. The default is "". Leave empty to forego saving the image.
    title : String, optional
        Title of the image. The default is "".
    name_addition: String. 
        Suffix for filename. The default is "".
    Returns
    -------
    None.

    """
    
    nrows = len(populations)
    ncols = 1

    fig, axes = plt.subplots(
        nrows = nrows,
        ncols = ncols,
        figsize = (12*ncols, 4*nrows),
        sharex = True, 
        sharey = True,
        tight_layout = True)

    
    bin_boundaries = define_bin_boundaries()

    if colormap == None:
        colormap = co.cm.dense
    
    
    if nrows > 1:
        for pop, ax in zip(populations, axes.ravel()):
            bins = [col for col in num.columns if pop in col]
            colors = colormap(np.linspace(0, 1, len(bins)))
            
            for ind, col in enumerate(num[bins].columns):
                col_split = re.split("(?<=[a-zA-Z])(?=\d)", col)
                lower_boundary = bin_boundaries[col_split[0]][int(col_split[1])-1]*1e9
                upper_boundary = bin_boundaries[col_split[0]][int(col_split[1])]*1e9
                # Add >6 before the periods to have them left-align. Maybe looks better?
                label = f"{lower_boundary:.2f} - {upper_boundary:.2f} nm" if upper_boundary < 1000 \
                        else f"{lower_boundary*1e-3:.2f} - {upper_boundary*1e-3:.2f} $\\rm \mu$m"
                # ax.set_yscale("log")
                
                if ind == 0:
                    ax.plot(num.index, num[col], color = colors[ind])
                    ax.fill_between(num.index, num[col], 0, color = colors[ind], label = label)
    
                else:
                    sum_so_far = num[num[bins].columns[range(ind)]].sum(axis = 1)
                    ax.plot(num.index, num[col] + sum_so_far, color = colors[ind])
                    ax.fill_between(num.index, num[col] + sum_so_far, sum_so_far, color = colors[ind], label = label)
    
            ax.legend(title = "Bins", bbox_to_anchor = (1.0, 1.02))
            ax.set_title(f"Population {pop}")
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)
            
    else:
        for pop in populations:
            ax = axes
            bins = [col for col in num.columns if pop in col]
            colors = colormap(np.linspace(0, 1, len(bins)))
            
            for ind, col in enumerate(num[bins].columns):
                col_split = re.split("(?<=[a-zA-Z])(?=\d)", col)
                lower_boundary = bin_boundaries[col_split[0]][int(col_split[1])-1]*1e9
                upper_boundary = bin_boundaries[col_split[0]][int(col_split[1])]*1e9
                # Add >6 before the periods to have them left-align. Maybe looks better?
                label = f"{lower_boundary:.2f} - {upper_boundary:.2f} nm" if upper_boundary < 1000 \
                        else f"{lower_boundary*1e-3:.2f} - {upper_boundary*1e-3:.2f} $\\rm \mu$m"
                # ax.set_yscale("log")
                
                if ind == 0:
                    ax.plot(num.index, num[col], color = colors[ind])
                    ax.fill_between(num.index, num[col], 0, color = colors[ind], label = label)
    
                else:
                    sum_so_far = num[num[bins].columns[range(ind)]].sum(axis = 1)
                    ax.plot(num.index, num[col] + sum_so_far, color = colors[ind])
                    ax.fill_between(num.index, num[col] + sum_so_far, sum_so_far, color = colors[ind], label = label)
    
            ax.legend(title = "Bins", bbox_to_anchor = (1.0, 1.02))
            ax.set_title(f"Population {pop}")
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)

    if title != "":
        fig.suptitle(title, y = 1, x = 0.435)

    plt.tight_layout()
    # y-label
    fontsize = 13
    fig.text(-0.01, 0.5, "# particles m$^{-3}$", va = "center", rotation = "vertical", fontsize = fontsize)
    ax.set_xlabel("Time (s)", fontsize = fontsize)
    ax.legend(title = "Bins", bbox_to_anchor = (1.0, 1.02))  

    if exp_name != "":
        figure_name = f"stacked_distribution_timeseries_{name_addition}.png" if name_addition != "" else "stacked_distribution_timeseries.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)

    

def plot_variation_scatter(numdict, metadict, exp_name = "", binName = "1a1", time = 0, name_addition = ""):
    """
    Generates sensitivity scatters of a environmental variable variation series.

    Parameters
    ----------
    numdict : Dictionary
        Dictionary of Pandas Dataframes (like num in other functions).
    metadict : Dictionary
        Dictinary of metadata strings.
    exp_name : String, optional
        Name of the experiment. Leave empty to forego saving the image. The default is "".
    binName : String, optional
        Name of the bin to be plotted. The default is "1a1".
    time : Integer, optional
        Time to be plotted. The default is 0.
    name_addition : String, optional
        File name addition, for easy lookup. The default is "".

    Returns
    -------
    None.

    """
    colormap = plt.cm.viridis
    fig = plt.figure(figsize = (10,6))
    plt.tight_layout()
    ax = fig.add_subplot(projection = "3d")

    numlist = [num[binName].iloc[time] for num in numdict.values()]
    norm = mc.Normalize(vmin = min(numlist), vmax = max(numlist))
    for key, num in numdict.items():
        metadata = metadict[key]
        environmental_vals = [float(line.split("=")[1].strip()) for line in metadata.split("\n") if line.startswith(("pt", "pqm1", "pap"))]

        pt, pqm1, pap = environmental_vals

        sc = ax.scatter(pt, pqm1, pap, c = num[binName].iloc[time], cmap = colormap, norm = norm)

    cbar = plt.colorbar(sc)
    cbar.set_label("# particles m$^{-3}$")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel("q $(\\frac{kg}{kg})$")
    ax.set_zlabel("Air pressure (Pa)")
    # ax.view_init(30, 90)
    plt.title(f"{binName} at {time} s")

    figure_name = f"sensitivity_scatter_{name_addition}.png" if name_addition != "" else "sensitivity_scatter.png"
    savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

    if not os.path.exists(savepath):
        os.makedirs(savepath)

    full_savepath = os.path.join(savepath, figure_name)
    plt.tight_layout()
    
    if exp_name != "":
        figure_name = f"Sensitivity_scatter_{name_addition}.png" if name_addition != "" else "sensitivity_scatter.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)


def plot_variation_surface(numdict, metadict, exp_name = "", binName = "1a1", time = 0, elev = 20, azi = 110, name_addition = "",
                           fig = None, ax = None, colormap = None):
    """
    Generates sensitivity surfaces of environmental variable variation runs.

    Parameters
    ----------
   numdict : Dictionary
       Dictionary of Pandas Dataframes (like num in other functions).
   metadict : Dictionary
       Dictinary of metadata strings.
   exp_name : String, optional
       Name of the experiment. Leave empty to forego saving the image. The default is "".
   binName : String, optional
       Name of the bin to be plotted. Entering multitple sums the bins together. The default is "1a1".
   time : Integer, optional
       Time to be plotted. The default is 0.
    elev : Float, optional
        Elevation of the perspective. The default is 20.
    azi : Float, optional
        Azimuth of the perspective. The default is 110.
    name_addition : String, optional
        File name addition, for easy lookup. The default is "".
    fig : Matplotlib Figure, optional
        Figure to add multiple surfaces in one. The default is None.
    ax : Matplotlib Axes, optional
        Axis to add multiple surfaces in one. The default is None.
    colormap : Matplotlib Colormap, optional
        Colormap to color surfaces. The default is None.

    Returns
    -------
    fig : Matplotlib Figure, optional
        Figure to add multiple surfaces in one.
    ax : Matplotlib Axes, optional
        Axis to add multiple surfaces in one.

    """
    if fig is None and ax is None:
        fig = plt.figure(figsize = (10,6))
        ax = fig.add_subplot(projection = "3d")

    numlist = []
    environmental_vals = []
    
    if type(binName) == str:
        for key, num in numdict.items():
            metadata = metadict[key]
            environmentals = [float(line.split("=")[1].strip()) for line in metadata.split("\n") if line.startswith(("pt", "pqm1", "pap"))]
    
            numlist.append(num[binName].iloc[time])
            environmental_vals.append(environmentals)
      
    # If we input a list of bins, assume we want the total (sum) of it
    elif type(binName) == list:
        for key, num in numdict.items():
            metadata = metadict[key]
            environmentals = [float(line.split("=")[1].strip()) for line in metadata.split("\n") if line.startswith(("pt", "pqm1", "pap"))]
    
            numlist.append(num[binName].sum(axis = 1).iloc[time])
            environmental_vals.append(environmentals)

    environmental_df = pd.DataFrame(data = environmental_vals)

    gridpoints = 1000
    temp_grid = np.linspace(environmental_df[0].min(), environmental_df[0].max(), num = gridpoints)
    humi_grid = np.linspace(environmental_df[1].min(), environmental_df[1].max(), num = gridpoints)
    pap_grid = np.linspace(environmental_df[2].min(), environmental_df[2].max(), num = gridpoints)

    constant_var = [column for column in environmental_df.columns if environmental_df[column].nunique() == 1][0]

    if constant_var == 0:
        x_grid, y_grid = np.meshgrid(humi_grid, pap_grid)
        num_grid = griddata((environmental_df[1].values, environmental_df[2].values), numlist, (x_grid, y_grid), method = "linear")
        x_grid = x_grid*1e3
        y_grid = y_grid*1e-2

        ax.set_xlabel("q $(\\frac{kg}{kg})$ $\\times 10^{-3}$")
        ax.set_ylabel("Ambient pressure (hPa)")
        ax.yaxis.labelpad = 5
        constant_title = f"T = {environmental_df[0].iloc[0]:.2f} K"

    elif constant_var == 1:
        x_grid, y_grid = np.meshgrid(temp_grid, pap_grid)
        num_grid = griddata((environmental_df[0].values, environmental_df[2].values), numlist, (x_grid, y_grid), method = "linear")
        y_grid = y_grid*1e-2
        
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Ambient pressure (hPa)")
        constant_title = f"q = {environmental_df[1].iloc[0]:.3e}" + " $(\\frac{kg}{kg})$"

    elif constant_var == 2:
        x_grid, y_grid = np.meshgrid(temp_grid, humi_grid)
        num_grid = griddata((environmental_df[0].values, environmental_df[1].values), numlist, (x_grid, y_grid), method = "linear")
        y_grid = y_grid*1e3

        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("q $(\\frac{kg}{kg})$ $\\times 10^{-3}$")
        constant_title = f"p = {environmental_df[2].iloc[0]:.3e} Pa"

    max_order = order_of_magnitude(np.abs(num_grid).max())

    num_grid = num_grid * (10 ** -max_order)

    
    if colormap == None:
        colormap = co.cm.dense_r
        
    norm = mc.Normalize(vmin = num_grid.min(), vmax = num_grid.max())
    bin_boundaries = define_bin_boundaries()
    
    if type(binName) == str:
        binName_split = re.split("(?<=[a-zA-Z])(?=\d)", binName)
        lower_boundary = bin_boundaries[binName_split[0]][int(binName_split[1])-1]
        upper_boundary = bin_boundaries[binName_split[0]][int(binName_split[1])]
        title = f"Bin {lower_boundary*1e9:.2f} - {upper_boundary*1e9:.2f} nm at {time} s \n[{constant_title}]"
    elif type(binName) == list:
        lowerbinName_split = re.split("(?<=[a-zA-Z])(?=\d)", binName[0])
        upperbinName_split = re.split("(?<=[a-zA-Z])(?=\d)", binName[-1])
        lower_boundary = bin_boundaries[lowerbinName_split[0]][int(lowerbinName_split[1])-1]
        upper_boundary = bin_boundaries[upperbinName_split[0]][int(upperbinName_split[1])]
        title = f"Sum of bins {lower_boundary*1e9:.2f} - {upper_boundary*1e9:.2f} nm at {time} s \n[{constant_title}]"
    
    fig.suptitle(title)
    ax.plot_surface(x_grid, y_grid, num_grid, norm = norm, cmap = colormap)
    ax.set_zlabel(f"# particles m$^{{{-3}}}\\times 10^{{{max_order}}}$")
    ax.zaxis.labelpad = 5
    ax.set_box_aspect(aspect = None, zoom = 0.9)
    fig.tight_layout()
    ax.view_init(elev, azi)

    if exp_name != "":
        figure_name = f"Sensitivity_surface_{binName}_{time}s_{name_addition}.png" if name_addition != "" else f"sensitivity_surface_{binName}_{time}s.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)

    return fig, ax