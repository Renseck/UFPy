# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:28:14 2024

The functions plot_size_dist, define_bin_boundaries, plot_size_dist_evolution were part of the original script, and
only slightly modified by me (Rens van Eck, BSc). The other functions were written by me.

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

from utils import order_of_magnitude, lognormal, parse_metadata

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")

title_fontsize = 15
label_fontsize = 14
tick_fontsize = 12 


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
        fig.text(0.47, -0.0, "Diameter (nm)", ha="center")

    if ncols > 1:
        label = kwargs.get("label", "placeholder")
        if label != "placeholder":
            kwargs.pop("label")
            
        for pop, ax in zip(populations, axes.ravel()):
            bins = [col for col in rdry.columns if pop in col]
            bins = [binName for binName in bins if binName in num.columns]
    
            for n in rows:
    
                r_row = rdry.iloc[n][bins]
                N_row = num.iloc[n][bins]

                ax.plot(r_row*1e9, N_row, label=f"{label + ' '}{n} s" if label == "" else f"{label}", **kwargs)
                
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
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
                
                if n == 1:
                    label = "Model init"
                    ax.plot(r_row*1e9, N_row, label=label, **kwargs)
                else:
                    ax.plot(r_row*1e9, N_row, label=f"{label + ' '}{n} s" if label == "Model" else f"{label}", **kwargs)
                
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
            # ax.set_xscale('log')
            # ax.set_yscale('log')
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)

            ax.legend(title="Data", bbox_to_anchor=(0.99, 0.99))
                
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
        fig=None, axes=None,
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

    if fig == None and axes == None:
        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(9*ncols, 6*nrows),
            sharex=True,
            sharey=True,
            dpi = 75
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
    
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
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
    
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
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
        colormap = None,
        highlights = None, highlight_colors = None):
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
    colormap : Colormap, optional
        Colormap of the data to be shown. The default is None.
    highlights : List, optional
        List of times to highlight by vertical lines. The default is None.
    highlight_colors : List, optional
        List of colors to highlight the times by. The default is None.

    Returns
    -------
    fig : TYPE
        Matplotlib figure object.
    axes : TYPE
        Matplotlib axes object.

    """
    num = num*1e-6
    
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
    
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], title = "Bins", bbox_to_anchor = (1.0, 1.02))
            
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
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
    
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles[::-1], labels[::-1], title = "Bins", bbox_to_anchor = (1.0, 1.02))
            
            if pop == "a":
                axtitle = "Population A: soluble particles"
            elif pop == "b":
                axtitle = "Population B: insoluble particles"
                
            ax.set_title(axtitle)
            ax.set_xlim(left=xmin, right=xmax)
            ax.set_ylim(bottom=ymin, top=ymax)
            
            if highlights != None:
                try:
                    iter(highlights)
                except Exception:
                    highlights = [highlights]
                
                for highlight, highlight_color in zip(highlights, highlight_colors):
                    ax.vlines(highlight, 0, num.max().sum(), linestyle = "dashed", color = highlight_color)

    if title != "":
        fig.suptitle(title, y = 1, x = 0.435)

    plt.tight_layout()
    # y-label
    fontsize = 13
    fig.text(-0.01, 0.5, "# particles cm$^{-3}$", va = "center", rotation = "vertical", fontsize = fontsize)
    ax.set_xlabel("Time (s)", fontsize = fontsize)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], title = "Bins", bbox_to_anchor = (1.0, 1.02))    

    if exp_name != "":
        figure_name = f"stacked_distribution_timeseries_{name_addition}.png" if name_addition != "" else "stacked_distribution_timeseries.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.savefig(full_savepath, bbox_inches = "tight", dpi = 150)
        
    return fig, axes

    

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
    environmental_keys = ["pt", "pqm1", "pap"]
    metadata_vals = []
    
    if type(binName) == str:
        for key, num in numdict.items():
            metadata = metadict[key]
            meta_parsed = parse_metadata(metadata)
            environmentals = [meta_parsed[key] for key in environmental_keys]
            flux = meta_parsed["particle_flux"]
            disp = meta_parsed["dispersion_rate"]
            
            environmentals.append(flux)
            environmentals.append(disp)
    
            numlist.append(num[binName].iloc[time])
            metadata_vals.append(environmentals)
      
    # If we input a list of bins, assume we want the total (sum) of it
    elif type(binName) == list:
        for key, num in numdict.items():
            metadata = metadict[key]
            meta_parsed = parse_metadata(metadata)
            environmentals = [meta_parsed[key] for key in environmental_keys]
            flux = meta_parsed["particle_flux"]
            disp = meta_parsed["dispersion_rate"]
    
            environmentals.append(flux)
            environmentals.append(disp)
    
            numlist.append(num[binName].sum(axis = 1).iloc[time])
            metadata_vals.append(environmentals)

    metadata_df = pd.DataFrame(data = metadata_vals)

    gridpoints = 1000
    temp_grid = np.linspace(metadata_df[0].min(), metadata_df[0].max(), num = gridpoints)
    humi_grid = np.linspace(metadata_df[1].min(), metadata_df[1].max(), num = gridpoints)
    pap_grid = np.linspace(metadata_df[2].min(), metadata_df[2].max(), num = gridpoints)
    
    
    # As a selection scheme, check if the environmentals are constant
    # If they're all constant, we must be varying the input flux/dispersion
    # If they're not all constant, we are varying environmentals.
    temp_constant = len(np.unique(temp_grid)) == 1
    humi_constant = len(np.unique(humi_grid)) == 1
    pap_constant = len(np.unique(pap_grid)) == 1
    all_env_constant = temp_constant and humi_constant and pap_constant
    
    # Not all environmentals are constant, which means we're varying them. Select axes as needed.
    if not all_env_constant:
        print("Detected environmental variation")
        constant_var = [temp_constant, humi_constant, pap_constant].index(True)
    
        # If temperature is constant
        if constant_var == 0:
            x_grid, y_grid = np.meshgrid(humi_grid, pap_grid)
            num_grid = griddata((metadata_df[1].values, metadata_df[2].values), numlist, (x_grid, y_grid), method = "linear")
            x_grid = x_grid*1e3
            y_grid = y_grid*1e-2
    
            ax.set_xlabel("q $(\\frac{kg}{kg})$ $\\times 10^{-3}$")
            ax.set_ylabel("Ambient pressure (hPa)")
            ax.yaxis.labelpad = 5
            constant_title = f"T = {metadata_df[0].iloc[0]:.2f} K"
    
        # If humidity is constant
        elif constant_var == 1:
            x_grid, y_grid = np.meshgrid(temp_grid, pap_grid)
            num_grid = griddata((metadata_df[0].values, metadata_df[2].values), numlist, (x_grid, y_grid), method = "linear")
            y_grid = y_grid*1e-2
            
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("Ambient pressure (hPa)")
            constant_title = f"q = {metadata_df[1].iloc[0]:.3e}" + " $(\\frac{kg}{kg})$"
    
        # If pressure is constant
        elif constant_var == 2:
            x_grid, y_grid = np.meshgrid(temp_grid, humi_grid)
            num_grid = griddata((metadata_df[0].values, metadata_df[1].values), numlist, (x_grid, y_grid), method = "linear")
            y_grid = y_grid*1e3
    
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("q $(\\frac{kg}{kg})$ $\\times 10^{-3}$")
            constant_title = f"p = {metadata_df[2].iloc[0]:.3e} Pa"
            
    # All environmentals _are_ constant, so we're varying the influx/dispersion.
    # Note there are about 10239210392130 ways of varying the input;
    # - are we varying dispersion?
    # - are we keeping the distribution shape equal but only increasing number?
    # - are we shifting the distribution peak around?
    elif all_env_constant:
        print("Detected flux/dispersion variation")
        dispersion_grid = np.linspace(metadata_df[4].min(), metadata_df[4].max(), num = gridpoints)
        
        # Don't ask me to explain how this works, it may be implicitly superfluous anyway
        flux_constant = len(set(frozenset(item) for item in metadata_df[3].values)) == 1
        disp_constant = len(np.unique(dispersion_grid)) == 1
        
        # I'm maybe going to only implement differentiation between two kinds of non-constant flux:
        # 1. The shape of the distribution is invariant, but the position (peak) shifts
        #    around (e.g. from 20nm to 100nm)
        # 2. The position/shape is invariant, but the number of of particles increases by some factor A
        #    in the way of A * lognormal(Dp)
        # Dispersion is assumed to be varied always too, or we can't have 2D surfaces. 
        
    
    # num_grid = num_grid*1e-6  # Convert from #/m^-3 to #/cm^-3

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
    
    ax.set_title(title, y = 0.95)
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
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)

    return fig, ax

def plot_dispersion_variation(numdict, metadict, rdry, exp_name = "", binName = "2a1", time = 1500, name_addition = ""):
    """
    Use this function to plot variation series where ONLY dispersion was changed

    Parameters
    ----------
    numdict : Dictionary
        Dictionary of Pandas Dataframes (like num in other functions).
    metadict : Dictionary
        Dictinary of metadata strings.
    rdry: DataFrame
        DataFrame containing bin boundaries per time step.
    exp_name : String, optional
        Name of the experiment. Leave empty to forego saving the image. The default is "".
    binName : String, optional.
        Name of bin to plot in dispersion relation thing. Default is "2a1".
    time : Integer, optional
        Time to be plotted. The default is 1500.
    name_addition : String, optional
        File name addition, for easy lookup. The default is "".

    Returns
    -------
    None.

    """
    salsa_bins = ["1a1", "1a2", "1a3", "2a1", "2a2", "2a3", "2a4", "2a5"]
    
    environmental_keys = ["pt", "pqm1", "pap"]
    metadata_vals = []
    numlist = []
    
    for key, num in numdict.items():
        metadata = parse_metadata(metadict[key])
        env = [metadata[key] for key in environmental_keys]
        flux = metadata["particle_flux"]
        disp = metadata["dispersion_rate"]
        
        env.append(flux)
        env.append(disp)
        
        numlist.append(num.iloc[time])
        metadata_vals.append(env)
        
    metadata_df = pd.DataFrame(data = metadata_vals)
    numdata_df = pd.DataFrame(data = numlist)
    numdata_df.index = range(len(numlist))

    total_df = numdata_df.join(metadata_df).rename(columns = {0: "pt", 1: "pqm1", 2: "pap", 3: "particle_flux", 4: "dispersion"})
    total_df = total_df.sort_values(by = "dispersion")

    colormap = co.cm.dense
    colors = colormap(np.linspace(0, 1, len(numdict)))

    # Plot full distributions for various dispersions
    plt.figure(figsize = (10,6))

    for ind, dispersion in enumerate(total_df["dispersion"].values):
        plt.plot(rdry[salsa_bins].iloc[time]*1e9, total_df[salsa_bins].iloc[ind]*1e-6, label = f"d = {dispersion*100:.2f}%", color = colors[ind],  alpha = 1)
        
    plt.legend(title = "Dispersion")
    plt.xlim((0, 400))
    plt.xlabel("Diameter (nm)", fontsize = label_fontsize)
    plt.ylabel("# particles cm$^{-3}$", fontsize = label_fontsize)
    
    T = total_df["pt"].iloc[0]
    q = total_df["pqm1"].iloc[0]
    p = total_df["pap"].iloc[0]
    plt.title(f"Number-size distribution after {time} s \n [T = {T:.1f} K, q = {q:.2e} $\\frac{{kg}}{{kg}}$, p = {p/100:.0f} hPa]",
              fontsize = title_fontsize)
    plt.tick_params(axis = "both", labelsize = tick_fontsize) 
    
    if exp_name != "":
        figure_name = f"Dispersion_sensitivity_{time}s_{name_addition}.png" if name_addition != "" else f"Dispersion_sensitivity_{time}s.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)
        
    # Plot the relation between dispersion (x-axis) and population of a certain bin (y-axis)
    bin_boundaries = define_bin_boundaries()
    binName_split = re.split("(?<=[a-zA-Z])(?=\d)", binName)
    lower_boundary = bin_boundaries[binName_split[0]][int(binName_split[1])-1]
    upper_boundary = bin_boundaries[binName_split[0]][int(binName_split[1])]
    
    plt.figure(figsize = (10,6))
    # The line
    plt.plot(total_df["dispersion"]*100, total_df[binName]*1e-6)
    # See if there's a maximum in the graph
    differences = total_df[binName].diff()
    max_index = total_df.index[differences.eq(0)]
    
    # There's not a literal 0 in the gradient, so try to find one by hand - not perfect.
    if max_index.empty:
        max_index = differences[(differences.shift(1) > 0) & (differences <= 0)].index
        
    max_dispersion = total_df["dispersion"].iloc[max_index]
    max_bincount = total_df[binName].iloc[max_index]

    if max_index.empty == False:
        plt.scatter(max_dispersion*100, max_bincount*1e-6, color = "red")
        plt.vlines(max_dispersion*100, ymin = 0.8*total_df[binName].min()*1e-6, ymax = total_df[binName].max()*1e-6, color = "red", linestyle  = "dashed")
        plt.hlines(max_bincount*1e-6, xmin = 0.8*total_df["dispersion"].min()*100, xmax = max_dispersion*100, color ="red", linestyle = "dashed")
        plt.text(1.01*max_dispersion*100, 0.9*max_bincount*1e-6,
                 s = f"{max_dispersion.values[0]*100:.3f}%", color = 'red', fontsize = label_fontsize,
                 bbox = dict(facecolor = "None", edgecolor = "black", boxstyle = 'round,pad=0.3'))
        
    
    plt.title(f"Dispersion vs particle count in bin {lower_boundary*1e9:.2f} - {upper_boundary*1e9:.2f} nm after {time} s\n [T = {T:.1f} K, q = {q:.2e} $\\frac{{kg}}{{kg}}$, p = {p/100:.0f} hPa]",
              fontsize = title_fontsize)
    plt.xlabel("Dispersion rate $d$ (%)", fontsize = label_fontsize)
    plt.ylabel("# particles cm$^{-3}$", fontsize  = label_fontsize)
    plt.tick_params(axis = "both", labelsize = tick_fontsize)
    plt.ylim((0.99*total_df[binName].min()*1e-6, 1.01*total_df[binName].max()*1e-6))
    plt.xlim((0.99*total_df["dispersion"].min()*100, 1.01*total_df["dispersion"].max()*100))
    
    if exp_name != "":
        figure_name = f"Dispersion_vs_{binName}_{time}s_{name_addition}.png" if name_addition != "" else f"Dispersion_vs_{binName}_{time}s.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)

def plot_flux_peak_variation(numdict, metadict, exp_name = "", binName = "", time = 1500, name_addition = ""):
    """
    Use this function to plot variation series where ONLY the peak of the particle influx was shifted

    Parameters
    ----------
    numdict : Dictionary
        Dictionary of Pandas Dataframes (like num in other functions).
    metadict : Dictionary
        Dictinary of metadata strings.
    exp_name : String, optional
        Name of the experiment. Leave empty to forego saving the image. The default is "".
    time : Integer, optional
        Time to be plotted. The default is 1500.
    binName : String, optional
        Name of the bin to be plotted. Entering multitple sums the bins together. The default is "1a1".
    name_addition : String, optional
        File name addition, for easy lookup. The default is "".

    Returns
    -------
    None.

    """
    salsa_bins = ["1a1", "1a2", "1a3", "2a1", "2a2", "2a3", "2a4"]

    environmental_keys = ["pt", "pqm1", "pap"]
    metadata_vals = []
    numlist = []

    for key, num in numdict.items():
        metadata = parse_metadata(metadict[key])
        env = [metadata[key] for key in environmental_keys]
        flux = metadata["particle_flux"]
        disp = metadata["dispersion_rate"]

        env.append(flux)
        env.append(disp)

        numlist.append(num.iloc[time])
        metadata_vals.append(env)

    metadata_df = pd.DataFrame(data = metadata_vals)
    numdata_df = pd.DataFrame(data = numlist)
    numdata_df.index = range(len(numlist))

    total_df = numdata_df.join(metadata_df).rename(columns = {0: "pt", 1: "pqm1", 2: "pap", 3: "particle_flux", 4: "dispersion"})

    # It has come to my attention that the loop above jumbles the order of the runs, which is a serious problem
    # when you're relying on the order being conserved. This little snippet of code attempts to resolve the order
    # but I cannot guarantee that it'll work for every case, nor have I tested it because I didn't care to :)
    input_list = [np.argmax(flux) for flux in metadata_df[3]]
    output_list = np.zeros_like(input_list)

    pos = 0
    last_number = 0

    for k in range(2):
        for i, number in enumerate(input_list):
            if number == last_number + 1:
                for j in range(i, len(input_list)):
                    if input_list[j] == number and output_list[j] == 0:
                        output_list[j] = pos
                        pos += 1
                last_number += 1
        last_number -= 1

    total_df["center_x"] = output_list
    total_df = total_df.sort_values(by = "center_x")
    
    plt.figure(figsize = (10,6))

    if binName == "":
        for binName in salsa_bins:
            plt.plot(total_df["center_x"], total_df[binName]*1e-6, label = f"{binName}", alpha = 0.75)
        title = f"Sensitivity to shifting input flux peak after {time} s"
        plt.legend(title = "Bin")

    else:
        plt.plot(total_df["center_x"], total_df[binName]*1e-6)
        title = f"Sensitivity to shifting input flux peak after {time} s in bin {binName}"

    plt.xlabel("Peak of lognormal distribution (nm)")
    plt.ylabel("# particles cm$^{-3}$")
    plt.title(title)
    
    if exp_name != "":
        figure_name = f"Distribution_peak_sens_{time}s_{name_addition}.png" if name_addition != "" else f"Distribution_peak_sens_{time}s.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)

def plot_flux_scale_variation(numdict, metadict, exp_name = "", binName = "", time = 1500, name_addition = ""):
    """
    Use this function to plot variation series where ONLY the scale of the particle influx was shifted

    Parameters
    ----------
    numdict : Dictionary
        Dictionary of Pandas Dataframes (like num in other functions).
    metadict : Dictionary
        Dictinary of metadata strings.
    exp_name : String, optional
        Name of the experiment. Leave empty to forego saving the image. The default is "".
    time : Integer, optional
        Time to be plotted. The default is 1500.
    binName : String, optional
        Name of the bin to be plotted. Entering multitple sums the bins together. The default is "1a1".
    name_addition : String, optional
        File name addition, for easy lookup. The default is "".

    Returns
    -------
    None.

    """
    salsa_bins = ["1a1", "1a2", "1a3", "2a1", "2a2", "2a3", "2a4"]

    environmental_keys = ["pt", "pqm1", "pap"]
    metadata_vals = []
    numlist = []

    for key, num in numdict.items():
        metadata = parse_metadata(metadict[key])
        env = [metadata[key] for key in environmental_keys]
        flux = metadata["particle_flux"]
        disp = metadata["dispersion_rate"]

        env.append(flux)
        env.append(disp)

        numlist.append(num.iloc[time])
        metadata_vals.append(env)

    metadata_df = pd.DataFrame(data = metadata_vals)
    numdata_df = pd.DataFrame(data = numlist)
    numdata_df.index = range(len(numlist))

    total_df = numdata_df.join(metadata_df).rename(columns = {0: "pt", 1: "pqm1", 2: "pap", 3: "particle_flux", 4: "dispersion"})
    total_df["scale"] = [int(total_df["particle_flux"].iloc[ind][1] / total_df["particle_flux"].iloc[0][1]) for ind in range(0, len(total_df))]
    total_df = total_df.sort_values(by = "scale")

    plt.figure(figsize = (10,6))

    if binName == "":
        for binName in salsa_bins:
            plt.plot(total_df["scale"], total_df[binName]*1e-6, label = f"{binName}", alpha = 0.75)
        title = f"Sensitivity to scale input flux magnitude after {time} s"
        plt.legend(title = "Bin")

    else:
        plt.plot(total_df["scale"], total_df[binName]*1e-6)
        title = f"Sensitivity to scaling input flux magnitude after {time} s in bin {binName}"

    plt.xlabel("Scaling factor A")
    plt.ylabel("# particles cm$^{-3}$")
    plt.title(title)
    
    if exp_name != "":
        figure_name = f"Distribution_peak_sens_{time}s_{name_addition}.png" if name_addition != "" else f"Distribution_peak_sens_{time}s.png"
        savepath = os.path.join(MODEL_PLOT_FOLDER, exp_name)

        if not os.path.exists(savepath):
            os.makedirs(savepath)

        full_savepath = os.path.join(savepath, figure_name)
        plt.tight_layout()
        plt.savefig(full_savepath, bbox_inches='tight', pad_inches=0, dpi = 150)

def show_model_lognormal_flux(sigma = 1.68, center_x = 90, scale = 1, label = "", file_addition = ""):
    """
    Generates a plot of the lognormal distribution, with correpsonding bars for what it looks like translated
    into SALSA2.0 bins. 

    Parameters
    ----------
    sigma : Float, optional
        Geometric standard deviation of the lognormal distribution. The default is 1.68.
    center_x : Float, optional
        X-coord of the peak of the distribution. The default is 90.
    scale : Float, optional
        Factor to multiply the whole distribution with. The default is 1.
    label : String, optional
        Label for the plot, shown in legend. The default is "".
    file_addition : String, optional
        Label for the file, for better tracking in file system. The default is "".

    Returns
    -------
    None.

    """
    if scale == 1:
        air_density = 1
    else:
        air_density = 1.2
        
    x = np.linspace(1, 1000, 10000)

    bin_boundaries = define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)

    y = lognormal(x, sigma, center_x = center_x, scale = scale)
    particle_flux = np.interp(salsa_bins, x, y)

    bar_widths = np.append(np.diff(salsa_bins), 0)
    bar_positions = salsa_bins + bar_widths / 2

    plt.plot(x, y/air_density, label = f"Lognormal{' ' + label}")
    plt.bar(bar_positions, particle_flux/air_density, width = bar_widths, label = "SALSA2.0 flux", color = "orange",
            alpha = 0.6, align = 'center', edgecolor = "black"  )
    #plt.plot(salsa_bins, diesel_flux)
    #plt.scatter(salsa_bins, diesel_flux, color = "orange" )
    plt.xlim((3, 1000))
    plt.xscale('log')
    
    plt.title("Model emissions input")
    plt.xlabel("Diameter (nm)")
    # plt.ylabel("n")
    plt.ylabel("# particles cm$^{-3}$ s$^{-1}$")

    plt.legend()
    
    file_name = f"Model_lognormal_flux_bars_{file_addition}.jpg" if file_addition != "" else "Model_lognormal_flux_bars.jpg"
    plt.savefig(os.path.join(RESULTS_FOLDER, f"Measurement Figures/{file_name}"), dpi = 150)
    plt.show()
    
def show_normalised_lognormal_flux(sigma = 1.68, center_x = 90, scale = 1, title = "", label = "", file_addition = ""):        
    """
    Generates the normalized version of the lognormal distribution as reported by Harris and Maricq, 2001.
    With correpsonding bars for what it looks like translated into SALSA2.0 bins. 

    Parameters
    ----------
    sigma : Float, optional
        Geometric standard deviation of the lognormal distribution. The default is 1.68.
    center_x : Float, optional
        X-coord of the peak of the distribution. The default is 90.
    scale : Float, optional
        Factor to multiply the whole distribution with. The default is 1.
    title : String, optional
        Title of the figure. The default is "". 
    label : String, optional
        Label for the plot, shown in legend. The default is "".
    file_addition : String, optional
        Label for the file, for better tracking in file system. The default is "".

    Returns
    -------
    None.

    """
    x = np.linspace(1, 1000, 10000)

    bin_boundaries = define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)

    y = lognormal(x, sigma, center_x = center_x, scale = scale)
    particle_flux = np.interp(salsa_bins, x, y)

    bar_widths = np.append(np.diff(salsa_bins), 0)
    bar_positions = salsa_bins + bar_widths / 2

    plt.plot(x, y, label = f"Lognormal{' ' + label}")
    plt.bar(bar_positions, particle_flux, width = bar_widths, label = "SALSA2.0 flux", color = "orange",
            alpha = 0.6, align = 'center', edgecolor = "black"  )
    #plt.plot(salsa_bins, diesel_flux)
    #plt.scatter(salsa_bins, diesel_flux, color = "orange" )
    plt.xlim((3, 1000))
    plt.xscale('log')
    
    if title == "":
        plt.title("Model emissions input")
    else:
        plt.title(title)
        
    plt.xlabel("Diameter (nm)")
    plt.ylabel("n")
    # plt.ylabel("# particles cm$^{-3}$ s$^{-1}$")

    plt.legend()
    
    file_name = f"Normalised_lognormal_flux_bars_{file_addition}.jpg" if file_addition != "" else "Normalised_lognormal_flux_bars.jpg"
    plt.savefig(os.path.join(RESULTS_FOLDER, f"Measurement Figures/{file_name}"), dpi = 150)
    plt.show()
    
def show_model_flux(particle_flux, file_addition = ""):
    """
    Shows any model flux given as argument, shown as bars. 

    Parameters
    ----------
    particle_flux : List
        List of particle fluxes.
    file_addition : String, optional
        Label for the file, for better tracking in file system. The default is "".

    Returns
    -------
    None.

    """
    air_density = 1.2
    
    bin_boundaries = define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    bar_widths = np.append(np.diff(salsa_bins), 0)
    bar_positions = salsa_bins + bar_widths / 2

    plt.bar(bar_positions, particle_flux/air_density, width = bar_widths, label = "SALSA2.0 flux", color = "orange",
            alpha = 0.6, align = 'center', edgecolor = "black"  )
    
    plt.xlim((3, 1000))
    plt.xscale('log')
    
    plt.title("Model emissions input")
    plt.xlabel("Diameter (nm)")
    # plt.ylabel("n")
    plt.ylabel("# particles cm$^{-3}$ s$^{-1}$")

    plt.legend()
    
    file_name = f"Model_flux_bars_{file_addition}.jpg" if file_addition != "" else "Model_flux_bars.jpg"
    plt.savefig(os.path.join(RESULTS_FOLDER, f"Measurement Figures/{file_name}"), dpi = 150)
    plt.show()
