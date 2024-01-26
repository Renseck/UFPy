# =============================================================================
# This file functions to analyse the output of the *SMALL* SALSA2.0 Standalone model, which DOES NOT include
# atmospheric chemistry and cloud formation. There is no GitHub repository or documentation for this distribution.
# Investigation of the source code (driver.f90) seems to indicate that number concentration are
# currently written to output in num.dat.
# 
# The functions plot_size_dist, define_bin_boundaries, plot_size_dist_evolution were part of the original script, and
# only slightly modified by me (Rens van Eck, BSc). The other functions were written by me
# =============================================================================

import os
import subprocess
import shutil
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors as mc
import pandas as pd
import cmocean as co
import netCDF4 as nc
import xarray as xr

BASE_FOLDER = "../../HAM_box_OpenIFS"
INPUT_FOLDER = os.path.join(BASE_FOLDER, "input")
DATA_FOLDER = os.path.join(BASE_FOLDER, "data")

def read_input(filename):
    if filename in ["orig", "original"]:
        dataset = nc.Dataset(os.path.join(INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"))

    elif filename in ["kappa"]:
        dataset = nc.Dataset(os.path.join(BASE_FOLDER, "lut_kappa.nc"))
        
    # This bit is for safekeeping the original input files, just incase that ends up being necessary
    # I'm just going to put a copy of it in the /data/backup folder, idc if it's bad practice
    if not os.path.exists(os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc")):
        shutil.copy(os.path.join(INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"), os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc"))
    if not os.path.exists(os.path.join("../data/Backup", "lut_kappa.nc")):
        shutil.copy(os.path.join(BASE_FOLDER, "lut_kappa.nc"), os.path.join("../data/Backup", "lut_kappa.nc"))

    return dataset

def edit_input(original_dataset, new_filename):
    # This will probably also require us to delete the original input file, but I'm foregoing that for the moment.
    
    edited_dataset = nc.Dataset(os.path.join(INPUT_FOLDER, new_filename), "w", format = "NETCDF4")
    
    # Copy dimensions from the original to the edited dataset
    for dim_name, dim_obj in original_dataset.dimensions.items():
        edited_dataset.createDimension(dim_name, len(dim_obj))

    # Copy variables and their attributes from the original to the edited dataset
    for var_name, var_obj in original_dataset.variables.items():
        edited_var = edited_dataset.createVariable(var_name, var_obj.dtype, var_obj.dimensions)
        edited_var[:] = var_obj[:]  # Copy variable values
        for attr_name in var_obj.ncattrs():
            edited_var.setncattr(attr_name, var_obj.getncattr(attr_name))  # Copy variable attributes

    # Modify the values of a variable in the edited dataset
    # This bit is going to include some Latin Hypercube Sampling (LHS) method to explore parameter space,
    # Which will involve rewriting every variable one-by-one. 
    edited_dataset.variables['your_variable'][:] = np.ones_like(edited_dataset.variables['your_variable'][:]) * 42

    return edited_dataset

def run_linux_command(command, verbose = True):
    # This is currently only setup to work with running either the "salsa_box" or the "ham_box" scripts.
    # I also have no clue what the difference between these two is.
    try: 
        command = f"cd ../.. && cd HAM_box_OpenIFS && ./{command}"
        result = subprocess.run(["wsl", "bash", "-c", command], check = True,
                                stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if verbose:
           print(f"{result.stderr.decode('ascii').strip()}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing Linux command: {e}")
        

def plot_size_dist(
    rdry, num, rows = [0], populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    name_addition = ""
    ):

    ## make sure that row_nr is a list-like object
    try:
        iter(rows)
    except Exception:
        rows = [rows]

    nrows = 1
    ncols = len(populations)
    
    fig, axes = plt.subplots(
        nrows = nrows,
        ncols = ncols,
        figsize = (6*ncols, 4*nrows),
        sharex = True,
        sharey = True,
    )
    plt.tight_layout()
    fig.text(-0.01, 0.5, "# particles cm$^{-3}$", va = "center", rotation = "vertical")
    fig.text(0.5, 0.04, "Diameter (m)", ha = "center")

    for pop,ax in zip(populations, axes.ravel()):
        bins = [col for col in rdry.columns if pop in col]

        for n in rows:

            r_row = rdry.iloc[n][bins]
            N_row = num.iloc[n][bins]

            ax.plot(r_row, N_row, label=n)
        ax.set_title(f"Population {pop}")
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(left=xmin, right=xmax)
        ax.set_ylim(bottom=ymin, top=ymax)
        if pop == "b":
            ax.legend(title = "Time", bbox_to_anchor = (1.2, 1.02))
            
    figure_name = "size_distribution{addition}.png".format(addition = f"_{name_addition}" if name_addition != "" else "")
    plt.savefig(os.path.join("../results", figure_name), bbox_inches = 'tight', pad_inches = 0)
    plt.show()


def define_bin_boundaries(populations = ['1a', '2a', '2b']):

    return {
        pop : (
            np.logspace(np.log10(3e-9), np.log10(50e-9), 4) if pop[0]=='1' else
            np.concatenate([
                np.logspace(np.log10(50e-9), np.log10(700e-9), 5, endpoint=True)[:-1],
                np.logspace(np.log10(700e-9), np.log10(10000e-9),4, endpoint=True),
            ])
        )
        for pop in populations
    }


def plot_size_dist_evolution(
    num, populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    vmin = None, vmax = None,
    name_addition = ""
    ):

    bin_boundaries = define_bin_boundaries()

    nrows = len(populations)
    ncols = 1
    
    fig, axes = plt.subplots(
        nrows = nrows,
        ncols = ncols,
        figsize = (12*ncols, 4*nrows),
        sharex = True,
        sharey = True,
    )

    ## computing common color bar bounds:
    vmin = num.values.min() if vmin is None else vmin
    vmax = num.values.max() if vmax is None else vmax

    for pop,ax in zip(populations, axes.ravel()):
        bins = [col for col in rdry.columns if pop in col]

        ## combining bin boundaries (omitting double values)
        rbounds = [bounds for k,bounds in sorted(bin_boundaries.items()) if pop in k]
        rbounds = np.concatenate([rb[:-1] for rb in rbounds[:-1]]+[rbounds[-1]])

        ## time bounds
        tbounds = np.arange(num.shape[0]+1)

        ## generating meshgrid
        t,r = np.meshgrid(tbounds,rbounds)
      
        norm = mc.LogNorm(vmin = vmin, vmax = vmax)
        cls = ax.pcolormesh(t, r, num[bins].T, norm = norm, cmap = co.cm.dense)

        ax.set_title(f"Population {pop}")
        ax.set_yscale('log')
        ax.set_ylabel("Diameter (m)")
        ax.set_xlim(left=xmin, right=xmax)
        ax.set_ylim(bottom=ymin, top=ymax)
        
    fig.colorbar(cls, ax=axes.ravel().tolist())
    figure_name = "size_distribution_LES_box{addition}.png".format(addition = f"_{name_addition}" if name_addition != "" else "")
    plt.savefig(os.path.join("../results", figure_name), bbox_inches = 'tight', pad_inches = 0)
    ax.set_xlabel("Time")
    # plt.close()
    plt.show()

        
if __name__ == '__main__':
    num  = pd.read_csv(os.path.join(DATA_FOLDER, "num.dat"),  sep=r"\s+")
    rdry = pd.read_csv(os.path.join(DATA_FOLDER, "rdry.dat"), sep=r"\s+")
    plot_size_dist(rdry, num, rows=[1,200,1000, 2000, 4000, 7080], ymin=1)
    plot_size_dist_evolution(num, vmin=1)
