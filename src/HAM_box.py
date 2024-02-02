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
import re
import pandas as pd
import cmocean as co
import netCDF4 as nc

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")


def read_input(filename):
    if filename in ["orig", "original"]:
        dataset = nc.Dataset(os.path.join(HAM_INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"))

    elif filename in ["kappa"]:
        dataset = nc.Dataset(os.path.join(HAM_BASE_FOLDER, "lut_kappa.nc"))
        
    # This bit is for safekeeping the original input files, just incase that ends up being necessary
    # I'm just going to put a copy of it in the /data/backup folder, idc if it's bad practice
    if not os.path.exists(os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc")):
        shutil.copy(os.path.join(HAM_INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"), os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc"))
    if not os.path.exists(os.path.join("../data/Backup", "lut_kappa.nc")):
        shutil.copy(os.path.join(HAM_BASE_FOLDER, "lut_kappa.nc"), os.path.join("../data/Backup", "lut_kappa.nc"))

    return dataset

def edit_input(original_dataset, new_filename):
    # This will probably also require us to delete the original input file, but I'm foregoing that for the moment.
    
    edited_dataset = nc.Dataset(os.path.join(HAM_INPUT_FOLDER, new_filename), "w", format = "NETCDF4")
    
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

def run_model(command = "./ham_box", recompile = True, verbose = True):
    # This is currently only setup to run commands in the HAM_box_OpenIFS directory.
    if recompile:
        command = f"cd ../.. && cd HAM_box_OpenIFS && make dest_dir=/src/*/ && {command}"
    else:
        command = f"cd ../.. && cd HAM_box_OpenIFS && {command}"
        
    try: 
        result = subprocess.run(["wsl", "bash", "-c", command], check = True,
                                stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        if verbose:
           print(f"{result.stderr.decode('ascii').strip()}")
    except subprocess.CalledProcessError as e:
        print(f"Error executing Linux command: {e}")

def copy_model_metadata(destination_folder):
    mo_ham = os.path.join(HAM_SRC_FOLDER, "mo_ham.f90")
    salsactl = os.path.join(HAM_SRC_FOLDER, "mo_ham_salsactl.f90")
    files = [mo_ham, salsactl]

    params = [{"name": "nham_subm", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "nseasalt", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "npist", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "naerorad", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "ndrydep", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "nwetdep", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "ndust", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "lscond", "type": "bool", "file": "mo_ham", "pattern": ""},
              {"name": "lscoag", "type": "bool", "file": "mo_ham", "pattern": ""},
              {"name": "nsoa", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "iaeroham", "type": "num", "file": "mo_ham", "pattern": ""},
              {"name": "nsnucl", "type": "num", "file": "mo_ham_salsactl", "pattern": ""},
              {"name": "locgas ", "type": "bool", "file": "mo_ham_salsactl", "pattern": ""},
              {"name": "lsol2b", "type": "bool", "file": "mo_ham_salsactl", "pattern": ""},
              {"name": "lhydr2b", "type": "bool", "file": "mo_ham_salsactl", "pattern": ""},
              ]
    
    # These are basic regex patterns, which we'll be formatting to catch those terms we're after
    base_num_pattern = r"INTEGER\s*(?:,\s*(PUBLIC|PARAMETER))?\s*::\s*{}\s*=\s*(\d+)"
    # base_bool_pattern = r"LOGICAL\s*(?:,\s*(PUBLIC|PARAMETER))?\s*::\s*{}\s*=\s*(\.TRUE\.|\.FALSE\.)"
    base_bool_pattern = r"(?:LOGICAL\s*)?(?:,\s*(PUBLIC|PARAMETER))?\s*(?:::)?\s*{}\s*=\s*(\.TRUE\.|\.FALSE\.)"
    patterndict = {"num": base_num_pattern,
                   "bool": base_bool_pattern}
    
    # Generate the specific regex patterns and write them into the dicts
    for paramdict in params:
        paramdict["pattern"] = patterndict[paramdict["type"]].format(paramdict["name"])

    # Start writing the metadata in a loop here - split these up between the various control files
    # Could maybe figure out a way to write that into a clever loop?
    mo_ham_params = [item for item in params if item["file"] == "mo_ham"]
    salsa_params = [item for item in params if item["file"] == "mo_ham_salsactl"]
    metadata = ""
                
    # Read mo_ham.f90 first
    with open(mo_ham) as file:
        lines = file.readlines()
    
    # Filter out commented lines because those are really messing with my regex
    lines = [line.strip() for line in lines if line.strip().startswith("!") == False]
    content = "\n".join(lines)
        
    for paramdict in mo_ham_params:
        paramstate = re.search(paramdict["pattern"], content).group(2)
        metadata += f"{paramdict['name']} = {paramstate}\n"
        
    # Read mo_ham_salsactl.f90 next
    with open(salsactl) as file:
        lines = file.readlines()
    
    # Filter out commented lines because those are really messing with my regex
    lines = [line.strip() for line in lines if line.strip().startswith("!") == False]
    content = "\n".join(lines)
        
    for paramdict in salsa_params:
        paramstate = re.search(paramdict["pattern"], content).group(2)
        metadata += f"{paramdict['name']} = {paramstate}\n"
        
    # Metadata is collected - bang it into a file
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)
    with open(os.path.join(full_destination_path, "metadata.txt"), "w") as metafile:
        metafile.write(metadata)
    

def copy_model_data(destination_folder):
    # Start by making sure the input is a string, to forego any funny business
    destination_folder = str(destination_folder)
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)
    
    # Check if the path already exists - if no, make it
    if not os.path.exists(full_destination_path):
        os.makedirs(full_destination_path)
        
    # Path exists now, so take the num.dat file from the data/ folder in HAM_box, and copy it to the new folder
    # We're assuming that if the path DOES exist, we're intending to overwrite it. This may come to bite us in the rear
    shutil.copy(os.path.join(HAM_DATA_FOLDER, "num.dat"), full_destination_path)
    copy_model_metadata(destination_folder)
        
def read_model_data(destination_folder):
    # Start by making sure the input is a string, to forego any funny business
    destination_folder = str(destination_folder)
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)
    
    # Check if the folder exist. If not, boot out immediately. Otherwise, read out the data.
    if not os.path.exists(full_destination_path):
        print("That folder does not exist! Check in the results directory, please.")
        return None
    
    else:
        num  = pd.read_csv(os.path.join(full_destination_path, "num.dat"),  sep=r"\s+")
        return num
        
    
def find_keyword(keyword, directory_path = HAM_SRC_FOLDER):
    # This function exists because I am sick and tired of looking through source files myself
    # Ensure the directory path is valid
    results = []
    if not os.path.exists(directory_path):
        print(f"Error: Directory '{directory_path}' not found.")
        return

    # Iterate through all files in the directory
    for filename in os.listdir(directory_path):
        if filename.endswith(".f90"):
            file_path = os.path.join(directory_path, filename)

            # Open and read the file
            with open(file_path, 'r') as file:
                # Iterate through each line in the file
                for line_number, line in enumerate(file, start=1):
                    # Check if the keyword is present in the line
                    if keyword in line:
                        result = {"filename": filename,
                                  "line_number": line_number,
                                  "line": line}
                        results.append(result)

    return results


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
    rdry, num, populations = ['a', 'b'],
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
    num  = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "num.dat"),  sep=r"\s+")
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry.dat"), sep=r"\s+")
    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works 
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100
    
    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    
    experiment_name = "NUCL0"
    plot_size_dist(rdry, num, rows=[1,200,1000, 2000, 4000, 7080], ymin=1, name_addition = experiment_name)
    plot_size_dist_evolution(rdry, num, vmin=1, name_addition = experiment_name)
