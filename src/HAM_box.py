# =============================================================================
# This file functions to analyse the output of the *SMALL* HAM box model, mostly using the salsa2.0 package, 
# which DOES NOT include atmospheric chemistry and cloud formation.
# There is no GitHub repository or documentation for this distribution.
# Investigation of the source code (driver.f90) seems to indicate that number concentration are
# currently written to output in num.dat.
#
# The functions plot_size_dist, define_bin_boundaries, plot_size_dist_evolution were part of the original script, and
# only slightly modified by me (Rens van Eck, BSc). The other functions were written by me.
# =============================================================================

import json
import os
# from matplotlib import pyplot as plt
# from matplotlib import colors as mc
# from scipy.interpolate import griddata
import re
# from itertools import product
import shutil
import subprocess

import netCDF4 as nc
import numpy as np
import pandas as pd
from tqdm import tqdm

import HAM_plot as hp

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")


def read_input(filename):
    # maybe pointless/deprecated function
    if filename in ["orig", "original"]:
        dataset = nc.Dataset(os.path.join(HAM_INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"))

    elif filename in ["kappa"]:
        dataset = nc.Dataset(os.path.join(HAM_BASE_FOLDER, "lut_kappa.nc"))

    # This bit is for safekeeping the original input files, just incase that ends up being necessary
    # I'm just going to put a copy of it in the /data/backup folder, idc if it's bad practice
    if not os.path.exists(os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc")):
        shutil.copy(os.path.join(HAM_INPUT_FOLDER, "HAM_box_inp_200007.01_activ.nc"),
                    os.path.join("../data/Backup", "HAM_box_inp_200007.01_activ.nc"))
    if not os.path.exists(os.path.join("../data/Backup", "lut_kappa.nc")):
        shutil.copy(os.path.join(HAM_BASE_FOLDER, "lut_kappa.nc"), os.path.join("../data/Backup", "lut_kappa.nc"))

    return dataset


def edit_input(original_dataset, new_filename):
    # maybe pointless function
    # This will probably also require us to delete the original input file, but I'm foregoing that for the moment.

    edited_dataset = nc.Dataset(os.path.join(HAM_INPUT_FOLDER, new_filename), "w", format="NETCDF4")

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


def run_model(experiment_name, recompile=True, verbose=True):
    """
    Runs the HAM_box_OpenIFS model. Make sure it is contained within the folder containing the UFPy folder.

    Parameters
    ----------
    experiment_name : STR
        Name of the folder to which the data should be written (in results/).
    recompile : BOOL, optional
        Recompile the model before running. The default is True.
    verbose : BOOL, optional
        Output error messages encountered in terminal. The default is True.

    Returns
    -------
    None.

    """
    # Clear out any old num.dat files first, so that there's no annoying stragglers
    num_files = [file for file in os.listdir(HAM_DATA_FOLDER) if file.startswith("num")]
    num_files.remove("num_m7.dat")
    for file in num_files:
        os.remove(os.path.join(HAM_DATA_FOLDER, file))
    
    # This is currently only setup to run commands in the HAM_box_OpenIFS directory.
    run_command = "./ham_box"
    if recompile:
        command = f"cd ../.. && cd HAM_box_OpenIFS && make dest_dir=/src/*/ && {run_command}"
    else:
        command = f"cd ../.. && cd HAM_box_OpenIFS && {run_command}"

    try:
        result = subprocess.run(["wsl", "bash", "-c", command], check=True,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if verbose:
            print(f"OUTPUT: \n{result.stdout.decode('ascii').strip()}")
            print(f"ERRORS: \n{result.stderr.decode('ascii').strip()}")

        copy_model_data(experiment_name)

    except subprocess.CalledProcessError as e:
        print(f"Error executing Linux command: {e}")


def run_variation(variation_name, nested_environmental_values):
    """
    Execute multiple model runs, with varied environental values

    Parameters
    ----------
    variation_name : String
        Name of the "overarching" experiment.
    nested_environmental_values : List of list
        Nested list / matrix of environmental values.

    Returns
    -------
    None.

    """
    defaults = [298, 0.0058535, 101325]

    # Make the "parent" folder for the data to be written into
    if not os.path.exists(os.path.join(MODEL_L0_FOLDER, variation_name)):
        os.mkdir(os.path.join(MODEL_L0_FOLDER, variation_name))

    with tqdm(total=len(nested_environmental_values), desc="Processing", leave=True, position=0) as pbar:
        for run, environmental_values in enumerate(nested_environmental_values):
            write_environmental_data(environmental_values)
            pt, pqm1, pap = environmental_values

            description =  f"Processing [{pt:.2f}, {pqm1:.3e}, {pap:.3e}]"
            description = description.ljust(45, " ")
            pbar.set_description(desc = description)
            experiment_name = os.path.join(variation_name, str(run))
            # Recompile on the first (zeroth) run just in case
            if run == 0:
                run_model(experiment_name, recompile=True, verbose=False)
            else:
                run_model(experiment_name, recompile=False, verbose=False)

            pbar.update(1)
    # After the loop, reset to defaults
    write_environmental_data(defaults)


def gen_params():
    """
    This function reads the params.json file, to see which parameters need to be kept track of

    Returns
    -------
    paramlist : LIST
        List of dictionaries containing regex patterns for each parameter.
    files : LIST
        List of files containing the parameter info.

    """
    with open("params.json", "r") as paramfile:
        paramjson = json.load(paramfile)

    paramlist = paramjson["params"]
    filenames = list(set(dict["file"] for dict in paramlist))
    files = [os.path.join(HAM_SRC_FOLDER, filename + ".f90") for filename in filenames]

    # These are basic regex patterns, which we'll be formatting to catch those terms we're after
    base_num_pattern = r"(?:INTEGER\s*)?(?:,\s*(PUBLIC|PARAMETER))?\s*(?:::)?\s*{}\s*=\s*(\d+)"
    base_bool_pattern = r"(?:LOGICAL\s*)?(?:,\s*(PUBLIC|PARAMETER))?\s*(?:::)?\s*{}\s*=\s*(\.TRUE\.|\.FALSE\.)"
    patterndict = {"num": base_num_pattern,
                   "bool": base_bool_pattern
                   }

    # Generate the specific regex patterns and write them into the dicts
    for item in paramlist:
        item["pattern"] = patterndict[item["type"]].format(item["name"])

    return paramlist, files


def write_environmental_data(environmental_vals):
    """
    Writes values to the environmental.dat file **IN ORDER**
    [pt, pqm1, pap]

    Parameters
    ----------
    environmental_vals : LIST
        Contains the values of the environmental values **IN ORDER**.

    Returns
    -------
    None.

    """
    # The order of these variables is absolutely crucial - take care
    # Ambient temperature, specific humidity, ambient pressure
    
    # Check if the list is multi-dimensional, to allow for non-uniform environmentals.
    # This checks positively - if the list contains lists...
    content = ""
    if any(isinstance(el, list) for el in environmental_vals):
        for line in environmental_vals:
            content += " ".join(str(val) for val in line) # Same functions, just add a space at the end
            content += "\n"
    
    else:
         content = " ".join(str(val) for val in environmental_vals)

    with open(os.path.join(HAM_INPUT_FOLDER, "environmental.dat"), "w") as outfile:
        outfile.write(content)


def read_environmental_data():
    """
    Reads the environmental.dat file, to see which values have been used to run the model last.

    Returns
    -------
    metadata : STRING
        Formatted string of model environmental data.

    """
    # The order of the variables is absolutely crucial - take care
    # Ambient temperature, specific humidity, ambient pressure
    environmental_vars = ["pt", "pqm1", "pap"]
    metadata = ""

    with open(os.path.join(HAM_INPUT_FOLDER, "environmental.dat"), "r") as infile:
        content = infile.read()

    environmental_vals = content.split(" ")

    for var, val in zip(environmental_vars, environmental_vals):
        metadata += f"{var} = {val}\n"

    return metadata


def read_model_metadata():
    """
    Reads CURRENT metadata from all files found in the params.json file

    Parameters
    ----------
    file_list : LIST
        List of names of files (without .f90) to be read from.
    paramdict_list : LIST
        List of dictionaries containing parameter info to be read (see copy_model_metadata() for more info).

    Returns
    -------
    metadata : STRING
        Formatted string of model metadata.

    """

    metadata = ""

    paramdict_list, file_list = gen_params()

    for file in file_list:

        with open(file) as f:
            lines = f.readlines()

        lines = [line.strip() for line in lines if line.strip().startswith("!") == False]
        content = "\n".join(lines)

        for paramdict in paramdict_list:
            matches = re.search(paramdict["pattern"], content)
            if matches == None:
                continue
            else:
                paramstate = matches.group(2)
                metadata += f"{paramdict['name']} = {paramstate}\n"

    # Read environmental variables from .dat file, because we can't retrieve those with regex
    metadata += read_environmental_data()

    return metadata


def copy_model_metadata(destination_folder):
    """
    Copies the model metadata to metadata.txt.

    Parameters
    ----------
    destination_folder : STR
        Name of the folder to which metadata.txt should be written.

    Returns
    -------
    None.

    """
    # Read out the metadata from the various files
    metadata = read_model_metadata()

    # Metadata is collected - bang it into a file
    # Start by making sure the input is a string, to forego any funny business
    destination_folder = str(destination_folder)
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)
    with open(os.path.join(full_destination_path, "metadata.txt"), "w") as metafile:
        metafile.write(metadata)


def check_metadata():
    """
    Loops through all model result folders and reads their metadata.

    Returns
    -------
    metadict : DICT
        Dictionary, with keys the name of every result folder, and values being the parameters of each model run folder.

    """
    metadict = {}
    for folder in os.listdir(MODEL_L0_FOLDER):
        full_path = os.path.join(MODEL_L0_FOLDER, folder)
        try:
            with open(os.path.join(full_path, "metadata.txt"), "r") as metafile:
                metadata = metafile.read()
                metadict[folder] = metadata
        except FileNotFoundError:
            pass

    return metadict

def parse_metadata(metadata):
    """
    Parse metadata string to dictionary for easier variable reading

    Parameters
    ----------
    metadata : STR
        String of metadata as outputted by read_model_metadata.

    Returns
    -------
    parsed : DICT
        Dictionary of metadata.

    """
    parsed = {}
    for line in metadata.split("\n"):
        if "=" in line:
            name, value = line.split("=")
            try:
                parsed[name.strip()] = float(value)
                if name.strip() not in ["pt", "pqm1", "pap"]:
                    parsed[name.strip()] = int(value)
            except:
                parsed[name.strip()] = value.strip()
        else:
            break
            
    return parsed

def copy_model_data(destination_folder):
    """
    Copies model output to the destination folder.

    Parameters
    ----------
    destination_folder : STR
        Name of the folder to which num.dat should be copied.

    Returns
    -------
    None.

    """
    # Start by making sure the input is a string, to forego any funny business
    destination_folder = str(destination_folder)
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)

    # Check if the path already exists - if no, make it
    if not os.path.exists(full_destination_path):
        os.makedirs(full_destination_path)

    # Path exists now, so take the num.dat file from the data/ folder in HAM_box, and copy it to the new folder
    # We're assuming that if the path DOES exist, we're intending to overwrite it. This may come to bite us in the rear
    data_files = [file for file in os.listdir(HAM_DATA_FOLDER) if file.endswith("dat")]
    data_files.remove("rdry_orig.dat")
    data_files.remove("num_m7.dat")

    copy_model_metadata(destination_folder)

    for data_file in data_files:
        shutil.copy(os.path.join(HAM_DATA_FOLDER, data_file), full_destination_path)


def read_model_data(destination_folder):
    """
    Reads the num.dat file within the destination_folder.

    Parameters
    ----------
    destination_folder : STR
        Name of the folder in results/ to be read from.

    Returns
    -------
    num : DataFrame / dictionary of DataFrames
        Pandas DataFrame of the num.dat file. Dictionary of DataFrames when given a directory containing subdirectories.

    metadata : string / dictionary of strings
        String of model metadata. Dictionary of strings when given a directory containing subdirectories.
    """
    # Start by making sure the input is a string, to forego any funny business
    destination_folder = str(destination_folder)
    full_destination_path = os.path.join(MODEL_L0_FOLDER, destination_folder)

    # First check if this folder has subfolders - if yes, it's data from a variation experiment
    has_subfolders = any(os.path.isdir(os.path.join(full_destination_path, item))
                         for item in os.listdir(full_destination_path))

    # Check if the folder exist. If not, boot out immediately. Otherwise, read out the data.
    if not os.path.exists(full_destination_path):
        print("That folder does not exist! Check in the results directory, please.")
        return None, None

    else:
        # If it doesnt have subfolders, just read the data as normal:
        if has_subfolders == False:
            num = pd.read_csv(os.path.join(full_destination_path, "num.dat"),  sep=r"\s+")
            try:
                with open(os.path.join(full_destination_path, "metadata.txt"), "r") as metafile:
                    metadata = metafile.read()

            except FileNotFoundError:
                pass

            return num, metadata

        # Otherwise if it does, read every single subfolder's data, and put every dataframe into a dictionary
        else:
            subfolders = os.listdir(full_destination_path)
            full_subfolder_paths = [os.path.join(full_destination_path, subfolder) for subfolder in subfolders]

            numdict = {}
            metadict = {}
            for subfolder, full_subfolder_path in zip(subfolders, full_subfolder_paths):
                with open(os.path.join(full_subfolder_path, "metadata.txt"), "r") as metafile:
                    metadata = metafile.read()

                metadict[str(subfolder)] = metadata
                numdict[str(subfolder)] = pd.read_csv(os.path.join(full_subfolder_path, "num.dat"), sep=r"\s+")

            return numdict, metadict


def find_keyword(keyword, directory_path=HAM_SRC_FOLDER):
    """This function exists because I am sick and tired of looking through source files myself"""
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

if __name__ == '__main__':
    experiment_name = "test"
    
    environmental_data = [291, 0.0096788, 100901]
    write_environmental_data([environmental_data]*6)
    
    run_model(experiment_name=experiment_name, recompile=True)
    
    bin_boundaries = hp.define_bin_boundaries()
    bin_names = [f"{key}{i+1}" for key, array in bin_boundaries.items() for i, _ in enumerate(array[:-1])]
    num, metadata = read_model_data(experiment_name)
    num5 = pd.read_csv(os.path.join(os.path.join(MODEL_L0_FOLDER, experiment_name), "num_5.dat"), sep = r"\s+")
    metadata = parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # As a sort of blueprint: First check, by metadata, if a model has already been run. If yes, don't run it again
    # but return the data that's already present. If no, go ahead and run it, and copy the data into a new folder.
    fig, axes = hp.plot_size_dist(rdry, num*1e-6, rows=[1, 100, 500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 1)", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    hp.plot_size_dist(rdry, num5*1e-6, rows=[100, 500], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution (cell 5)", populations = ["a"],
                      fig = fig, axes = axes, label = "Model")
    # hp.plot_size_dist_evolution(rdry, num, vmin=1, exp_name = experiment_name, title = "Size distribution evolution")
    hp.stacked_timeseries_plot(num, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 1)")
    hp.stacked_timeseries_plot(num5, populations = ["a"], ymin = 1, exp_name = experiment_name, title = "Size distribution evolution (cell 5)")
    # copy_model_data(experiment_name)
