# =============================================================================
# This file functions to analyse the output of the *SMALL* HAM box model, mostly using the salsa2.0 package, 
# which DOES NOT include atmospheric chemistry and cloud formation.
# There is no GitHub repository or documentation for this distribution.
# Investigation of the source code (driver.f90) seems to indicate that number concentration are
# currently written to output in num.dat.
#
# =============================================================================

import json
import os
from matplotlib import pyplot as plt
# from matplotlib import colors as mc
# from scipy.interpolate import griddata
import re
from itertools import product
import shutil
import subprocess

import netCDF4 as nc
import numpy as np
import pandas as pd
from tqdm import tqdm

import HAM_plot as hp
from utils import lognormal, parse_metadata

HAM_BASE_FOLDER = "../../HAM_box_OpenIFS"
HAM_INPUT_FOLDER = os.path.join(HAM_BASE_FOLDER, "input")
HAM_DATA_FOLDER = os.path.join(HAM_BASE_FOLDER, "data")
HAM_SRC_FOLDER = os.path.join(HAM_BASE_FOLDER, "src\\src_HAM\\")
RESULTS_FOLDER = "../results"
MODEL_L0_FOLDER = os.path.join(RESULTS_FOLDER, "Model L0")
MODEL_PLOT_FOLDER = os.path.join(RESULTS_FOLDER, "Model Figures")

smps_bins_float = [3.0, 11.5, 15.4, 20.5, 27.4, 36.5, 48.7, 64.9, 86.6, 115.5, 154.0, 205.4, 273.8, 365.2]
smps_close = [1669, 1981, 816, 882, 1165, 1360, 1455, 1408, 1069, 507, 11, 0, 0]
smps_far = [474, 634, 431, 711, 1063, 1345, 1479, 1401, 1031, 482, 9, 0, 0]

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
        print(f"Error executing Linux command:\n{e}")


def run_variation(variation_name, nested_environmental_values, particle_flux, dispersion_rate):
    """
    Execute multiple model runs, with varied environental values

    Parameters
    ----------
    variation_name : STRINMG
        Name of the "overarching" experiment.
    nested_environmental_values : List of list
        (Nested) list / matrix of environmental values.
    particle_flux : LIST
        (Nested) list / matrix of particle fluxes.
    dispersion_rate : FLOAT
        (Nested) Dispersion rate per timestep.

    Returns
    -------
    None.

    """
    defaults = [298, 0.0058535, 101325]
    
    # If no flux input is given
    if type(particle_flux) == type(None):
        particle_flux = [1664.5948668000997, 1664.5948668000997, 173524.05010339778,
                         6880411.6196088055, 12947338.526942546, 4834694.831920029,
                         358238.33226261917, 16317.238454753566, 16317.238454753566]
        
    if type(dispersion_rate) == type(None):
        dispersion_rate = 0.01
            
    meta = parse_metadata(read_model_metadata())
    dimension = meta["kproma"]

    # Make the "parent" folder for the data to be written into
    if not os.path.exists(os.path.join(MODEL_L0_FOLDER, variation_name)):
        os.mkdir(os.path.join(MODEL_L0_FOLDER, variation_name))
        
    env_is_nested = any(isinstance(el, list) for el in nested_environmental_values)
    flux_is_nested = any(isinstance(el, list) for el in particle_flux)

    # Environmental variation
    if env_is_nested and not flux_is_nested:
        print("Starting environmental variation")
        write_particle_input_data(particle_flux = particle_flux, dispersion_rate = dispersion_rate)
        
        with tqdm(total=len(nested_environmental_values), desc="Processing", leave=True, position=0) as pbar:
            for run, environmental_values in enumerate(nested_environmental_values):
                if dimension == 1:
                    write_environmental_data(environmental_values)
                else:
                    write_environmental_data([environmental_values]*dimension)
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
            
    
    # Flux/dispersion variation
    elif flux_is_nested and not env_is_nested:
        print("Starting flux/dispersion variation")
        if dimension == 1:
            write_environmental_data(nested_environmental_values)
        else:
            write_environmental_data([nested_environmental_values]*dimension)
        
        # Create a "grid" of particle flux and dispersion rate values
        combinations = list(product(particle_flux, dispersion_rate))
        
        with tqdm(total = len(combinations), desc = "Processing", leave = True, position = 0) as pbar:
            for run, combination in enumerate(combinations):
                flux, dispersion = combination
                write_particle_input_data(particle_flux = flux, dispersion_rate = dispersion)
                
                flux_print = [float(f"{x:.2f}") for x in flux[0:2]]
                description = f"Processing [{flux_print} {dispersion}]"
                description = description.ljust(45, " ")
                pbar.set_description(desc = description)
                experiment_name = os.path.join(variation_name, str(run))
                
                # Recompile on the first (zeroth) run just in case
                if run == 0:
                    run_model(experiment_name, recompile=True, verbose=False)
                else:
                    run_model(experiment_name, recompile=False, verbose=False)
    
                pbar.update(1)
    
    # If both are nested, just throw an error and let the user think about what they've done
    elif flux_is_nested and env_is_nested:
        print("Can't specify variation for both environmental variables and flux/dispersion. Please revise")
        


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

def write_particle_input_data(particle_flux = [4e6, 6e6, 3.6e6, 1.8e5, 3.8e4, 4e3, 0, 0, 0], dispersion_rate = 0.01):
    """
    Writes values to the particle_input.dat file **IN ORDER**
    [input per bin (sep by space) per second]
    [dispersion_rate]

    Parameters
    ----------
    particle_flux : LIST, optional
        Contains particle flux per bin. The default is [4e6, 6e6, 3.6e6, 1.8e5, 3.8e4, 4e3].
    dispersion_rate : FLOAT, optional
        Dispersion rate per timestep. The default is 0.01.

    Returns
    -------
    None.

    """
    if particle_flux is None:
        particle_flux = [4e6, 6e6, 3.6e6, 1.8e5, 3.8e4, 4e3, 0]
    
    content = " ".join(map(str, particle_flux)) + "\n" + str(dispersion_rate)
    
    with open(os.path.join(HAM_INPUT_FOLDER, "particle_input.dat"), "w") as outfile:
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
        content = infile.readlines()

    content = content[0].strip()
    environmental_vals = content.split(" ")

    for var, val in zip(environmental_vars, environmental_vals):
        metadata += f"{var} = {val}\n"

    return metadata

def read_particle_input_data():
    """
    Reads the particle_input.dat file, to see which values have been used to run the model last.

    Returns
    -------
    metadata : STRING
        Formatted string of model environmental data.

    """
    particle_input_names = ["particle_flux", "dispersion_rate"]
    metadata = ""
    with open(os.path.join(HAM_INPUT_FOLDER, "particle_input.dat"), "r") as infile:
        content = infile.read()
        
    particle_input_vals = content.split("\n")
    
    for var, val in zip(particle_input_names, particle_input_vals):
        if var  == "particle_flux":
            val = [float(value) for value in val.split(" ")]
            
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
    metadata += read_particle_input_data()

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


def read_model_data(destination_folder, gridcell = 1):
    """
    Reads the num.dat file within the destination_folder.

    Parameters
    ----------
    destination_folder : STR
        Name of the folder in results/ to be read from.
    gridcell : INT
        Index of gridcells whose data to read. Default is 1.

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
    
    if gridcell == 1:
        num_file = "num.dat" 
    elif gridcell > 1:
        num_file = f"num_{gridcell}.dat"

    # Check if the folder exist. If not, boot out immediately. Otherwise, read out the data.
    if not os.path.exists(full_destination_path):
        print("That folder does not exist! Check in the results directory, please.")
        return None, None

    else:
        # If it doesnt have subfolders, just read the data as normal:
        if has_subfolders == False:
            num = pd.read_csv(os.path.join(full_destination_path, num_file),  sep=r"\s+")
            try:
                with open(os.path.join(full_destination_path, "metadata.txt"), "r") as metafile:
                    metadata = metafile.read()

            except FileNotFoundError:
                pass

            return num, metadata

        # Otherwise if it does, read every single subfolder's data, and put every dataframe into a dictionary
        else:
            subfolders = sorted(os.listdir(full_destination_path), key = float)
            full_subfolder_paths = [os.path.join(full_destination_path, subfolder) for subfolder in subfolders]

            numdict = {}
            metadict = {}
            for subfolder, full_subfolder_path in zip(subfolders, full_subfolder_paths):
                with open(os.path.join(full_subfolder_path, "metadata.txt"), "r") as metafile:
                    metadata = metafile.read()

                metadict[str(subfolder)] = metadata
                numdict[str(subfolder)] = pd.read_csv(os.path.join(full_subfolder_path, num_file), sep=r"\s+")

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
    experiment_name = "gasoline_emissions_with_background_test"
    
    environmental_data = [293, 0.0096610, 99890]
    write_environmental_data([environmental_data]*6)
    
    x = np.linspace(10, 1000, 10000)
    sigma = 1.7
    scale = 8e6
    
    bin_boundaries = hp.define_bin_boundaries()
    salsa_bins = np.unique(np.concatenate(list(bin_boundaries.values()), 0)[0:10] * 1e9)
    
    main_lognormal = lognormal(x, sigma, center_x = 20, scale = scale)
    main_flux = np.interp(salsa_bins, x, main_lognormal)
    
    # secondary_flux = np.interp(salsa_bins, smps_bins_float[1:], smps_far)
    secondary_flux = np.array([0., 0., 0., 1355., 1271., 177., 0., .0, 0.])

    particle_flux = main_flux + 3000*secondary_flux
    # particle_flux = 2*np.array([6e5, 7e5, 7.5e5, 5e6, 7e6, 0])
    dispersion_rate = 0.012
    write_particle_input_data(particle_flux = 2*particle_flux, dispersion_rate = dispersion_rate)
    
    run_model(experiment_name=experiment_name, recompile=True)

    num, metadata = read_model_data(experiment_name, gridcell = 1)
    num5, _ = read_model_data(experiment_name, gridcell = 5)
    metadata = parse_metadata(metadata)
    rdry = pd.read_csv(os.path.join(HAM_DATA_FOLDER, "rdry_orig.dat"), sep=r"\s+")

    # rdry has radii which are off by 2 orders of magnitudes, because SALSA works
    # with cm, "for some reason". Divide everything by 100 to make it SI compliant.
    rdry = rdry/100

    # Plotting
    fig, axes = hp.plot_size_dist(rdry, num, rows=[1], ymin=1, xmin = -20, xmax = 400,
                      exp_name = experiment_name, title = "Size distribution", populations = ["a"],
                      linestyle = "dashed", label = "Model")
    

    fig, axes = hp.plot_size_dist(rdry, num5, rows=[100, 600, 1500], ymin=1, xmin = -20, xmax = 400,
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

    # copy_model_data(experiment_name)
    
    
