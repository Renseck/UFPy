
<h1 align = "center"> UFPy </h1> <br>

This file will serve as a complete and thorough documentation of the (main) functions contained in this project, their arguments and how they work together.

1. [run_model](#run_model)
2. [Model metadata](#model-metadata)
3. [Model data](#model-data)
4. [Plotting](#plotting)

---

## run_model
This commands runs the model. Make sure that the HAM_box_OpenIFS model is also contained within your GitHub folder (needed for relative pathing).

```python
def run_model(experiment_name, recompile = True, verbose = True):
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
```
Model is run (and recompiled first if set to `True`) by sending a command through the WSL terminal. Output  data is copied and saved to a folder named `experiment_name` in the `results/` folder (see [copy_model_data](#copy_model_data)).

---

## Model metadata

### read_model_metadata
```python
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
```
No arguments. Uses [gen_params](#gen_params) to read the Fortran source files for the parameters set in `params.json`, and extracts their values (through regex). These are all written into a string, with each line containing one parameter and its value. Then, [read_enviromental_data](#read_environmental_data) is called, and added to the string of metadata. Primarily built to be used in [copy_model_metadata](#copy_model_metadata).

### copy_model_metadata
```python
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
```
Set a value to `desination_folder` to copy the model metadata of that run to a folder in `results` of that name. This function calls [read_model_metadata](#read_model_metadata), and writes it into the `metadata.txt` file, to keep track of the settings of each model run individually.

### gen_params
```python
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
```
No arguments. Reads the params.json file, which is used to document which parameters of the model need to be kept track of. Each of these requires:
- A name, which should be identical to how it is found in the `.f90` source files of the model,
- An indication whether they are numerical (`num`) or boolean (`bool`),
- The name of the file in which the parameter is found (without .`f90`).

### read_environmental_data
```python
def read_environmental_data():
    """
    Reads the environmental.dat file, to see which values have been used to run the model last.
    
    Returns
    -------
    metadata : STRING
        Formatted string of model environmental data.

    """
```
No arguments. Reads the `environmental.dat` file, which contains, in order, the values to be fed into the models input for ambient temperature, specific humidity and ambient pressure (resp. pt, pqm1, pap). 

### write_environmental_data
```python
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
```
Arugment `environmental_vals` should be given as `list`, containing the values of the model environmental variables in the following order:
\[Ambient temperature, specific humidity, ambient pressure\], or in the terms used in the `.f90` source files: \[pt, pqm1, pap\].
These are written to the `environmental.dat` file in the `input` folder of the HAM_box_OpenIFS folder, to be read by the model. Instructions on how to edit the model source files may follow. This function is primarily built to modify environmental values without having to recompile the model, saving time.

### check_metadata
```python
def check_metadata():
    """
    Loops through all model result folders and reads their metadata.
    
    Returns
    -------
    metadict : DICT
        Dictionary, with keys the name of every result folder, and values being the parameters of each model run folder.

    """
```
No arguments. Collects all the metadata available from each model run so far, and outputs them into a dictionary. Can be used to reduce running time dynamically, by checking whether a run has already been performed, and skipping a certain set of settings if so. Not yet in use.

---

## Model data

### copy_model_data
```python
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
```
Set a value to `desination_folder` to copy the model output of that run to a folder in `results` of that name. This function also calls [copy_model_metadata](#copy_model_metadata), and writes this data into the same folder.

### read_model_data
```python
def read_model_data(destination_folder):
    """
    Reads the num.dat file within the destination_folder.
    
    Parameters
    ----------
    destination_folder : STR
        Name of the folder in results/ to be read from.

    Returns
    -------
    num : DataFrame
        Pandas DataFrame of the num.dat file.

    """
```
Primary way of reading back the output of model runs. Avoids having to rerun the model over and over again, at the expense of some storage space.

---

## Plotting

These functions were supplied graciously by [Harri Kokkola](#https://en.ilmatieteenlaitos.fi/cv-harri-kokkola). They have been edited minutely to allow for dynamic saving.

### define_bin_boundaries
```python
def define_bin_boundaries(populations = ['1a', '2a', '2b']):
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
```

### plot_size_dist
```python
def plot_size_dist(
    rdry, num, rows = [0], populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    exp_name = "", title = "",
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
        Top y-axis likmit. The default is None.
    exp_name : String, optional
        Name of the experiment. The default is "". Leave empty to forego saving the image.
    title : String, optional
        Title of the image. The default is "".

    Returns
    -------
    None.

    """
```

### plot_size_dist_evolution
```python
def plot_size_dist_evolution(
    rdry, num, populations = ['a', 'b'],
    xmin = None, xmax = None,
    ymin = None, ymax = None,
    vmin = None, vmax = None,
    exp_name = "", title = "",):
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

    Returns
    -------
    None.

    """
```