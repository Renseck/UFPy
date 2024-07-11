<h1 align = "center"> UFPy </h1> <br>

<p align = "center">
 <img src = "header.png" alt= "Image header"/>
</p>

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Acknowledgments](#acknowledgments)
- [License](#license)

## Introduction

Master thesis project, in which we try to model ultrafine particle (UFP) concentrations using the [HAM](https://redmine.hammoz.ethz.ch/projects/hammoz) microphyics package, along with the [SALSA2.0](https://gmd.copernicus.org/articles/11/3833/2018/) module. 
Model data will be compared to measurement data. Measurements were made by [RIVM](https://www.rivm.nl/). Results of this may be used to improve performance of the LOTOS-EUROS model, used by the RIVM to forecast air quality. **NOTE**: this model works in conjuction with a very specific distribution of the HAM/SALSA2.0 box model, which is not freely sharable.

## Installation

Start by downloading the UFPy folder, and make sure you place the HAM_box_OpenIFS (the folder must have this exact name, not freely available) model into the folder which contains the UFPy folder. Next, create a conda environment using the file supplied in the `src/` folder, by running `conda env create -f path/to/UFPy/src/environment.yml`, from wherever you want to install this environment. Then activate this environment by running `conda activate HAMSALSA`. To run this model in a Windows environment, one first needs to install a Linux Subsystem, such as WSL (which was used by me; installation instructions can be found [here](https://learn.microsoft.com/en-us/windows/wsl/install)). Running the model requires packages like `netcdf`, `hdf4`, `hdf5` to be installed first. I haven’t kept a neat list of what those packages are, but they will come up as you’re trying to `make` the model - install as needed.

Crucially, **glibc** ≥ 2.29 needs to be installed in your Ubuntu, or it doesn’t seem to work at all. This can be tricky. Updating the current glibc version may very well break your OS, but there is a [way](https://unix.stackexchange.com/a/299665) to install an *additional* glibc, which can then be linked to. This overcame any issues for me, and hasn’t yet presented any new ones (02/02/2024).

## Usage

See the documentation files in `src/` for some example runs, as well as detailed descriptions of the functions in the py-files. A simple `python main.py` will run the model, provided you've activated the correct environment.

## Acknowledgments

Supervised by [Prof. Guus Velders](https://www.uu.nl/medewerkers/GJMVelders).

## License

[MIT](https://choosealicense.com/licenses/mit/) license, Copyright (c) 2024 Rens van Eck.