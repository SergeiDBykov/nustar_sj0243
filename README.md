#  NuSTAR Pipeline for X-ray pulsar data analysis

This are  scripts and notebooks which are used to run the pipeline for the **NuSTAR** X-ray observatory for the pulsar data. The emphasis  is on  the **phase-averaged** and **phase-resolved** products. 

The data analysis is illustrated on the X-ray pulsar **Swift J0243.6+6124**, and  these results I used in my [paper]().

------

## The structure of the pipeline is as follows:

- `./nustar_data`  
 the data are stored in this folder. The folders are named after the observation ID, e.g. `./90302319004`. Data is downloaded from the [NASA HEASARC archive](https://heasarc.gsfc.nasa.gov).

- `./nustar_products`  
 the products are stored in this folder. The products are organized in folders for each observation. The folders are named after the observation ID with the prefix 'out', e.g. `./out_90302319004`


- `./nustar_scripts`  
 the scripts are stored in this folder. It containt a few scripts to run the pipeline.
    1. `./nustar_scripts/pulsar_init.py`  is used to set up pathe and information about pulsar/observations. See comments in the script for more information.
    2. `./nustar_scripts/nu_class.py` contains functions and classes for the pipeline (e.g. scripts  requesting `nuproducts` ftools commands). 
    3. `./nustar_scripts/nu_pyxspec.py` contains functions to work with `PyXspec`, a python wrapper for `Xspec` spectral fitting package.
    4. `./nustar_scripts/storage.py` contains utility functions to work with X-ray  spectral data.
    5. `./nustar_scripts/nu_utils.py`  contains some utility functions



- `./pipeline_notebooks/`  
 Stores the notebooks to reduce the raw science data into high-level products. The  instructions are given in the notebooks. I highly recommend to create a separate notebooks for each observation. `*_ph_averaged_*` notebooks are used to create the phase-averaged products. `*_ph_resolved_*` notebooks are used to create the phase-resolved products.

- `./data_analysis`  
 the notebooks for the data analysis are stored in this folder.

    1. `./nustar_scripts/model_*/` contains notebooks with spectral fitting of science product. Again I advise to create a separate notebook for each observation, and for every spectral model (e.g. `cutoffpl`) create an appropriate folder (e.g. `./nustar_scripts/model_cuutoffpl/`).
    2. `./nustar_scripts/plots_and_tables/` contains (example) notebooks with the final results report: plots, tables, etc
    
-----

# The pipeline is run as follows:
exampe is given  for observation `90302319004` and spectral model `relxilllp`.


0. Download the data from the HEASARC archive. Exampe is given in the script `./nustar_data/download.sh`.
1. Fill `./nustar_scripts/pulsar_init.py` script to set up the path and information about the observations.
2. Open notebook `./pipeline_notebooks/00_ph_averaged_90302319004.ipynb` to set up and run the pipeline for phase-averaged products. Follow the instructions in the notebook.
3. Open notebook `./pipeline_notebooks/01_ph_resolved_90302319004.ipynb` to set up and run the pipeline for phase-resolved products. Follow the instructions in the notebook.
4. Open notebook `./data_analysis/model_relxilllp/relxilllp_90302319004.ipynb` to set up and run the  spectral models fitting (phase-averaged and phase-resolved spectra). Follow the instructions in the notebook.
5. Use example notebooks from `./data_analysis/plots_and_tables/` to create your own result visualisation (plots, tables, etc).
