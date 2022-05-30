from typing import Tuple
import astropy.io.fits as fits
from scipy.optimize import curve_fit
import numpy as np
import os
import pandas as pd
import matplotlib
import seaborn as sns
import matplotlib.pyplot as plt



### matplitlib settings

# matplotlib.use('MacOSX') 
rc = {
    "figure.figsize": [10, 10],
    "figure.dpi": 100,
    "savefig.dpi": 300,
    # fonts and text sizes
    #'font.family': 'sans-serif',
    #'font.family': 'Calibri',
    #'font.sans-serif': 'Lucida Grande',
    'font.style': 'normal',
    "font.size": 15,
    "axes.labelsize": 15,
    "axes.titlesize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,

    # lines
    "axes.linewidth": 1.25,
    "lines.linewidth": 1.75,
    "patch.linewidth": 1,

    # grid
    "axes.grid": True,
    "axes.grid.which": "major",
    "grid.linestyle": "--",
    "grid.linewidth": 0.75,
    "grid.alpha": 0.75,

    # ticks
    "xtick.top": True,
    "ytick.right": True,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.width": 1.25,
    "ytick.major.width": 1.25,
    "xtick.minor.width": 1,
    "ytick.minor.width": 1,
    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 4,
    "ytick.minor.size": 4,

    'lines.markeredgewidth': 0, #1.5,
    "lines.markersize": 10,
    "lines.markeredgecolor": "k",
    'axes.titlelocation': 'left',
    "axes.formatter.limits": [-3, 3],
    "axes.formatter.use_mathtext": True,
    "axes.formatter.min_exponent": 2,
    'axes.formatter.useoffset': False,
    "figure.autolayout": False,
    "hist.bins": "auto",
    "scatter.edgecolors": "k",
}

def set_mpl(palette = 'shap', desat = 0.8):
    matplotlib.rcParams.update(rc)
    if palette == 'shap':
        #colors from shap package: https://github.com/slundberg/shap
        cp = sns.color_palette( ["#1E88E5", "#ff0d57", "#13B755", "#7C52FF", "#FFC000", "#00AEEF"])
        sns.set_palette(cp, color_codes = True, desat = desat)
    elif palette == 'shap_paired':
        #colors from shap package: https://github.com/slundberg/shap, + my own pairing of colors
        cp = sns.color_palette( ["#1E88E5", "#1e25e5", "#ff0d57", "#ff5a8c",  "#13B755", "#2de979","#7C52FF", "#b69fff", "#FFC000", "#ffd34d","#00AEEF", '#3dcaff'])
        sns.set_palette(cp, color_codes = True, desat = desat)
    else:
        sns.set_palette(palette, color_codes = True, desat = desat)
set_mpl()


def ratio_error(a, b, da, db):
    #calc  the error on ratio of two variables a and b with their errors da and db
    f = a / b
    sigma = np.abs(f) * np.sqrt((da / a) ** 2 + (db / b) ** 2)
    return f, sigma


### OS commands


def run_command(cmd: str, cmd_name: str, rewrite: bool = True) -> str:
    """
    run_command creates a executable file (.sh) with the command. It DOES NOT run this command. This shoud be run in terminal.

    Args:
        cmd (str): command to execute
        cmd_name (str): name of the .sh file
        rewrite (bool, optional): whether to delete existing file with the same name before writing. Defaults to True.

    Returns:
        str: path to command
    """
    #print("Creating command command:", cmd)
    print("Writing to file: ", cmd_name+'.sh')

    if rewrite:
        os.system(f"rm -f {cmd_name}.sh")

    os.system(f"echo '\n {cmd}' >> {cmd_name}.sh")
    os.system(f"chmod +x {cmd_name}.sh")
    return os.path.abspath(cmd_name)


def create_dir(dir: str):
    """
    create_dir creates a directory with given name. If exists, does nothing

    Args:
        dir (str): directory to create
    """
    os.system(f"mkdir -p {dir}")


### useful functions for the pipeline


def start_stop(a: np.ndarray, trigger_val: float) -> np.ndarray:
    """
    finds indexes ff beginnings and ends of sequence of trigger values.
    Used in GTI creation in phase-resolved spectroscopy (trigger_val = phase bin number)
    source https://stackoverflow.com/questions/50465162/numpy-find-indeces-of-mask-edges

    Args:
        a (np.ndarray): input array
        trigger_val (float): value to trigger masking. Defaults to 1.

    Returns:
        np.ndarray: indeces array
    """
    # "Enclose" mask with sentients to catch shifts later on
    mask = np.r_[False, a == trigger_val, False]
    # Get the shifting indices
    idx = np.flatnonzero(mask[1:] != mask[:-1])

    return idx.reshape(-1, 2) - [0, 1]



def reduce_list(list: list) -> list:
    flat_list = [item for sublist in list for item in sublist]
    return flat_list

