import pandas as pd
from glob import glob
import seaborn as sns
from matplotlib import rc
import matplotlib
import numpy as np
import pickle
import matplotlib.pyplot as plt
from xspec import Spectrum, Model, Fit, Plot, AllData, AllModels, Xset
import xspec
import os
Xset.parallel.error = 20
Xset.chatter = 5
Fit.query = "yes"

def create_dir(dir: str):
    """
    create_dir creates a directory with given name. If exists, does nothing

    Args:
        dir (str): directory to create
    """
    os.system(f'mkdir -p {dir}')


def make_figure(log=1, residuals='rat'):
    fig = plt.figure(figsize=(14, 6))
    rows = 3
    cols = 3
    ax_spe = plt.subplot2grid((rows, cols), (0, 0), rowspan=2, colspan=3)
    ax_ra = plt.subplot2grid((rows, cols), (2, 0),
                             rowspan=1, colspan=3, sharex=ax_spe)
    plt.subplots_adjust(hspace=0)
    ax_ra.set_xlabel('E, keV')  # type: ignore
    if residuals == 'rat':
        ax_ra.axhline(1, color='k')  # type: ignore
        ax_ra.set_ylabel('$data/model$')  # type: ignore
    if residuals == 'del':
        ax_ra.axhline(0, color='k')  # type: ignore
        ax_ra.set_ylabel('$[data-model]/error$')  # type: ignore

    ax_spe.set_ylabel(  # type: ignore
        '$EF_E, keV^2 (phot\, cm^{-2} s^{-1} keV^{-1})$')  # type: ignore
    if log:
        for ax in [ax_spe]:
            ax.set_xscale('log')  # type: ignore
            ax.set_yscale('log')  # type: ignore
        ax_ra.set_yscale('linear')  # type: ignore

    return fig, ax_spe, ax_ra


def showmodel(m):
    # stolen from https://github.com/evandromr/pyxspec_utils
    """Print current model information
        Display a formated view of current model information
        such as the one produced by `model.show()` on pyXspec
        or by `show par` on Xspec.
        The errors are taken from the `xspec.Fit.error` calculation
        if that was not performed errros will be zero.
        Parameters
        ----------
        m: Xspec.Model
            The model from which you want information
        Returns
        --------
        Print output as a formated table
        Example
        --------
        >>> import xspec
        >>> import pyxspec_utils as pu
        >>> m1 = xspec.Models("wabs(powerlaw+mekal)")
        >>> pu.printmodel(m1)
        Model: wabs(powerlaw + mekal)
        P#   C#   Component    Parameter  Unit    Value        Errors
        -------------------------------------------------------------
        1    1    wabs         nH         10^22   1.0     (0.0 , 0.0)
        2    2    powerlaw     PhoIndex           1.0     (0.0 , 0.0)
        3    2    powerlaw     norm               1.0     (0.0 , 0.0)
        4    3    mekal        kT         keV     1.0     (0.0 , 0.0)
        5    3    mekal        nH         cm-3    1.0     (0.0 , 0.0)
        6    3    mekal        Abundanc           1.0     (0.0 , 0.0)
        7    3    mekal        Redshift           0.0     (0.0 , 0.0)
        8    3    mekal        switch             1.0     (0.0 , 0.0)
        9    3    mekal        norm               1.0     (0.0 , 0.0)
    """
    print("Model: {}".format(m.expression))
    print()
    print("{:4} {:4} {:12} {:10} {:7} {:15} {:12}".format("P#",
                                                          "C#",
                                                          "Component",
                                                          "Parameter",
                                                          "Unit",
                                                          "Value",
                                                          "Errors"))
    print("--"*38)
    pid = 1
    for cid, component in enumerate(m.componentNames):
        for parameter in eval("m.{}.parameterNames".format(component)):
            u = eval("m.{}.{}.unit".format(component, parameter))
            val = eval("m.{}.{}.values[0]".format(component, parameter))
            err = eval("m.{}.{}.error[:2]".format(component, parameter))
            print("{:<4} {:<4} {:<12} {:<10} {:<7} {:<10.5} \
                   ({:<10.5}, {:<10.5})".format(pid, cid + 1, component,
                                                parameter, u, val,
                                                err[0], err[1]))
            pid += 1


'''
#in case you need components name
for comp_name in m904.componentNames:
    comp = getattr(m904, comp_name)
    for par_name in comp.parameterNames:
        par = getattr(comp, par_name)
        par.frozen = True
        print(comp.name, par.name)
'''
