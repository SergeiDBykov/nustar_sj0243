import astropy.io.fits as fits
import numpy as np
import pandas as pd
from glob import glob
import os
from scipy.optimize import curve_fit

import matplotlib
from matplotlib import rc
import seaborn as sns
import matplotlib.pyplot as plt
from xspec import AllModels, Fit, Plot,  AllData
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

    'lines.markeredgewidth': 1.5,
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

def set_mpl(palette = 'vaporwave'):
    matplotlib.rcParams.update(rc)
    if palette=='mallsoft':
        sns.set_palette(["#fbcff3", "#f7c0bb", "#acd0f4", "#8690ff", "#30bfdd", "#7fd4c1"], color_codes = True) 
    elif palette=='vaporwave':
        sns.set_palette(['#94D0FF', "#966bff",'#FF6AD5', '#ff6a8b' ,'#8bde8b', '#20de8b'], color_codes = True) 
    else:
        sns.set_palette(palette, color_codes = True)
set_mpl()


def ratio_error(a, b, da, db):
    f = a / b
    sigma = np.abs(f) * np.sqrt((da / a) ** 2 + (db / b) ** 2)
    return f, sigma


# %% SYSTEM WORK


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
    print("Creating command command:", cmd)
    print("Writing to file: ", cmd_name)

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


# def open_dir_in_term():
#     appscript.app("Terminal").do_script(f"cd {os.getcwd()}")


# %% functions


def start_stop(a: np.ndarray, trigger_val: float) -> np.ndarray:
    """
    start_stop finds indexes if beginnings and ends of sequence of trigger values.
    Used in GTI creation (trigger_val = phase bin number)
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


def gauss(t: np.ndarray, t0: float, sigma: float, N: float) -> np.ndarray:
    """
    gaussian function  normalized to N at t0
    """
    return N * np.exp(-((t - t0) ** 2) / (2 * sigma**2))


def fit_efsearch_data(
    efsearcf_fits_file: str, savefig: bool = True, fit_data: bool = True
):
    """
    fit_efsearch_data tries to fit efsearch curve with a gaussian shape

    Args:
        efsearcf_fits_file (str): fits file with efsearch results
        savefig (bool, optional): whether to save figure. Defaults to True.
        fit_data (bool, optional): whether to fit data. Defaults to True.

    Returns:
        p, perr: best fit period and error
    """

    efsearch = fits.open(efsearcf_fits_file)
    period = efsearch[1].data["period"]  # type: ignore
    chisq = efsearch[1].data["chisqrd1"]  # type: ignore
    sigma = (max(period) - min(period)) / 5
    p0 = [period[np.argmax(chisq)], sigma, max(chisq)]

    fig, ax = plt.subplots()
    ax.plot(period, chisq)
    ax.plot(period, gauss(period, *p0), "r-.")

    try:
        popt, perr = curve_fit(gauss, period, chisq, p0=p0)  # type: ignore
        popt = np.array(popt)
        perr = np.array(perr)
        perr = np.sqrt(np.diag(perr))
        ax.plot(period, gauss(period, *popt), "k:")
    except:
        popt = np.array([0])
        perr = np.array([0])

    ax.set_title("Period=" + str(popt[0]) + "  sigma=" + str(popt[1]) + "\n")
    ax.set_xlabel("Period")
    ax.set_ylabel("chi^2")
    if savefig:
        fig.savefig(f"{efsearcf_fits_file}.png")
        plt.close(fig)
    return popt[0], perr[0]


def reduce_list(list: list) -> list:
    flat_list = [item for sublist in list for item in sublist]
    return flat_list


# %% nustar stuff


def PIfromE(E):
    """
    PI channel number from energy
    """
    return (E - 1.6) / 0.04


def EfromPI(PI):
    """
    energy from PI channel number
    """

    return PI * 0.04 + 1.6


def scan_phase_resolved_products(n_phase: int, prodpath: str = "phase_resolved"):
    """
    scan_phase_resolved_products for a given path to phase resolved path, it creates a list of spectra and lightcurves in which out can iterate over every FMP and bin

    Args:
        n_phase (int): number of phase bins
        prodpath (str, optional): product path for phase resolved stuff. Defaults to 'phase_resolved'.

    Returns:
        Tuple: list of light curves (original, bary and orbitally corrected), list of spectra and a function to return list of any combination in form bla-bla-bla_A_binX_sr.format
    """
    binnum = np.arange(1, n_phase + 1)

    def propname(modes, bin, postfix):
        # return [f"phase_resolved_{mode}_bin{bin}_sr.{postfix}" for mode in modes]
        return [f"phase_resolved_bin{bin}{mode}_sr.{postfix}" for mode in modes]

    def propname_per_bin(modes, postfix):
        return [propname(modes, bin, postfix) for bin in binnum]

    lclist = propname_per_bin(["A", "B"], "lc")
    lclist_bary = propname_per_bin(["A", "B"], "lc_bary")
    lclistorb_corr = propname_per_bin(["A", "B"], "lc_bary_orb_corr")

    spelist = propname_per_bin(["A", "B"], "pha")

    return lclist, lclist_bary, lclistorb_corr, spelist, propname_per_bin


def make_grppha_and_wd(
    folder,
    pha_files,
    NCHAN_1="5",
    NCHAN_2="15",
    en_lo="4.",
    en_hi="30.",
    en_split=15,
    labels=None,
    title="",
):
    if folder is not None:
        os.chdir(folder)
    else:
        pass
    enarr = []
    enerrarr = []
    dataarr = []
    dataerrarr = []

    if labels is None:
        labels = pha_files

    for infile in pha_files:
        stemname = infile.split(".")[0:-1]
        stemname = ".".join(stemname)
        outfile = stemname + ".pi_linbin"
        split_chan = int(PIfromE(en_split))
        # if os.path.exists(outfile)  os.remove(outfile) else None
        # grppha = f'''grppha infile="{infile}" outfile="{outfile}"  comm="group 1 4096 {NCHAN} & exit & exit" clobber=yes'''
        grppha = f"""grppha infile="{infile}" outfile="{outfile}"  comm="group 1 {split_chan} {NCHAN_1} {split_chan+1} 4096 {NCHAN_2} & exit" clobber=yes""" #4096 instead of 4095 were here before 
        print(grppha)
        os.system(grppha)

        from xspec import Plot
        from xspec import AllData

        AllData.clear()
        if folder is not None:
            AllData(f"1:1 {outfile}")
        else:
            path, filename = os.path.split(infile)
            stemname = filename.split(".")[0]
            filename = stemname + ".pi_linbin"

            print(infile, path, filename)
            os.chdir(path)
            AllData(f"1:1 {filename}")
        AllData.ignore(f"**-{en_lo} {en_hi}-**")
        # AllData.ignore("bad")
        Plot.device = "/null"
        Plot("da")
        Plot.xAxis = "keV"
        en = Plot.x()
        data = Plot.y()
        en_err = Plot.xErr()
        data_err = Plot.yErr()

        enarr.append(np.array(en))
        dataarr.append(np.array(data))
        dataerrarr.append(np.array(data_err))
        enerrarr.append(np.array(en_err))

        os.remove(outfile)

    enarr = np.array(enarr)
    dataarr = np.array(dataarr)
    dataerrarr = np.array(dataerrarr)
    enerrarr = np.array(enerrarr)
    # labels[0] = labels[0]+' (ref)'

    fig, [ax_spe, ax_rat] = plt.subplots(2, figsize=(8, 8), sharex=True)  # type: ignore
    plt.subplots_adjust(hspace=0, wspace=0)
    color = "red"
    ax_spe.plot(enarr[0], dataarr[0], "o", label=labels[0], alpha=0.7, color=color)
    ax_spe.errorbar(
        enarr[0], dataarr[0], dataerrarr[0], fmt="none", alpha=0.4, ecolor=color
    )

    for i, (d, e) in enumerate(zip(dataarr[1:], dataerrarr[1:]), 1):
        assert len(d) == len(dataarr[0]), "different size of data"
        # d = np.array(d)
        # e = np.array(e)
        # dataarr[0] = np.array(dataarr[0])
        # dataerrarr[0] = np.array(dataerrarr[0])
        try:
            rat, rat_err = ratio_error(np.array(d), np.array(dataarr[0]), np.array(e), np.array(dataerrarr[0]))
        except:
            raise Exception(f"ratio error for {type(d)} and {type(dataarr[0])}")
        norm_factor = 1 / rat[0] * (1 + i * 0.2)
        rat = rat * norm_factor
        rat_err = rat_err * norm_factor

        ax_rat.step(
            enarr[0],
            rat,
            "o-",
            label=f"{labels[i]}/{labels[0]}  x  {norm_factor:=.2f}",
            where="mid",
            alpha=0.7,
        )
        color = ax_rat.get_lines()[-1].get_color()
        ax_rat.errorbar(
            enarr[0], rat, rat_err, enerrarr[0], ecolor=color, fmt="none", alpha=0.4
        )

        ax_spe.plot(enarr[0], d, "o", label=labels[i], alpha=0.7, color=color)
        ax_spe.errorbar(enarr[0], d, e, fmt="none", alpha=0.4, ecolor=color)

    for axt in [ax_rat, ax_spe]:  # , ax_zoom]:
        axt.set_xscale("log")
        axt.axvline(6.4, color="k", ls="--", lw=0.5)
        axt.axvline(6.7, color="k", ls="--", lw=0.5)
        axt.axvline(7.1, color="k", ls="--", lw=0.5)
        axt.set_xlabel("Energy (keV)")

    ax_spe.set_title(title)
    ax_spe.set_yscale("log")
    ax_spe.legend(loc="lower left")
    ax_rat.legend(loc="upper right")
    ax_spe.set_ylabel("Photons/ (keV s cm^2)")
    ax_rat.set_ylabel("Ratio to ref, arb. units")
    # ax.legend(bbox_to_anchor=(1, 0.5), loc='center left', fontsize=8)
    return fig



# def plot_spe_ratio(
#     model,
#     ph_res_folder,
#     bins_number,
#     labels=None,
#     title="",
#     min_sig = 75,
#     min_bin = 75,
#     zoom_en = [4, 30],
#     zoom_rat = [0.8, 1.3],
#     colors = None,
#     ax_rat = None,
# ):
#     if ph_res_folder is not None:
#         os.chdir(ph_res_folder)
#     else:
#         pass

#     #pha_files = glob(f"{ph_res_folder}/phase_resolved_bin{bins_number}A_sr.pi")
#     pha_files = [f"{ph_res_folder}/phase_resolved_bin{bin}A_sr.pi" for bin in bins_number]
#     pha_files = pha_files + [f"{ph_res_folder}/phase_resolved_bin{bin}B_sr.pi" for bin in bins_number]
#     pha_files = [x.rsplit('/', 1)[1] for x in pha_files]
#     if labels is None:
#         #labels = sorted([f'Bin {i}' for i in bins_number])
#         labels = [f'Bin {i}' for i in bins_number]
#     print(pha_files)

#     enarr = []
#     dataarr = []
#     dataerrarr = []
#     enerrarr = []
#     for infile in pha_files:
#         from xspec import Plot
#         from xspec import AllData

#         AllData.clear()
#         #AllData(' '.join([f'{i}:{i} {file}' for i,file in enumerate(spelist, 1)]))
#         AllData(f"1:1 {infile}")
#         AllData.ignore(f"*: **-4. 79.-**")
#         AllData.ignore("bad")
#         m1 = model
#         #m1.setPars({1: '1. -1'})
#         #m2 = AllModels(2)
#         #m2.setPars({1: '1. '})

#         Fit.query = 'yes'
#         Fit.statMethod = "chi"
#         Fit.query = "yes"
#         Fit.perform()
#         print('fitting done')

#         Plot.device = "/null"
#         Plot.setRebin(min_sig, min_bin)
#         Plot("rat")
#         Plot.xAxis = "keV"
#         en = Plot.x()
#         data = Plot.y()
#         en_err = Plot.xErr()
#         data_err = Plot.yErr()

#         enarr.append(np.array(en))
#         dataarr.append(np.array(data))
#         dataerrarr.append(np.array(data_err))
#         enerrarr.append(np.array(en_err))


#     enarr = enarr
#     dataarr = dataarr
#     dataerrarr = dataerrarr
#     enerrarr = enerrarr
#     # labels[0] = labels[0]+' (ref)'
#     ms, alpha = 7, 0.7
#     if ax_rat is None:
#         fig, ax_rat = plt.subplots( figsize=(12, 8), sharex=True)  # type: ignore
#         plt.subplots_adjust(hspace=0, wspace=0)
#     else:
#         pass
#     for i,(en, en_err, data, data_err) in enumerate(zip(enarr, enerrarr, dataarr, dataerrarr)):
#         #norm_factor = 1 / data[0] #* (1 + i * 0.2)
#         norm_factor = 1
#         data = data * norm_factor
#         data_err = data_err * norm_factor
#         if colors is None:
#             ax_rat.plot(en, data, 'o', label=None, alpha = alpha, ms = ms)   
#             color =  ax_rat.get_lines()[-1].get_color()         
#         #if colors is None:
#         #    color = plt.rcParams['axes.prop_cycle'].by_key()['color'][##(bins_number[i])%6]
#         #else:
#         #    color = colors[bins_number[i]]
        
#         #ax_rat.plot(en, data, 'o', label=labels[i], color=color, alpha = alpha, ms = ms)

#         ax_rat.errorbar(en, data, data_err, en_err,
#                         fmt='none', ecolor=color, alpha=alpha)
#         #ax_rat.fill_between(en, data - data_err, data+data_err, alpha=0.5, #facecolor = color, edgecolor = 'k', label=labels[i])


#     for axt in [ax_rat]:  # , ax_zoom]:
#         axt.set_xscale("log")
#         axt.axvline(6.4, color="k", ls="--", lw=0.5)
#         axt.axvline(6.7, color="k", ls="--", lw=0.5)
#         axt.axvline(7.1, color="k", ls="--", lw=0.5)
#         axt.set_xlabel("Energy (keV)")
#     #ax_rat.set_xlim(zoom_en)
#     #ax_rat.set_ylim(zoom_rat)

#     ax_rat.set_title(title)
#     ax_rat.legend(loc="upper right")
#     ax_rat.set_ylabel("Ratio to ref, arb. units")
#     # ax.legend(bbox_to_anchor=(1, 0.5), loc='center left', fontsize=8)
#     return ax_rat



# bins = ['2', '10', '4', '7']
# make_grppha_and_wd(folder='/Users/sdbykov/work/xray_pulsars/sj0243_nu/results/out90302319004/products/phase_resolved_shift_1/', pha_files=[
#     f'phase_resolved_bin{bin}A_sr.pha' for bin in bins], NCHAN_1=5, NCHAN_2=25, en_split=10., labels=[f'bin {bin}' for bin in bins], title='90302319004')

# make_grppha_and_wd(folder='/Users/sdbykov/work/xray_pulsars/sj0243_nu/results/out90302319004/products/phase_resolved_shift_1/', pha_files=[
#                     'phase_resolved_bin2B_sr.pha', 'phase_resolved_bin10B_sr.pha', 'phase_resolved_bin4B_sr.pha', 'phase_resolved_bin7B_sr.pha', 'phase_resolved_bin3B_sr.pha'], NCHAN_1=5, NCHAN_2=25, en_split=10.)


# make_grppha_and_wd(folder='/Users/sdbykov/work/xray_pulsars/sj0243_nu/results/out90302319006/products/phase_resolved/', pha_files=[
#                    'phase_resolved_bin1A_sr.pha', 'phase_resolved_bin2A_sr.pha', 'phase_resolved_bin3A_sr.pha', 'phase_resolved_bin4A_sr.pha',  'phase_resolved_bin5A_sr.pha', 'phase_resolved_bin7A_sr.pha', 'phase_resolved_bin9A_sr.pha'], NCHAN=5)


# bins = ['9', '2', '3', '4', '7', '1']
# make_grppha_and_wd(folder='/Users/sdbykov/work/xray_pulsars/sj0243_nu/results/out90401334002/products/phase_resolved/', pha_files=[
#     f'phase_resolved_bin{bin}A_sr.pha' for bin in bins], NCHAN_1=10, NCHAN_2=40, en_split=10., en_hi='50.', labels=[f'bin {bin}' for bin in bins], title='90401334002')

# make_grppha_and_wd(folder='/Users/sdbykov/work/xray_pulsars/sj0243_nu/results/out90401334002/products/phase_resolved', pha_files=['phase_resolved_bin9A_sr.pha', 'phase_resolved_bin2A_sr.pha', 'phase_resolved_bin3A_sr.pha', 'phase_resolved_bin4A_sr.pha', 'phase_resolved_bin7A_sr.pha', 'phase_resolved_bin1A_sr.pha'], NCHAN='40', en_hi = '79.')
