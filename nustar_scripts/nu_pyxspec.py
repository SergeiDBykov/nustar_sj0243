from typing import Optional, Tuple

from nustar_scripts.nu_class import NustarObservation
from .nu_utils import set_mpl, create_dir
set_mpl()
import os
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
from matplotlib import gridspec
import xspec
from xspec import Model, Fit, Plot, AllData, AllModels,  Xset
from .storage import Container, Storage
from typing import Optional



### (Py)Xspec settings ###
AllData.clear()
Xset.parallel.error = 20
Xset.chatter = 5
Fit.query = "yes"
Plot.device = '/null'


### displaying functions ###

def make_figure():
    """Make a figure with a grid of subplots for plotting spectra in eeufs format"""

    fig, [ax_spe, ax_ra] =  plt.subplots(nrows=2, ncols = 1, sharex = True, gridspec_kw = {'hspace':0, 'height_ratios': [2,1]}, figsize = (14,6)) # type: ignore

    ax_ra.set_xlabel('E, keV')  # type: ignore

    ax_ra.axhline(0, color='k')  # type: ignore
    ax_ra.set_ylabel('$[data-model]/error$')  # type: ignore

    ax_spe.set_ylabel(  # type: ignore
        '$EF_E, keV^2 (phot cm^{-2} s^{-1} keV^{-1})$')  # for eeufs plots type: ignore

    for ax in [ax_spe]:
        ax.set_xscale('log')  # type: ignore
        ax.set_yscale('log')  # type: ignore
    ax_ra.set_yscale('linear')  # type: ignore

    return fig, ax_spe, ax_ra


def showmodel(m):
    # stolen blatantly from https://github.com/evandromr/pyxspec_utils
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


def plot_spe_nu(Plot,
    ax_spe: Optional[plt.Axes]=None,
    ax_res: Optional[plt.Axes]=None,
    fig: Optional[plt.Figure] = None,  # type: ignore
    label:  str='',
    min_sig: str='50', min_bin: str='50',
    color1:  str = 'C0', color2:  str = 'C1',
    plotmode:  str='eeufs',
    figname=None):
    """Plot spectra in eeufs format for NuSTAR data (i.e. two detectors)"""


    ms = 5
    alpha = 0.4
    if ax_spe is None:
        fig, ax_spe, ax_res = make_figure()
    colors = [color1, color2]
    for group in [1, 2]:
        Plot.device = '/null'
        Plot(plotmode)
        Plot.setRebin(min_sig, min_bin)
        Plot.xAxis = 'keV'

        en = Plot.x(plotGroup=group)
        data = Plot.y(plotGroup=group)
        model = Plot.model(plotGroup=group)

        en_err = Plot.xErr(plotGroup=group)
        data_err = Plot.yErr(plotGroup=group)

        ax_spe.plot(en, data, 'o', label=label, color=colors[group-1], alpha = alpha, ms = ms)

        color = ax_spe.get_lines()[-1].get_color()
        ax_spe.plot(en, model, ':', color=color, alpha=alpha)
        ax_spe.errorbar(en, data, data_err, en_err,
                        fmt='none', ecolor=color, alpha=alpha)

        Plot.device = '/null'

        Plot('del')
        Plot.setRebin(min_sig, min_bin)
        Plot.xAxis = 'keV'

        rat = Plot.y(plotGroup=group)
        rat_err = Plot.yErr(plotGroup=group)

        ax_res.plot(en, rat, 'o', label=label, color=color, alpha=alpha, ms = ms)
        ax_res.errorbar(en, rat, rat_err, en_err,
                        fmt='none', ecolor=color, alpha=alpha)

        if figname is not None:
            fig.savefig(figname)

    return None

def load_nu_spe(spe1: str, spe2: str, en_lo='4.', en_hi='79.'):
    AllData.clear()
    AllData(f"1:1 {spe1} 2:2 {spe2}")
    AllData.ignore(f'**-{en_lo} {en_hi}-**')
    AllData.ignore("bad")
    s1, s2 = AllData(1), AllData(2)
    print(f' loaded {spe1} and {spe2} from {os.getcwd()}')
    return s1, s2


def load_bin_spe(bin_num: str, en_lo='4.', en_hi='79.'):
    spe1 = f'phase_resolved_bin{bin_num}A_sr.pi'
    spe2 = f'phase_resolved_bin{bin_num}B_sr.pi'
    return load_nu_spe(spe1=spe1, spe2=spe2, en_lo=en_lo, en_hi=en_hi)



def fit_spectra(model: xspec.model.Model,
                prefix: str,  
                model_name: str,
                dataset: str = 'spe_and_lc',
                error_conf: str = '2.71',
                en_lo: str='4.', en_hi: str='79.',
                min_sig: str = '50', min_bin: str ='50',
                rewrite: bool = True,
                calc_errors: bool=True,
                ignore_comp_errors:Optional[list] = None,
                eqw_comps: Optional[list] = None,
                eqw_cl: str='90') -> Storage:
    """
    fit_spectra loads given spectra and fits it with the given model

    Args:
        model (xspec.model.Model): xspec model to fit spectra with. If string, the model is  loaded from the model file (if any).
        prefix (str): prefix for the model name, e.g. Observation number
        model_name (str): name of a model, e.g. 'cutoffplaw'
        dataset (str, optional): spectral files to load. adds _srA/_srB to it.  Defaults to 'spe_and_lc'. Function should be called in the folder with the spectra and responses.
        error_conf (str, optional): Xspec error command. 2.71 is 90% probability interval. Defaults to '2.71'.
        en_lo (str, optional): lower energy of the spectral fits, in keV. Defaults to '4.'.
        en_hi (str, optional): upper energy of the spectral fits, in keV. Defaults to '79.'.
        min_sig (str, optional): setpl  for plotting: minimum significance of the bin. Defaults to '50'.
        min_bin (str, optional): setpl for plotting: maximum number of adjacent bins . Defaults to '50'.
        rewrite (bool, optional): whether to rewrite output files. Defaults to True.
        calc_errors (bool, optional): whether to call error command in Xspec. Defaults to True.
        ignore_comp_errors (Optional[list], optional): List of components for which errors are not calculated. If none, calculates every error. Defaults to None.
        eqw_comps (Optional[list], optional): List of components for which equivalent width is calculated. If none, calculates nothing. Defaults to None.
        eqw_cl (str, optional): confidence interlval for eqw calculation, in per cents. Defaults to '90'.

    Returns:
        fit Storage  object
    """


    Fit.query = "yes"
    assert '_' not in model_name,  'model_name cannot contain "_"' 
    data_name = prefix + '_'+model_name
    create_dir(f'xspec')
    create_dir(f'xspec/{model_name}')
    create_dir(f'xspec/{model_name}/xcm')
    if rewrite:
        os.system(f'rm -rf xspec/{model_name}/*{data_name}*')
        os.system(f'rm -rf xspec/{model_name}/xcm/*{data_name}*')
        print(
            f"deleted xspec/{model_name}/*{data_name}* and xspec/{model_name}/xcm/*{data_name}* files")

    try:
        s = Storage().from_pikle(f'xspec/{model_name}/{data_name}.storage')
        print(f'loading storage xspec/{model_name}/{data_name}.storage...')
        return s
    except:
        print('no storage found. fitting...')
        pass

    spe1 = f'{dataset}A_sr.pi'
    spe2 = f'{dataset}B_sr.pi'
    S1, S2 = load_nu_spe(spe1=spe1, spe2=spe2, en_lo=en_lo, en_hi=en_hi)

    #save xspec data infromation for quick loading via xspec's @<path to xcm file> command.
    xspec.Xset.save(f'xspec/{model_name}/xcm/{data_name}.xcm_data', info='f')

    if isinstance(model, xspec.model.Model):
        m1 = model
        # untie cross-calibration constant
        m1.setPars({1: '1. -1'})
        m2 = AllModels(2)
        m2.setPars({1: '1. '})
    else:
        try:
            xspec.Xset.restore(model)
            m1 = AllModels(1)
            m2 = AllModels(2)
            model = AllModels(1)
            print('model restored from file')

        except:
            raise Exception(
                f'model not understood. should be either xspec.model.Model or xspec .xcm file with model parameters. \n Given: {model}')

    #save initial model for quick access
    xspec.Xset.save(
        f'xspec/{model_name}/xcm/{data_name}.xcm_model_init', info='m')

    Fit.statMethod = "chi"

    Fit.query = 'yes'
    Fit.perform()
    Fit.query = 'yes'
    print('fitting done')
    if calc_errors:
        if ignore_comp_errors is None:
            ignore_comp_errors = []
        print(f'skipping errors for: {ignore_comp_errors}')
        err_idx = [str(m2.startParIndex)]
        for comp_name in m1.componentNames:
            comp = getattr(m1, comp_name)
            for par_name in comp.parameterNames:
                par = getattr(comp, par_name)
                if not par.frozen and par not in ignore_comp_errors:
                    err_idx.append(str(par.index))
        for err_ix in err_idx:
            try:
                err_commdns = f"nonew maximum 20. {error_conf} {err_ix}"
                Fit.error(err_commdns)
            except:
                print(f"Errorr fail with {err_ix}")

    else:
        print('skip errors')
        pass

    print(
        f"Fit done;  chi2 = {Fit.statistic} for {Fit.dof} dof, chi2_red = {Fit.statistic/Fit.dof}, H0 prob = {Fit.nullhyp}")
    #save model and model+data for quick access
    xspec.Xset.save(
        f'xspec/{model_name}/xcm/{data_name}.xcm_model', info='m')
    xspec.Xset.save(f'xspec/{model_name}/xcm/{data_name}.xcm', info='a')

    try:
        plot_spe_nu(Plot, min_bin=min_bin, min_sig=min_sig,
                    figname=f'xspec/{model_name}/{data_name}.png')
    except:
        plt.close('all')
        plot_spe_nu(Plot,  min_bin=min_bin, min_sig=min_sig,
                    figname=f'xspec/{model_name}/{data_name}.png')

    #container objects from storage.py by P. Medvedev
    cA = Container(data_name+'_FPMA')
    cA.get_session(group=1)
    cB = Container(data_name+'_FPMB')
    cB.get_session(group=2)

    if eqw_comps is  not None:
        eqw_comps = [m1.componentNames.index(comp)+1 for comp in eqw_comps]
        for comp in eqw_comps:
            if eqw_cl==0:
                AllModels.eqwidth(comp)
                eqw = S1.eqwidth
            else:
                try:
                        AllModels.eqwidth(comp, err=True, number=100, level=eqw_cl)
                        eqw = S1.eqwidth
                except:
                        print('eqw error calculation failed')
                        AllModels.eqwidth(comp)
                        eqw = S1.eqwidth

            eqw_row = pd.Series(
                {'comp': m1.componentNames[comp-1],
                'par': 'eqw',
                'ipar': 0,
                'val': eqw[0],
                'error_l': eqw[1],
                'error_u': eqw[2],
                'er_status': 'FFFFFFFFF',
                'sigma': 0.0,
                'frozen': False,
                'link': ''})
            cA._params = pd.concat([cA.params, eqw_row.to_frame().T])
            #cA._params = cA.params.append(eqw_row, ignore_index=True)
            #cA._params = pd.concat([cA.params, eqw_row], ignore_index = True) 
    
    s = Storage()
    s(cA)
    s(cB)

    s.to_pickle(f'xspec/{model_name}/{data_name}.storage')
    print(s)
    return s


# def scan_containers_ph_ave(model):
#     tmp_list = []
#     for storage in glob(f'xspec/{model}/*storage'):
#         s = Storage().from_pikle(storage)
#         fpma = s[0].params
#         fpmb = s[1].params
#         fit = s[0].fit
#         # splitting a string like this: 90302319002_bbody_FPMA
#         ObsID, model, detA = s._srcID[0].split('_')
#         _, _,  detB = s._srcID[1].split('_')

#         for df in [fpma, fpmb]:
#             # df['ObsID']='obs'+ObsID
#             df['ObsID'] = ObsID
#             df['model'] = model
#             #df['statistic'] = fit.statistic
#             #df['dof'] = fit.dof

#         fpma['det'] = detA
#         fpmb['det'] = detB

#         tmp_list.append(fpma)
#         tmp_list.append(fpmb)
#         # tmp_list.append(fit)

#     ph_ave_results = pd.concat(tmp_list)
#     ph_ave_results = ph_ave_results.drop(
#         ['ipar', 'sigma', 'link'], axis=1)  # delete unneeded columns

#     ph_ave_results.index = range(len(ph_ave_results))

#     ph_ave_results = ph_ave_results.drop(ph_ave_results[(
#         ph_ave_results.det == 'FPMB') & (ph_ave_results.par != 'factor')].index)
#     ph_ave_results = ph_ave_results.drop(ph_ave_results[(
#         ph_ave_results.det == 'FPMA') & (ph_ave_results.par == 'factor')].index)
#     ph_ave_results = ph_ave_results.drop(
#         ['det'], axis=1)  # delete unneeded columns

#     ph_ave_results_reind = ph_ave_results.set_index(
#         ['ObsID','model', 'comp',  'par'])
#     chi2_str = f"chi2 {fit.statistic.iloc[0]:.2f}/{fit.dof.iloc[0]:.0f}"
#     ph_ave_results_reind.loc[ObsID, model, 'stat', 'chi2'] = chi2_str
#     ph_ave_results_reind.loc[ObsID, model, 'flux', 'flux'] = 'chi2 ---'
#     return ph_ave_results_reind


# def query_par(fit_res, ObsID, model, comp, par, shift=0):
    
#     df = fit_res.loc[(ObsID, shift, model,    comp, par)].sort_values('phase')
#     title = ('.').join([ObsID, model, comp, par])
#     return df, title



def scan_containers_ph_res(model_name:  str) -> pd.DataFrame:
    """
    scan_containers_ph_res scans xspec containers for the results of the  phase-resolved  spectrum fitting.

    Args:
        model_name (str): model name

    Returns:
        pd.DataFrame: DataFrame with  fitting result. Pivot = [component name, parameter name, phase]
    """

    tmp_list = []
    #scan for storage pickles
    for storage in glob(f'xspec/{model_name}/*storage'):
        s = Storage().from_pikle(storage)
        fpma = s[0].params
        fpmb = s[1].params
        # splitting a string like this: 90302319002_bin2_bbody_FPMA
        ObsID, binnum, model, detA = s._srcID[0].split('_')
        _, _, _, detB = s._srcID[1].split('_')

        for df in [fpma, fpmb]:
            # df['ObsID']='obs'+ObsID
            df['ObsID'] = ObsID
            df['binnum'] = int(binnum[3:])
            df['model'] = model

        fpma['det'] = detA
        fpmb['det'] = detB

        tmp_list.append(fpma)
        tmp_list.append(fpmb)

    #make a  data frame
    ph_res_results = pd.concat(tmp_list)
    
    #delete linked pars
    ph_res_results=ph_res_results[ph_res_results.link=='']
    # delete unneeded columns
    ph_res_results = ph_res_results.drop(
        ['ipar', 'sigma', 'link'], axis=1)  
    ph_res_results = ph_res_results.reset_index()
    
    #set phase 
    nph = ph_res_results['binnum'].max()
    ph_res_results['binnum'] = ph_res_results['binnum'] / nph - 0.5/(nph)
    ph_res_results = ph_res_results.rename(columns = {'binnum':'phase'})

    #drop fpmb info, including constant term
    ph_res_results = ph_res_results.drop(ph_res_results[(
        (ph_res_results.det == 'FPMB') & (ph_res_results.par != 'factor'))].index, axis = 0)

    #drop fpma constant info
    ph_res_results = ph_res_results.drop(ph_res_results[(
        ph_res_results.det == 'FPMA') & (ph_res_results.par == 'factor')].index)
    
    #drop columns
    ph_res_results = ph_res_results.drop(['det'], axis = 1)
    ph_res_results = ph_res_results.drop(['index'], axis = 1)

    #sort by phase
    ph_res_results = ph_res_results.sort_values(['phase'], ascending=True)

    #set pivot
    ph_res_results = ph_res_results.pivot_table(index = ['comp','par', 'phase'])
    return ph_res_results



def ph_res_param(
    df: pd.DataFrame,
    comp:str, par:str,
    funct=lambda x: x,
    plot:  bool =True, ax:Optional[plt.Axes]=None,
    **plt_kwargs) -> Tuple: 
    """
    ph_res_param plots a phase-resolved spectral parameters for a given data frame and components. Produces plots with phase from 0 to 2pi.

    Args:
        df (pd.DataFrame): dataframe with fitting results. Pivot as in scan_containers_ph_res
        comp (str): component name
        par (str): parameter name
        funct (function, optional): function to act on the values, e.g. lambda x: x/100. Defaults to lambdax:x.
        plot (bool, optional): whether to  plot  the result. Defaults to True.
        ax (Optional[plt.Axes], optional): axis to plot result on. Defaults to None.

    Raises:
        Exception: _description_

    Returns:
        Tuple: phase, value and error
    """

    ser = df.loc[pd.IndexSlice[comp, par, :]]
    def get_parr_array(ser, funct):

        mean = funct(ser.val)
        lo = funct(ser.error_l)
        hi = funct(ser.error_u)

        hi = (hi-mean)
        lo = (mean-lo)
        err = np.vstack((lo, hi))

        return mean, err

    phase = ser.index.values#ser.phase
    try:
        mean, err = get_parr_array(ser, funct)
    except:
        raise Exception(f'Error with get par {ser}')

    phase = np.concatenate((phase, phase+1))
    #dphase = np.diff(phase)[0]
    mean = np.concatenate((mean, mean))
    err = np.hstack((err, err))

    if plot:
        if ax == None:
            fig, ax = plt.subplots(figsize=(12, 4))
        
        ax.plot(phase, mean, alpha = 0, label ='_pass', color = 'white')
        ax.set_ylim()
        ax.errorbar(phase, mean, yerr=err, drawstyle='steps-mid', **plt_kwargs)

        ax.set_xlabel('phase')
        ax.set_ylabel(comp+':'+par)

    return phase, mean, err






def plot_ph_res_storage(df_ph_res: pd.DataFrame,  nu_obs: NustarObservation, prodpath_ph_res: str, nrows: int=3) -> plt.Figure:
    """
    plot_ph_res_storage plots the whole contained with phase-resolved parameters. Does not  plot frozen parameters.

    Args:
        df_ph_res (pd.DataFrame): DataFrame with phase-resolved results. Pivot as in scan_containers_ph_res
        nu_obs (NustarObservation): NustarObservation object needed for plotting pulse profiles.
        prodpath_ph_res (str): path to phase-resolved products (for plotting pulse profiles of all bins)
        nrows (int, optional): number or rows in a plot. Defaults to 3.
    
    returns:
        Figure: figure with phase-resolved plots
    """


    #reset index to get rid of phase column
    df_ph_res_reset = df_ph_res.reset_index()
    df_ph_res_reset  = df_ph_res_reset[df_ph_res_reset.frozen  == False]

    #https://stackoverflow.com/questions/29975835/how-to-create-pandas-groupby-plot-with-subplots
    #group by components
    grouped = df_ph_res_reset.groupby(['comp', 'par'])
    rowlength = int(grouped.ngroups/nrows)+1
    
    fig, axs = plt.subplots(figsize=(14,8), 
                            nrows=nrows, ncols=rowlength,    
                            gridspec_kw=dict(hspace=0.0),  squeeze=True) 
    axs = axs.flatten()
    targets = zip(grouped.groups.keys(), axs)
    for (key, ax) in targets:
        df_tmp = grouped.get_group(key).set_index(['comp', 'par', 'phase'])
        ph_res_param(df=df_tmp, comp=key[0], par=key[1], ax=ax)
        ax.set_ylabel(key)


    axs[0].clear()
    efolds = glob('*.efold')
    _, colors = nu_obs.plot_efolds_of_bins(prodpath=prodpath_ph_res,        efolds_files=efolds, ax_efold=axs[0], fig=fig,
                                               save=False, legend=False, phase_zero_efold_file=nu_obs.products_path+'/phase_resolved/'+'phase_resolved_bin1AB_sr.lc_bary_nphase_128.efold')

    plt.show()

    return fig

