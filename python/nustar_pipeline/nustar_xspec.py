# %% imports and definitions
from matplotlib import gridspec
from python_for_nustar.nu_core import set_mpl
set_mpl()
import os
from python_for_nustar.pyxspec_lib.definitions import plt, make_figure, create_dir, showmodel, xspec, Model, Fit, Plot, AllData, AllModels, np, pd, glob
from python_for_nustar.pyxspec_lib.storage import Container, Storage
Fit.query = "yes"
AllData.clear()
Plot.device = '/null'
SPE_TIMEOUT = 30  # after this much seconds of fitting, the function would produce an error due to timeout


def plot_spe_nu(Plot, ax_spe=None, ax_res=None, fig=None, label='',
                min_sig=50, min_bin=50,
                    color1 = 'C0',
                    color2 = 'C1',
                plotmode='eeufs', figname=None, residuals='rat'):
    ms = 5
    alpha = 0.4
    if ax_spe is None:
        fig, ax_spe, ax_res = make_figure(residuals=residuals)
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

        if residuals == 'rat':
            Plot('rat')
            Plot.setRebin(min_sig, min_bin)
            Plot.xAxis = 'keV'

            rat = Plot.y(plotGroup=group)
            rat_err = Plot.yErr(plotGroup=group)

            ax_res.plot(en, rat, 'o', label=label, color=color, alpha=alpha, ms = ms)
            ax_res.errorbar(en, rat, rat_err, en_err,
                            fmt='none', ecolor=color, alpha=alpha)
        if residuals == 'del':
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

    return en, en_err,  data, data_err, model, rat, rat_err

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
    print(f' loaded {spe1} and {spe2}')
    return load_nu_spe(spe1=spe1, spe2=spe2, en_lo=en_lo, en_hi=en_hi)


# @timeout_decorator.timeout(SPE_TIMEOUT)
def fit_spectra(model: xspec.model.Model,
                prefix: str,  model_name: str,
                dataset: str = 'spe_and_lc',
                error_conf: str = '2.71',
                en_lo='4.', en_hi='79.',
                min_sig=50, min_bin=50,
                rewrite: bool = True,
                calc_errors=True,
                perturb_fit_sigma=0.0,
                ignore_comp_errors=[],
                eqw_comps=['gaussian'],
                eqw_cl=90):
    Fit.query = "yes"
    assert '_' not in model_name,  'model_name cannot contain "_"'
    data_name = prefix + '_'+model_name
    create_dir(f'xspec/{model_name}')
    create_dir(f'xspec/{model_name}/xcm')
    if rewrite:
        os.system(f'rm -r xspec/{model_name}/*{data_name}*')
        os.system(f'rm -r xspec/{model_name}/xcm/*{data_name}*')
        print(
            f"deleted xspec/{model_name}/*{data_name}* and xspec/{model_name}/xcm/*{data_name}* files")

    try:
        print(f'Try loading storage xspec/{model_name}/{data_name}.storage')
        s = Storage().from_pikle(f'xspec/{model_name}/{data_name}.storage')
        print('loaded')
        return s
    except:
        print('no storage found. fitting...')
        pass

    spe1 = f'{dataset}A_sr.pi'
    spe2 = f'{dataset}B_sr.pi'
    S1, S2 = load_nu_spe(spe1=spe1, spe2=spe2, en_lo=en_lo, en_hi=en_hi)
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

    xspec.Xset.save(
        f'xspec/{model_name}/xcm/{data_name}.xcm_model_init', info='m')

    Fit.query = 'yes'
    Fit.statMethod = "chi"

    Fit.query = "yes"

    if perturb_fit_sigma != 0:
        print('perturbing ininital guess before fitting')
        for comp_name in m1.componentNames:
            comp = getattr(m1, comp_name)
            for par_name in comp.parameterNames:
                par = getattr(comp, par_name)
                if not par.frozen:
                    val, delta, _, _, _, _ = par.values
                    # print(comp.name, par.name)
                    par.values = np.random.normal(val, delta*perturb_fit_sigma)
        # showmodel(m1)

    Fit.perform()
    print('fitting done')
    if calc_errors:
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
                Fit.error(f"nonew maximum 20. {error_conf} {err_ix}")
            except:
                print(f"Errorr fail with {err_ix}")

    else:
        print('skip errors')
        pass

    print(
        f"Fit done \n chi2 = {Fit.statistic} for {Fit.dof} dof, chi2_red = {Fit.statistic/Fit.dof}, H0 prob = {Fit.nullhyp}")

    xspec.Xset.save(
        f'xspec/{model_name}/xcm/{data_name}.xcm_model', info='m')
    xspec.Xset.save(f'xspec/{model_name}/xcm/{data_name}.xcm', info='a')

    try:
        plot_spe_nu(Plot, residuals='del', min_bin=min_bin, min_sig=min_sig,
                    figname=f'xspec/{model_name}/{data_name}.png')
    except:
        plt.close('all')
        plot_spe_nu(Plot, residuals='del', min_bin=min_bin, min_sig=min_sig,
                    figname=f'xspec/{model_name}/{data_name}.png')

    cA = Container(data_name+'_FPMA')
    cA.get_session(group=1)
    cB = Container(data_name+'_FPMB')
    cB.get_session(group=2)

    eqw_comps = [m1.componentNames.index(comp)+1 for comp in eqw_comps]
    if len(eqw_comps)!=0:
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


def ph_res_param(ser,
                 funct=lambda x: x,
                 plot=1, ax=None, title='', colors=None,
                 delta_limits=0.05, plot_relative=False,
                 **plt_kwargs):

    def get_parr_array(ser, funct):

        mean = funct(ser.val)
        lo = funct(ser.error_l)
        hi = funct(ser.error_u)

        hi = (hi-mean)
        lo = (mean-lo)
        err = np.vstack((lo, hi))

        return mean, err

    phase = ser.phase
    try:
        mean, err = get_parr_array(ser, funct)
    except:
        raise Exception(f'Error with get par {ser}')

    phase = np.concatenate((phase, phase+1))
    dphase = np.diff(phase)[0]
    mean = np.concatenate((mean, mean))
    err = np.hstack((err, err))
    if colors is not None:
        colors  = colors[1:] + colors[:]

    if plot:
        if ax == None:
            fig, ax = plt.subplots(figsize=(12, 4))
            ax.set_title(title)
        
        ax.plot(phase, mean, alpha = 0, label ='_pass', color = 'white')
        ax.set_ylim()
        ax.errorbar(phase, mean, yerr=err, drawstyle='steps-mid', **plt_kwargs)
        if colors is not None:
            for phase_bin, color_bin in zip(phase, colors):
                ax.axvspan(phase_bin-dphase/2, phase_bin+dphase/2, color = color_bin, alpha = 0.2)
                #ax.axvspan(phase_bin-dphase/2+1, phase_bin+dphase/2+1, color = color_bin, alpha = 0.2)
        ax.set_xlabel('phase')
        ax.set_ylabel('')
        valmin, valmax = np.min(mean), np.max(mean)
        #valmean = np.mean(mean)
        #delta = valmax - valmin
        # ax.set_ylim(valmean-delta_limits*delta, valmean+delta_limits*delta)
        #ax.set_ylim((1-delta_limits)*valmin, (1+delta_limits)*valmax)
        if plot_relative:
            axt = ax.twinx()
            factor = np.ma.average(mean,
                                   weights=np.max(err, axis=0))
            axt.errorbar(
                phase, mean/factor, err/factor, alpha=0.3, fmt='-', color='k', ecolor='k')
            axt.grid(False)
            axt.set_ylabel('rel. change')
        ax.legend()
    return phase, mean, err


def scan_containers_ph_res(model_name, ph_res_folder = None):
    if ph_res_folder is not None:
        os.chdir(ph_res_folder)
    tmp_list = []
    for storage in glob(f'xspec/{model_name}/*storage'):
        s = Storage().from_pikle(storage)
        fpma = s[0].params
        fpmb = s[1].params
        fit = s[0].fit
        # splitting a string like this: 90302319002_bin2_shift_0_bbody_FPMA
        ObsID, binnum, shift, model, detA = s._srcID[0].split('_')
        _, _, _, _, detB = s._srcID[1].split('_')

        for df in [fpma, fpmb]:
            # df['ObsID']='obs'+ObsID
            df['ObsID'] = ObsID
            df['binnum'] = int(binnum[3:])
            df['shift'] = int(shift[5:])
            df['model'] = model
            df['statistic'] = fit.statistic
            df['dof'] = fit.dof

        fpma['det'] = detA
        fpmb['det'] = detB

        tmp_list.append(fpma)
        tmp_list.append(fpmb)
        # tmp_list.append(fit)

    ph_res_results = pd.concat(tmp_list)
    ph_res_results = ph_res_results.drop(
        ['ipar', 'sigma', 'link'], axis=1)  # delete unneeded columns
    ph_res_results.index = range(len(ph_res_results))
    nph = ph_res_results['binnum'].max()
    ph_res_results['phase'] = ph_res_results['binnum'] / nph - 0.5/(nph)
    ph_res_results = ph_res_results.drop(ph_res_results[(
        ph_res_results.det == 'FPMB') & (ph_res_results.par != 'factor')].index)
    ph_res_results = ph_res_results.drop(ph_res_results[(
        ph_res_results.det == 'FPMA') & (ph_res_results.par == 'factor')].index)
    ph_res_results = ph_res_results.drop(
        ['det', 'binnum'], axis=1)  # delete unneeded columns

    ph_res_results_reind = ph_res_results.set_index(
        ['ObsID', 'shift', 'model', 'comp',  'par'])
    ph_res_results_reind = ph_res_results_reind.sort_index()

    return ph_res_results_reind


def scan_containers_ph_ave(model = '*'):
    tmp_list = []
    for storage in glob(f'xspec/{model}/*storage'):
        s = Storage().from_pikle(storage)
        fpma = s[0].params
        fpmb = s[1].params
        fit = s[0].fit
        # splitting a string like this: 90302319002_bbody_FPMA
        ObsID, model, detA = s._srcID[0].split('_')
        _, _,  detB = s._srcID[1].split('_')

        for df in [fpma, fpmb]:
            # df['ObsID']='obs'+ObsID
            df['ObsID'] = ObsID
            df['model'] = model
            #df['statistic'] = fit.statistic
            #df['dof'] = fit.dof

        fpma['det'] = detA
        fpmb['det'] = detB

        tmp_list.append(fpma)
        tmp_list.append(fpmb)
        # tmp_list.append(fit)

    ph_ave_results = pd.concat(tmp_list)
    ph_ave_results = ph_ave_results.drop(
        ['ipar', 'sigma', 'link'], axis=1)  # delete unneeded columns

    ph_ave_results.index = range(len(ph_ave_results))

    ph_ave_results = ph_ave_results.drop(ph_ave_results[(
        ph_ave_results.det == 'FPMB') & (ph_ave_results.par != 'factor')].index)
    ph_ave_results = ph_ave_results.drop(ph_ave_results[(
        ph_ave_results.det == 'FPMA') & (ph_ave_results.par == 'factor')].index)
    ph_ave_results = ph_ave_results.drop(
        ['det'], axis=1)  # delete unneeded columns

    ph_ave_results_reind = ph_ave_results.set_index(
        ['ObsID','model', 'comp',  'par'])
    chi2_str = f"chi2 {fit.statistic.iloc[0]:.2f}/{fit.dof.iloc[0]:.0f}"
    ph_ave_results_reind.loc[ObsID, model, 'stat', 'chi2'] = chi2_str
    ph_ave_results_reind.loc[ObsID, model, 'flux', 'flux'] = 'chi2 ---'
    return ph_ave_results_reind


def query_par(fit_res, ObsID, model, comp, par, shift=0):
    
    df = fit_res.loc[(ObsID, shift, model,    comp, par)].sort_values('phase')
    title = ('.').join([ObsID, model, comp, par])
    return df, title


def query_par_ph_ave(fit_res, model, comp, par):
    df = fit_res.loc[(model,    comp, par)]  # .sort_values('ObsID')
    return df


def plot_ph_res_storage(ph_res_results,  nu_obs, prodpath_ph_res, cols=3):
    variable_pars = ph_res_results[(ph_res_results.frozen == False) & (
        ph_res_results.index.get_level_values('par') != 'factor')]
    frozen_pars = ph_res_results[ph_res_results.frozen == True]
    frozen_pars = [f"{par[3]}:{par[4]}" for par in frozen_pars.index.unique()]
    params = variable_pars.index.unique()
    npars = len(params)+1

    nplots = npars+cols
    rows = int(np.ceil(nplots / cols))

    gs = gridspec.GridSpec(rows, cols)
    fig = plt.figure(figsize=(14, 8))

    ax0 = fig.add_subplot(gs[0])
    # plt.text(x=0.1, y=0.94,
    #         s=f"frozen params: {frozen_pars}", fontsize=6, transform=fig.transFigure)
    ax0.set_title(
        f"frozen params: {frozen_pars[:len(frozen_pars)//2]} \n {frozen_pars[len(frozen_pars)//2:]}", fontsize=6,)
    for ii in range(cols):
        ax = fig.add_subplot(gs[ii], sharex=ax0)
        efolds = glob('*.efold')
        _, colors = nu_obs.check_efold_of_bins(prodpath=prodpath_ph_res, efolds_files=efolds, fiducial=None, ax_efold=ax, fig=fig,
                                               save=False, legend=False, phase_zero_efold='phase_resolved_bin1AB_sr.lc_bary_orb_corr_nphase_128.efold')

    for ii, par in enumerate(params, cols):
        ax = fig.add_subplot(gs[ii], sharex=ax0)
        ObsID, prod_shift, model_name, comp, parname = par
        df, title = query_par(fit_res=ph_res_results, ObsID=ObsID,
                              model=model_name,    comp=comp,     par=parname, shift=prod_shift)
        ph_res_param(df, label=f"{comp}:{parname}",  funct=lambda x: x,
                     alpha=0.6, color='k', lw=2,  ax=ax)
    fig.suptitle(f"{ObsID} {model_name}")

    ax = fig.add_subplot(gs[ii+1], sharex=ax0)
    ax.scatter(ph_res_results.phase, ph_res_results.statistic/ph_res_results.dof)
    ax.scatter(ph_res_results.phase+1, ph_res_results.statistic/ph_res_results.dof)
    ax.legend(['chi^2/dof'])
    fig.tight_layout()
    plt.subplots_adjust(hspace=0)
    return fig
