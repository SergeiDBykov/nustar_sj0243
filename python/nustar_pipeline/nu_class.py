
from python_for_nustar.nu_core import (
    np,
    plt,
    os,
    create_dir,
    glob,
    fits,
    pd,
    run_command,
    fit_efsearch_data,
    start_stop,
    fits,
    EfromPI,
    ratio_error,
)

fpm = "FPM"
modules = ["A", "B"]


class NustarObservation:
    def __init__(self, ObsID: str, nu_path: str):
        """
        __init__ initialize observation analysis: create directories for products and assign directory for data

        Args:
            ObsID (str): ObsID of NuSTAR observation
            nu_path (str): path to nustar data of a particular source (e.g. /Users/sdbykov/work/xray_pulsars/groj2058_nu/)
            obs_params (dict): dictionary of orbital parameters used for orbital correction
        """

        print("###")
        print(f"Observation {ObsID} loaded successfully")
        self.ObsID = ObsID
        self.data_path = nu_path + "data/" + ObsID + "/"
        os.chdir(nu_path + "results/")

        create_dir("out" + self.ObsID)
        os.chdir("out" + ObsID)
        out_path = os.getcwd()
        self.out_path = out_path
        create_dir("products")
        os.chdir("products")
        self.products_path = os.getcwd()


###### NUPIPELINE ########

    def nupipeline(self, ObsList_bright: list):
        """
        nupipeline creates a command fo nupipeline of observation. if ObsID is in list ObsList_brights, sets status expression as recommendeed for bright sources

        Args:
            ObsList_bright (List): list of observations when source was bright (>100 cts/s, see NuSTAR tutorial). You can assess it by looking at light curve pre-peoduced by NusTAR team in event_cl folder of observation: from data folder run $ open */*/*A*lc* an check light curves by eye

        """

        os.chdir(self.out_path)
        if len(glob("*A01_cl.evt")) != 0:
            raise Exception("nupipeline has already been launched!")

        nupipeline = f"""nupipeline indir={self.data_path} steminputs=nu{self.ObsID} outdir={self.out_path} obsmode=SCIENCE createexpomap=yes"""
        if self.ObsID in ObsList_bright:
            nupipeline += ' statusexpr="STATUS==b0000xxx00xxxx000"'
        run_command(cmd=nupipeline, cmd_name="nupipeline", rewrite=True)

    def make_regions(self):
        """
        make_regions creates a ds9 script which plots images from modules and you can easily set and check region positions for source and background estimation
        """

        ds9_set = f"ds9 -tile nu{self.ObsID}A01_cl.evt -scale log  -tile nu{self.ObsID}B01_cl.evt  -scale log -lock scale"
        run_command(cmd=ds9_set, cmd_name="ds9_set")

        ds9_check = f"ds9 -tile nu{self.ObsID}A01_cl.evt -regions srcA.reg -tile nu{self.ObsID}A01_cl.evt   -regions bkgA.reg -tile nu{self.ObsID}B01_cl.evt -region srcB.reg  -tile nu{self.ObsID}B01_cl.evt -regions bkgB.reg  -scale log  -lock scale -lock frame wcs"

        run_command(cmd=ds9_check, cmd_name="ds9_check")

    ## NUPRODUCTS ########

    def nuproducts(
        self,
        mode: str,
        outdir: str,
        stemout="DEFAULT",
        bkgextract: str = "yes",
        phafile: str = "DEFAULT",
        bkgphafile: str = "DEFAULT",
        runmkrmf: str = "yes",
        runmkarf: str = "yes",
        lcfile: str = "DEFAULT",
        bkglcfile: str = "DEFAULT",
        imagefile: str = "DEFAULT",
        pilow: str = "60",
        pihigh: str = "1935",
        lcenergy = '10',
        binsize: float = 0.01,
        usrgtifile: str = "NONE",
        usrgtibarycorr: str = "no",
        rewrite=True,
    ):
        """
        nuproducts  executes a nuproduct command for nustar scientific products of level 3
        https://heasarc.gsfc.nasa.gov/lheasoft/ftools/caldb/help/nuproducts.html


        Args:
            mode (str): module of the FPM, A or B
            outdir (str): output directory name and stem for data files. E.g. 'spe_and_lc'
            stemout (str): nuproducts argument. Defaults to 'DEFAULT'.
            bkgextract (str, optional): nuproducts argument. Defaults to 'yes'.
            phafile (str, optional): nuproducts argument. Defaults to 'DEFAULT'.
            bkgphafile (str, optional): nuproducts argument. Defaults to 'DEFAULT'.
            runmkrmf (str, optional): nuproducts argument. Defaults to 'yes'.
            runmkarf (str, optional): nuproducts argument. Defaults to 'yes'.
            lcfile (str, optional): nuproducts argument. Defaults to 'DEFAULT'.
            bkglcfile (str, optional): nuproducts argument. Defaults to 'DEFAULT'.
            imagefile (str, optional): nuproducts argument. Defaults to 'DEFAULT'.
            pilow (str, optional): nuproducts argument. Defaults to '60'.
            pihigh (str, optional): nuproducts argument. Defaults to '1935'.
            binsize ( , optional): nuproducts argument. Defaults to 0.01.
            usrgtifile (str, optional): nuproducts argument. Defaults to 'NONE'.
            usrgtibarycorr (str, optional): nuproducts argument. Defaults to 'no'.

        """

        if mode not in modules:
            raise Exception("choose module name from A,B")

        os.chdir(self.products_path)
        create_dir(outdir)

        ObsID = self.ObsID
        outdi = outdir

        if stemout != "DEFAULT":
            stemout = stemout + mode

        nuproducts = f"""
    nuproducts \
    indir={self.out_path} \
    instrument={fpm}{mode} \
    steminputs=nu{ObsID} \
    stemout={stemout} \
    outdir={outdi} \
    srcregionfile={self.out_path}/src{mode}.reg \
    bkgextract={bkgextract} \
    bkgregionfile={self.out_path}/bkg{mode}.reg \
    binsize={binsize} \
    lcfile={lcfile} \
    phafile={phafile} \
    bkglcfile={bkglcfile} \
    bkgphafile={bkgphafile} \
    imagefile={imagefile} \
    usrgtifile={usrgtifile} \
    pilow={pilow} pihigh={pihigh} \
    lcenergy={lcenergy} \
    usrgtibarycorr={usrgtibarycorr}\
    runmkarf={runmkarf}\
    runmkrmf={runmkrmf}"""

        run_command(cmd=nuproducts, cmd_name=outdir + mode, rewrite=rewrite)

    ##### SPECTRA ######

    def make_spe(self, **kwargs):
        """
        make_spe: nuproducts without creation of lightcurve or images, arguments apply
        """
        self.nuproducts(**kwargs, lcfile="NONE",
                        bkglcfile="NONE", imagefile="NONE")

    def grppha(self, infiles: list, prodpath: str, group_min=30):
        """
        grppha createes a grppha command to group photons into bins

        Args:
            infiles (List): list of spectra to apply barycorr to
            prodpath (str): path where all spectra lie
            group_min (int, optional): grppha argument. Defaults to 30.
        """

        os.chdir(self.products_path + "/" + prodpath)
        for infile in infiles:
            outfile = infile[0:-4] + ".pi"
            grppha = f"""grppha infile="{infile}" outfile="{outfile}"  comm="group min {group_min} & exit" clobber=yes"""
            run_command(cmd=grppha, cmd_name="grppha", rewrite=False)

    ##### LIGHT CURVES ######

    def make_lc(self, **kwargs):
        """
        make_lc: nuproducts without spectra and image,  arguments apply
        """
        self.nuproducts(
            **kwargs,
            phafile="NONE",
            bkgphafile="NONE",
            runmkrmf="NONE",
            runmkarf="NONE",
            imagefile="NONE",
        )

    def barycorr(
        self,
        infiles: list,
        prodpath: str,
        barytime: str = "no",
        cmd_name: str = "barycorr",
        rewrite=0,
    ):
        """
        barycorr applies barycorr to input file. NOTE: Ra and Dec of the source are readed from fits file. If you need, you may add keywords ra, dec to barocorr command belowself.

        Args:
            infile (List): input file names
            prodpath (str): path to products whene file lies
            barytime (str, optional): Argument of barycor ftools task. Defaults to 'no'.
            cmd_name (str, optional): command name. Defaults to 'lcmath'.
            rewrite (bool, optional): whether to rewrit. Defaults to True.

        """
        os.chdir(self.products_path + "/" + prodpath)
        for infile in infiles:
            outfile = infile + "_bary"
            barycorr = f"""      barycorr \
        infile={infile} \
        outfile={outfile}\
        orbitfiles={self.data_path}auxil/nu{self.ObsID}_orb.fits.gz\
        barytime={barytime} """
            run_command(barycorr, cmd_name=cmd_name, rewrite=rewrite)

    def lcmath(
        self,
        infiles: list,
        outfile: str,
        prodpath: str,
        cmd_name: str = "lcmath",
        rewrite: bool = True,
    ):
        """
        lcmath creates a script for running lcmath ftools routine to AVERAGEE two input lightcurves

        Args:
            infiles (List): list of two lightcurvs to average
            outfile (str): output filenamee
            prodpath (str): path to light curves
            cmd_name (str, optional): command name. Defaults to 'lcmath'.
            rewrite (bool, optional): whether to rewrit. Defaults to True.
        """

        os.chdir(self.products_path + "/" + prodpath)

        infile1, infile2 = infiles
        lcmath_cmd = f"""
        lcmath \
        infile={infile1} \
        bgfile = {infile2} \
        outfile = {outfile} \
        multi = 0.5 \
        multb = 0.5 \
        addsubr = yes
        """
        run_command(lcmath_cmd, cmd_name=cmd_name, rewrite=rewrite)

    def orb_correction_lc(self, filename: str, prodpath: str):
        """
        orb_correction_lc creates a new file of a lightcurve with times corrected for orbital motion

        I use the correction of LC file because the correction of event files
        would take a lot of time

        Args:
            filename (str): filename of a lighcurve to be corrected
            prodpath (str): path to lightcurve product folder
        """

        os.chdir(self.products_path + "/" + prodpath)

        #doppler_correction.correct_times(  # type: ignore
        #    fitsfile=filename, orb_params=None, time_orig_col="time"
        #)
        filename_path, filename_only = filename.rsplit('/', 1)
        new_filename = filename_path+'/'+filename_only+'_orb_corr'
        os.system(f'cp {filename} {new_filename}')

    ##### PERIOD SEARCH AND EPOCH FOLDING ######

    def make_efsearch(
        self,
        filename: str,
        prodpath: str,
        p0: str,
        nphase: str = "8",
        dres: str = "0.001",
        nper: str = "128",
        cmd_name: str = "efsearch",
        rewrite=True,
    ):
        """
        make_efsearch creates a script for running eefsearch for a given light curve

        Args:
            filename (str): filename of a light curve (barycentred)
            prodpath (str): [description]
            p0 (str): efsearch argument
            nphase (str, optional): efsearch argument. Defaults to '8'.
            dres (str, optional): efsearch argument. Defaults to '0.001'.
            nper (str, optional): efsearch argument. Defaults to '128'.
            prodpath (str): path to light curves
            cmd_name (str, optional): command name. Defaults to 'lcmath'.
            rewrite (bool, optional): whether to rewrit. Defaults to True.

        """
        os.chdir(self.products_path + "/" + prodpath)

        efsearch = f"""
        efsearch \
            cfile1="{filename}" \
                dper={p0} \
                nphase={nphase} \
                dres={dres} \
                nper={nper} \
                outfile="{filename}.efs" \
                window="-" \
                sepoch=INDEF \
                nbint=INDEF \
                plot=yes \
                plotdev="/xw"
                """
        run_command(efsearch, cmd_name=cmd_name, rewrite=rewrite)

    def fit_efsearch(self, filename: str, prodpath: str):
        """
        fit_efsearch fits efsearch file

        Args:
            filename (str): filename
            prodpath (str): path to file
        """
        os.chdir(self.products_path + "/" + prodpath)
        per, err = fit_efsearch_data(filename + ".efs")
        return per, err

    def make_efold(
        self,
        filename: str,
        prodpath: str,
        period: str,
        nphase: str = "32",
        cmd_name: str = "efold",
        rewrite=True,
    ):
        """
        make_efold makes efold script to fold loghtcurves with given period

        Args:
            filename (str): name of a lc
            prodpath (str): path to an lc
            period (str): period to fold with
            nphase (str, optional): number of phases. Defaults to '32'.
            cmd_name (str, optional): name of the command. Defaults to 'efold'.
            rewrite (bool, optional): whether to rewrite command file. Defaults to True.
        """

        os.chdir(self.products_path + "/" + prodpath)

        efold = f"""
            efold \
            nser = 1 \
            norm = 0 \
            cfile1 = {filename} \
            dper = {period} \
            nphase = {nphase} \
            nbint = INDEF \
            nintfm = INDEF \
            outfileroot = 'default' \
            window = '-' \
            sepoch = 0 \
            plot = no  \
            outfile = {filename}_nphase_{nphase}.efold"""

        run_command(efold, cmd_name=cmd_name, rewrite=rewrite)

    ###### PHASE RESOLVED SPECTRA AND GTIs #####

    def make_gti_from_lc(
        self,
        prodpath: str,
        mode: str,
        period: float,
        phase_bins: int = 10,
        outfolder: str = "gtis/",
        half_shift: int = 0,
    ):
        """
        make_gti_from_lc
        This task uses BARYCENTRED AND ORBIRALLY CORRECTED LIGHTCURVES to produce GTI with given period and timebins.
        Al long as original lc was barycorrected to different file,
        and this barycentred file was used to produce Orbital correction(see above),
        we juxtapose times from corected LC(TBD time) and original
        LC(times starts with zero). Timezero for the lightcurve is a time of
        the first event in an event file(hence we use time of original lc + timezero)

        Args:
            prodpath (str): path to light curve file from which GTIs are built
            mode (str): mode of FPM telescope. May be A, B or AB
            period (float): period for gti
            phase_bins (int, optional): number of phases. Defaults to 10.
            outfolder (str, optional): output folder for GTIs. Defaults to 'gtis/'.
            half_shift (int, optional): number of half bin size shifts to apply to zero time, If 0, makes no shift and starts with the first bin time as zero-time. Defaults to 0. Example:  If half_shift = 1, the bin that prebiously started at phi = 0.5, would start at phi = 0.5 - (phase_bin_width)/2

        """

        os.chdir(self.products_path + "/" + prodpath)
        gtis_created = glob(f"{outfolder}*gti*.fits")
        if gtis_created:
            print(gtis_created)
            if (
                input(
                    "GTIs for this mode have already been created. Delete folder to rebuild? [y/n]"
                )
                == "y"
            ):
                os.system(f"rm -r {outfolder}")
            else:
                raise Exception("Abort GTI creation")

        create_dir(outfolder)
        filename_orig = glob(f"*{mode}_sr.lc")[0]
        filename_corr = glob(f"*{mode}_sr.lc_bary_orb_corr")[0]

        ff_orig = fits.open(filename_orig)
        time_orig = ff_orig[1].data["time"]  # type: ignore
        time_orig_zero = ff_orig[1].header["timezero"]  # type: ignore
        # no need to add mjdref, because gti is in seconds since ref (see *gti* files in the root directory)
        time_orig += time_orig_zero

        time_corr = fits.open(filename_corr)[1].data["time"]  # type: ignore
        if min(np.diff(time_corr)) < 0:
            raise Exception("corrected times are not monotonous!")

        phase_duration = period / phase_bins
        t0 = time_corr[0]
        phases = (
            np.floor(
                (time_corr - t0 + half_shift * phase_duration / 2)
                / (period / phase_bins)
                % phase_bins
            )
            + 1
        )

        mindt = np.min(np.diff(time_orig))
        for ii in range(1, phase_bins + 1):
            ind = start_stop(phases, trigger_val=ii)
            gtis = time_orig[ind]

            tstart = gtis[:, 0] - mindt / 100
            tstop = gtis[:, 1] + mindt / 100

            hdu = fits.PrimaryHDU()
            c1 = fits.Column(name="START", format="1D", unit="s", array=tstart)
            c2 = fits.Column(name="STOP", format="1D", unit="s", array=tstop)
            cols = fits.ColDefs([c1, c2])
            datahdu = fits.BinTableHDU.from_columns(cols)
            datahdu.name = "GTI"
            datahdu.header["START"] = len(tstart)
            datahdu.header["STOP"] = len(tstop)
            hdulist = fits.HDUList([hdu, datahdu])
            if mode == "AB":
                hdulist.writeto(outfolder + f"gti_{ii}A.fits")
                hdulist.writeto(outfolder + f"gti_{ii}B.fits")
            else:
                hdulist.writeto(outfolder + f"gti_{ii}{mode}.fits")
            hdulist.close()
            print(
                f"GTIs have been created for module(s) {mode} with period {period} and {phase_bins} bins, phase bin number {ii}"
            )

        print(
            f"GTIs have been created for module(s) {mode} with period {period} and {phase_bins} bins"
        )

        # Plotting phase assigned-lightcurve. Need to check this plot before extracting products
        fig, [ax1, ax2] = plt.subplots(2, sharex=True)  # type: ignore

        ind_max = int(3500)
        ax1.plot(time_orig[:ind_max], phases[:ind_max])
        ax1.set_xlabel("original time column")
        ax1.set_ylabel("phase of a light curve bin")

        rate = ff_orig[1].data["rate"]  # type: ignore
        import matplotlib.cm as cm

        cmap = cm.get_cmap("Set3", phase_bins + 1)
        sc = ax2.scatter(
            time_orig[:ind_max], rate[:ind_max], c=phases[:ind_max], cmap=cmap, s=100
        )

        # fig.colorbar(sc, ax=ax2)
        ax2.set_xlabel("original time column")
        ax2.set_ylabel("count rate with color as  phase coding")

    def phase_resolved_spectra(
        self,
        mode: str,
        gtipath="spe_and_lc/gtis",
        folder="phase_resolved",
    ):
        """
        phase_resolved_spectra creates a script for A and B modules, in which nuproducts is run several times (for each phase bin) and extracts spectra AND lightcures for input GTI. Lightcurves are mandatory extracted so that one may check the folded light curves of each phase, and be sure that the spectra are extracted from correct time bins

        Args:
            mode (str): mode of FPM, A or B
            gtipath (str, optional): path to GTI files. Defaults to 'spe_and_lc/gtis'.
            folder (str, optional): folder when nuproduct will unpack spectra and lightcurves. Defaults to 'phase_resolved'.

        CHECK TERMINAL FOR ERRORS IN NUPRODUCTS (search 'nuproducts  error' in terminal window and relaunch problematic phase bins, having deleted beforehand all products created before error). I do not know why errors sometimes arrise, probably due to the conflicts in naming of temporal files and folders while I launch phase resolved for A and B modules in parallel terminal windows. When I launch phase-resolved products for FPMA and only after that FPMB, the errors do not appear (e.g. $ ./phase_resolvedA.sh; ./phase_resolvedB.sh ).

        The time it takes for one phase bin depends on the observation duration and on the lightcurve binning. It may take a few minutes per phase bin (lc+spe).
        """
        os.chdir(self.products_path)
        create_dir(folder)
        gtipath = os.path.abspath(gtipath)
        gtis = glob(f"{gtipath}/*{mode}.fits")
        if gtis == []:
            raise Exception("no GTI found for input parameters")

        phase_bins = len(gtis)
        gtis = [
            f"{gtipath}/gti_{ii}{mode}.fits" for ii in range(1, phase_bins + 1)
        ]  # this is because glob yields strange order of files, even sorted(glob(..)) does not help. It creates incorrect numeration of spectral bins.

        if (
            input(
                f"GTIs from {gtipath} for module {mode}: \n {[os.path.basename(x) for x in gtis]} \n Is this correct? [y/n]?:  "
            )
            == "y"
        ):
            os.chdir(folder)
        else:
            raise Exception("Stop phase resolved spectroscopy")

        for ph_num, gtifile in enumerate(gtis, 1):
            # usrgtibarycorr = no in the next line is very important!!!
            self.nuproducts(
                outdir=folder,
                usrgtifile=gtifile,
                mode=mode,
                stemout=f"phase_resolved_bin{ph_num}",
                usrgtibarycorr="no",
                rewrite=False,
            )

    def check_efold_of_bins(
        self,
        fiducial: str,
        efolds_files: list,
        phase_zero_efold=None,
        prodpath: str = "phase_resolved",
        ax_efold=None,
        fig=None,
        save=True,
        legend=True,
    ):
        """
        check_efold_of_bins builds period-folded lightcurves for each phase bin and compares with the folded lightcurve of the whole observation if necessary.

        Args:
            fiducial (str): efold fits file of the whole observation
            efolds_files (list): efold fits files of phase bins to plot
            prodpath (str, optional): path to phase-resolved products. Defaults to 'phase_resolved'.

        This is how to create a plot of RATIO between phase-bins efilds and fiducial lightcurve
        f0 = fits.open('../spe_and_lc/spe_and_lcAB_sr.lc_bary_orb_corr_nphase_128.efold')
        fig,ax = plt.subplots()

        for binnum in np.arange(1,11):
            fbin = fits.open(f'phase_resolved_bin{binnum}AB_sr.lc_bary_orb_corr_nphase_128.efold')
            ax.plot(f0[1].data['PHASE'], fbin[1].data['RATE1']/f0[1].data['RATE1'], label = binnum)

        ax.legend()


        """
        os.chdir(self.products_path + "/" + prodpath)

        assert len(glob("*A_bin_*sr.pha")) == len(
            glob("*B_bin_*sr.pha")
        ), "different number of spectra files for mode A and B"
        assert len(glob("*A_bin_*sr.lc")) == len(
            glob("*B_bin_*sr.lc")
        ), "different number of light curve files for mode A and B"

        if ax_efold is None or fig is None:
            fig, ax_efold = plt.subplots(figsize=(16, 5))
        else:
            pass

        if phase_zero_efold is None:
            ph_shift = 0
        else:
            ff_zero = fits.open(phase_zero_efold)
            phase_zero = ff_zero[1].data["PHASE"]
            rate_zero = ff_zero[1].data["RATE1"]
            idx = np.where(~np.isnan(rate_zero))[0][0]
            ph_shift = phase_zero[idx] * (-1)
            # print(f'#### PHASE SHIFT  =  {ph_shift}')
        if fiducial is not None:
            for fitsfile in [fiducial]:
                pulse_profile = fits.open(fitsfile)[1].data  # type: ignore
                phase = pulse_profile["PHASE"] + ph_shift
                rate = pulse_profile["RATE1"]
                error = pulse_profile["ERROR1"]
                color = "k"
                ax_efold.errorbar(
                    phase,
                    rate,
                    error,
                    drawstyle="steps-mid",
                    ls="-.",
                    alpha=0.3,
                    lw=5,
                    label="ref",
                    color=color,
                    ecolor=color,
                )
                ax_efold.errorbar(
                    phase + 1,
                    rate,
                    error,
                    drawstyle="steps-mid",
                    ls="-.",
                    alpha=0.3,
                    lw=5,
                    color=color,
                    ecolor=color,
                )
                ymin, ymax = 0.9 * np.min(rate), np.max(rate) * 1.1
        else:
            pass
        colors = []
        for fitsfile in efolds_files:
            # print(f'working with {fitsfile}')
            binnum = fitsfile.split("_")[2]
            pulse_profile = fits.open(fitsfile)[1].data  # type: ignore
            phase = pulse_profile["PHASE"] + ph_shift
            rate = pulse_profile["RATE1"]
            error = pulse_profile["ERROR1"]

            ax_efold.errorbar(
                phase, rate, error, drawstyle="steps-mid", ls="-", alpha=0.8, lw=3
            )

            # mean_phase = 0.9*np.median(phase[rate > 0])

            # binnum_numeral = binnum[3:-2]
            # ax_efold.text(mean_phase, 0.7*np.min(
            #     rate[rate > 0]), s=binnum_numeral, fontsize=6, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
            # ax_efold.text(mean_phase+1, 0.7*np.min(
            #     rate[rate > 0]), s=binnum_numeral, fontsize=6, bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))

            color = ax_efold.get_lines()[-1].get_color()
            colors.append(color)
            for dph in [-2, -1, 0, 1, 2]:
                ax_efold.fill_between(
                    phase + dph, 0, rate, fc=color, alpha=0.4, ec="k", label=f"{binnum}"
                )

            # if binnum == "bin1AB":
            #     for dph in [-2, -1, 0, 1, 2]:
            #         ax_efold.fill_between(
            #             phase + dph,
            #             0,
            #             rate,
            #             fc=color,
            #             alpha=0.4,
            #             ec="k",
            #             label=f"{binnum}",
            #             hatch="X",
            #         )

            for dph in [-2, -1, 0, 1, 2]:

                ax_efold.errorbar(
                    phase + dph,
                    rate,
                    error,
                    drawstyle="steps-mid",
                    ls="-",
                    alpha=0.8,
                    lw=3,
                    ecolor=color,
                    color=color,
                )

        if legend:
            ax_efold.legend()
        ax_efold.set_xlabel("Phase")
        ax_efold.set_ylabel("Rate")
        ax_efold.set_xlim(-0.1, 2.1)
        if fiducial is not None:
            ax_efold.set_ylim(ymin, ymax)  # type: ignore
        fig.tight_layout()
        if save:
            fig.savefig("efold_check.png")
        return fig, colors

    def check_lightcurve_of_bins(
        self, fiducial: str, lc_files: list, prodpath: str = "phase_resolved"
    ):
        """
        check_lightcurve_of_bins builds  lightcurves for each phase bin and compares with the  lightcurve of the whole observation if necessary.

        Args:
            fiducial (str): light curve fits file of the whole observation
            efolds_files (list): light curve fits files of phase bins to plot
            prodpath (str, optional): path to phase-resolved products. Defaults to 'phase_resolved'.

        This is how to create a plot of RATIO between light curves and fiducial lightcurve
        f0 = fits.open('../spe_and_lc/spe_and_lcAB_sr.lc_bary_orb_corr_nphase_128.efold')
        fig,ax = plt.subplots()

        for binnum in np.arange(1,11):
            fbin = fits.open(f'phase_resolved_bin{binnum}AB_sr.lc_bary_orb_corr_nphase_128.efold')
            ax.plot(f0[1].data['PHASE'], fbin[1].data['RATE1']/f0[1].data['RATE1'], label = binnum)

        ax.legend()


        """
        ind_max = int(1e3)
        os.chdir(self.products_path + "/" + prodpath)

        assert len(glob("*A_bin_*sr.pha")) == len(
            glob("*B_bin_*sr.pha")
        ), "different number of spectra files for mode A and B"
        assert len(glob("*A_bin_*sr.lc")) == len(
            glob("*B_bin_*sr.lc")
        ), "different number of light curve files for mode A and B"

        fig, ax_lc = plt.subplots(figsize=(16, 5))

        if fiducial is not None:
            for fitsfile in [fiducial]:
                light_curve = fits.open(fitsfile)[1].data  # type: ignore
                time = light_curve["TIME"]
                rate = light_curve["RATE"]
                error = light_curve["ERROR"]
                color = "k"
                ax_lc.errorbar(
                    time[: ind_max * (len(lc_files) + 2)],
                    rate[: ind_max * (len(lc_files) + 2)],
                    error[: ind_max * (len(lc_files) + 2)],
                    ls="-.",
                    alpha=0.2,
                    lw=10,
                    label="ref",
                    color=color,
                    ecolor=color,
                )

        else:
            pass

        for fitsfile in lc_files:
            print(f"working with {fitsfile}")
            binnum = fitsfile.split("_")[2]
            light_curve = fits.open(fitsfile)[1].data  # type: ignore
            time = light_curve["TIME"]
            rate = light_curve["RATE"]
            error = light_curve["ERROR"]

            ax_lc.errorbar(
                time[:ind_max],
                rate[:ind_max],
                error[:ind_max],
                ls="-",
                alpha=0.8,
                lw=3,
                label=binnum,
            )

            # color = ax_lc.get_lines()[-1].get_color()
            # ax_lc.fill_between(
            #    time[:ind_max], 0, rate[:ind_max], fc=color, alpha=0.4, ec='k', label=f"{binnum}")

        ax_lc.legend()
        ax_lc.set_xlabel("Time")
        ax_lc.set_ylabel("Rate")
        fig.tight_layout()
        fig.savefig("lc_check.png")
