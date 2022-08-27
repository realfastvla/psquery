import numpy as np
from matplotlib import pyplot as plt
import os.path

try:
    import sedpy
    from prospect.fitting import fit_model, lnprobfn
    from prospect.likelihood import lnlike_phot, lnlike_spec, write_log
    from prospect.models import priors, sedmodel
    from prospect.models.templates import TemplateLibrary
    from prospect.sources import CSPSpecBasis
    from prospect.utils.obsutils import fix_obs
    from prospect.plotting import sfh
except ImportError:
    print("sedpy or prospect not importing. cannot use sed modeling...")
# new imports
from astropy.coordinates import SkyCoord
from astropy import cosmology
from astroquery.mast import Catalogs
from dustmaps.sfd import SFDQuery
from extinction import fitzpatrick99

from . import irsaquery, psquery, noaoquery

cosmo = cosmology.Planck18

lamd = {
    "ps_g": 4866.0,
    "ps_r": 6215.0,
    "ps_i": 7545.0,
    "ps_z": 8679.0,
    "ps_y": 9633.0,
    "galex_FUV": 1549.0,
    "galex_NUV": 2305.0,
    "wise_w1": 33680.0,
    "wise_w2": 46180.0,
    "wise_w3": 120820.0,
    "wise_w4": 221940.0,
    'decam_g':4770,
    'decam_r':6371.3,
    'decam_z':9157.9
}


def extinct(radec, phot):
    """
    Galactic extinction for given position
    """
    sfd = SFDQuery()
    c = SkyCoord(radec[0], radec[1], unit="deg")
    EBV = sfd(c)
    new_phot = {}
    for i in phot:
        if i in lamd:
            wav = np.array([lamd[i]])
            new_phot[i] = phot[i] - sfd(c) * fitzpatrick99(wav, 3.1)[0]
        else:
            new_phot[i] = phot[i]
    return new_phot


def get_phot(radec, radius, legacy=False, galaxy=False, **kwargs):
    """
    Photometry should be extinction corrected, AB magnitudes
    Input ra, dec and search radius
    Output dictionary with all WISE, PS, and GALEX photometries for the closest source
    """
    # quering WISE and PS1
    wise = irsaquery.cone_wise(radec)
    if galaxy:
        ps1 = mastquery.cone_ps1_casjobs(radec, query='ForcedGalaxyShape')
        try:
            ps1 = ps1[0]
        except IndexError:
            print('No galaxy shape data')
            
    else:
        ps1 = psquery.query_radec(radec, table='stack', phot='Kron')
        ps1 = np.array([float(i) for i in ps1[2][:-1].split(",")])
    
    galex = Catalogs.query_object(
        "{0}  {1}".format(radec[0], radec[1]), radius=radius/3600, catalog="Galex"
    )
    try:
        noao = noaoquery.cone_legacy(radec, radius=1/3600)
    except:
        noao = []
        print('NO NOAO cannot using panstarrs...')

    # init phot with kwargs
    phot = kwargs.copy()

    if len(galex) != 0 and galex[0]["distance_arcmin"] * 60 < radius:
        if (
            type(galex[0]["fuv_mag"]) == np.float64
            and type(galex[0]["fuv_magerr"]) == np.float64
        ):
            phot["galex_FUV"] = galex[0]["fuv_mag"]
            phot["galex_FUV_err"] = galex[0]["fuv_magerr"]
        if (
            type(galex[0]["nuv_mag"]) == np.float64
            and type(galex[0]["nuv_magerr"]) == np.float64
        ):
            phot["galex_NUV"] = galex[0]["nuv_mag"]
            phot["galex_NUV_err"] = galex[0]["nuv_magerr"]
    
    if legacy and len(noao)!=0:
        noao = noao.iloc[0]
        phot['decam_g'] = noao.mag_g
        phot['decam_g_err'] = noao.mag_g_err
        phot['decam_r'] = noao.mag_r
        phot['decam_r_err'] = noao.mag_r_err
        phot['decam_z'] = noao.mag_z
        phot['decam_z_err'] = noao.mag_z_err
    elif galaxy:
        if ps1['gGalMag'] != -999.0 and ps1['gGalMagErr'] != -999.0:
            phot["ps_g"] = ps1['gGalMag']
            phot["ps_g_err"] = ps1['gGalMagErr']
        if ps1['rGalMag'] != -999.0 and ps1['rGalMagErr'] != -999.0:
            phot["ps_r"] = ps1['rGalMag']
            phot["ps_r_err"] = ps1['rGalMagErr']
        if ps1['iGalMag'] != -999.0 and ps1['iGalMagErr'] != -999.0:
            phot["ps_i"] = ps1['iGalMag']
            phot["ps_i_err"] = ps1['iGalMagErr']
        if ps1['zGalMag'] != -999.0 and ps1['zGalMagErr'] != -999.0:
            phot["ps_z"] = ps1['zGalMag']
            phot["ps_z_err"] = ps1['zGalMagErr']
        if ps1['yGalMag'] != -999.0 and ps1['yGalMagErr'] != -999.0:
            phot["ps_y"] = ps1['yGalMag']
            phot["ps_y_err"] = ps1['yGalMagErr']
    else:
        # adding ps1 data
        if ps1[3] != -999.0 and ps1[4] != -999.0:
            phot["ps_g"] = ps1[3]
            phot["ps_g_err"] = ps1[4]
        if ps1[5] != -999.0 and ps1[6] != -999.0:
            phot["ps_r"] = ps1[5]
            phot["ps_r_err"] = ps1[6]
        if ps1[7] != -999.0 and ps1[8] != -999.0:
            phot["ps_i"] = ps1[7]
            phot["ps_i_err"] = ps1[8]
        if ps1[9] != -999.0 and ps1[10] != -999.0:
            phot["ps_z"] = ps1[9]
            phot["ps_z_err"] = ps1[10]
        if ps1[11] != -999.0 and ps1[12] != -999.0:
            phot["ps_y"] = ps1[11]
            phot["ps_y_err"] = ps1[12]


    # adding WISE data (need to have WISE result and WISE position within ps1 radius)
    if len(wise) != 0 and wise[0]["dist"] < radius:
        if (
            type(wise[0]["w1mpro"]) == np.float64
            and type(wise[0]["w1sigmpro"]) == np.float64
        ):
            phot["wise_w1"] = wise[0]["w1mpro"] + 2.699
            phot["wise_w1_err"] = wise[0]["w1sigmpro"]
        if (
            type(wise[0]["w2mpro"]) == np.float64
            and type(wise[0]["w2sigmpro"]) == np.float64
        ):
            phot["wise_w2"] = wise[0]["w2mpro"] + 3.399
            phot["wise_w2_err"] = wise[0]["w2sigmpro"]
        #if (
        #    type(wise[0]["w3mpro"]) == np.float64
        #    and type(wise[0]["w3sigmpro"]) == np.float64
        #):
        #    phot["wise_w3"] = wise[0]["w3mpro"] + 5.174
        #    phot["wise_w3_err"] = wise[0]["w3sigmpro"]
        #if (
        #    type(wise[0]["w4mpro"]) == np.float64
        #    and type(wise[0]["w4sigmpro"]) == np.float64
        #):
        #    phot["wise_w3"] = wise[0]["w3mpro"] + 6.620
        #    phot["wise_w3_err"] = wise[0]["w3sigmpro"]

    phot = extinct(radec=radec, phot=phot)
    return phot


##
# Functions below drawn from prospector examples, modified by Jean Somalwar.
##


def run_fit(phot, hfile="results.h5", emcee=False, plot=True, **params):
    """Wrap the obs/model/sps generation and run for given photometry dict."""

    if "z" not in phot:
        phot["z"] = 0
    if "free_redshift" not in phot:
        phot["free_redshift"] = True
        print("Fitting for redshift...")

    # set default run_params, then overload
    run_params = {
        "object_redshift": phot["z"],
        "free_redshift": phot["free_redshift"],
        "zcontinuous": 1,
    }
    run_params["phot"] = phot
    if "name" in phot:
        run_params["name"] = str(phot["name"]).replace("_", "")
    else:
        run_params["name"] = ""
    run_params["model_template"] = "parasfh"  # default
    run_params["add_dust_emission"] = False  # include dust emission and fit WISE
    run_params["add_AGN_dust"] = False  # add AGN dust?
    run_params["fixed_metallicity"] = None  # fix metals?
    run_params["no_dust"] = False  # fix dust to zero
    run_params["verbose"] = True
    run_params["output_pickles"] = False
    run_params["dynesty"] = False
    run_params["optimize"] = True
    run_params["min_method"] = "lm"
    run_params["nmin"] = 10  # number of draws to start from
    for kk, vv in params:
        run_params[kk] = vv
    run_params["emcee"] = emcee

    # set it up
    obs = (build_obs(**run_params),)
    model = (build_model(**run_params),)
    sps = CSPSpecBasis(
        zcontinuous=1,
        dust_type=2,
        imf_type=1,
        add_neb_emission=run_params["add_dust_emission"],
        compute_vega_mags=False,
    )

    # not clear why these are tuples, but they should be dicts
    if isinstance(obs, tuple):
        obs = obs[0]
    if isinstance(model, tuple):
        model = model[0]

    print("Starting optimization....")
    output = fit_model(obs, model, sps, lnprobfn=lnprobfn, **run_params)
    print("Done optimization in {0:0.1f}s".format(output["optimization"][1]))
    (results, topt) = output["optimization"]

    # Find which of the minimizations gave the best result,
    # and use the parameter vector for that minimization
    ind_best = np.argmin([r.cost for r in results])
    theta_best = results[ind_best].x.copy()
    print(model.theta)
    print(theta_best)

    # quick dict for plotting
    fit_info = {}
    for i, k in enumerate(model.free_params):
        fit_info[k] = theta_best[i]
        if k == "mass":
            fit_info[k] = np.log10(theta_best[i] * sps.ssp.stellar_mass)
        if k == "fage_burst":
            fit_info[k] = (
                1 - theta_best[i]
            )  # convert from age of universe to age of the younger population
        print(k, fit_info[k])

    # generate model
    pspec, pphot, pfrac = model.mean_model(theta_best, obs=obs, sps=sps)

    # dict with latex labels for the plot legend
    theta_lbl_latex = {}
    downup = "${0[0]:0.1f}_{{{0[1]:0.1f}}}^{{{0[2]:0.1f}}}$"
    theta_lbl_latex["mass"] = r"log($M/M_\odot$)=" + downup
    theta_lbl_latex["zred"] = r"z=" + downup
    theta_lbl_latex["logzsol"] = r"$log(Z/Z_\odot$)=" + downup
    theta_lbl_latex["dust2"] = r"E(B-V)=" + downup
    theta_lbl_latex["tau"] = r"$\tau$=" + downup + " Gyr"
    theta_lbl_latex["tage"] = r"$t$=" + downup + " Gyr"
    theta_lbl_latex["fburst"] = r"$M_{{young}}/M$=" + downup
    theta_lbl_latex["fage_burst"] = r"$t_{{young}}/t$=" + downup
    theta_lbl_latex["fagn"] = r"$f_{{AGN}}$=" + downup

    # plot best-fit optimize model
    nice_label = ""
    for k in model.free_params:
        nice_label += theta_lbl_latex[k].format([fit_info[k], 0, 0]) + "  "
        if k == "dust2":
            nice_label += "\n"

    # plot Data, best fit model, and old models
    if "zred" in fit_info.keys():
        zplot = fit_info["zred"]
    else:
        zplot = phot["z"]

    wspec = sps.wavelengths

    if plot:
        plt.clf()
        uplims = obs["islim"]
        plt.errorbar(
            obs["phot_wave"][~uplims],
            obs["mags"][~uplims],
            yerr=obs["mags_unc"][~uplims],
            ecolor="k",
            marker="s",
            ls="",
            lw=1.5,
            alpha=1,
            mew=1.5,
            markerfacecolor="none",
            markeredgecolor="k",
            zorder=5,
        )
        plt.errorbar(
            obs["phot_wave"][uplims],
            obs["mags"][uplims],
            ecolor="k",
            marker="v",
            ls="",
            lw=1.5,
            alpha=1,
            mew=1.5,
            markerfacecolor="none",
            markeredgecolor="k",
            zorder=5,
        )

        line_model = plt.plot(
            wspec * (1 + zplot),
            mi2mg(pspec),
            lw=0.7,
            color="slateblue",
            alpha=0.7,
            label=nice_label,
        )

        plt.xlabel("Wavelength ($\AA$)")
        plt.ylabel("Magnitude (AB)")
        plt.xlim([obs["phot_wave"].min() / 1.2, obs["phot_wave"].max() * 1.3])
        plt.ylim(
            [
                np.nanmax(obs["mags"][np.isfinite(obs["mags"])]) + 3.5,
                np.nanmin(obs["mags"][np.isfinite(obs["mags"])]) - 1,
            ]
        )
        plt.xscale("log")
        plt.annotate(
            run_params["name"],
            (0.03, 0.94),
            xycoords="axes fraction",
            horizontalalignment="left",
        )

        plt.legend(loc=4, fontsize=10)
        plt.xscale("log")
        plt.show()

    if emcee:
        import prospect.io.read_results as reader
        from prospect.io import write_results as writer

        runit = False
        if os.path.exists(hfile):
            resp = input(f"Output file ({hfile}) already exists. \n 'y' to load those results or enter new file name.")
            if resp.lower() != 'y':
                hfile = resp
                print(f"Setting output file name to {hfile}.")
                runit = True
        else:
            runit = True

        run_params["optimize"] = False
        run_params["emcee"] = True
        run_params["nwalkers"] = 100  # Number of emcee walkers
        run_params["niter"] = 1000  # Number of iterations of the MCMC sampling
        run_params["nburn"] = [100]
        run_params["progress"] = False

        if runit:
            # Run emcee
            print("Starting MCMC")
            output = fit_model(obs, model, sps, lnprobfn=lnprobfn, **run_params)  # ,
            print("Done emcee in {0:0.1f}s".format(output["sampling"][1]))

            writer.write_hdf5(hfile, run_params, model, obs,
                  output["sampling"][0], output["optimization"][0],
                  tsample=output["sampling"][1],
                  toptimize=output["optimization"][1])

        # grab results (dictionary), the obs dictionary, and our corresponding models
        # When using parameter files set `dangerous=True`
        result, _obs, _model = reader.results_from(hfile, dangerous=False)

        # get the mean spectrum
        theta_medpos = np.median(result["chain"][:, -1, :], axis=0)
        mspec_medpos, mphot_medpos, _ = model.mean_model(theta_medpos, obs, sps=sps)

        # collect MCMC results
        specphot_obs = []

        run_params = result["run_params"]

        cred_int = 68.1
        fit_info = {}
        fit_info["medpos"] = {}  # median of posterior (default?)
        fit_info["maxp"] = {}  # max of likelihood
        fit_info[
            "weightp"
        ] = {}  # median of samples from poserior wegithed by likelihood

        theta_fit = {}
        theta_fit["medpos"] = np.zeros(
            len(result["theta_labels"])
        )  # to store the median of the posterior
        theta_fit["maxp"] = np.zeros(
            len(result["theta_labels"])
        )  # to store the max likelihood estimate
        theta_fit["weightp"] = np.zeros(
            len(result["theta_labels"])
        )  # to store the weighted mean

        # run FSPS to get best-fit model and surviving mass
        theta = model.theta.copy()
        _, _, _ = model.mean_model(theta, obs, sps=sps)
        sps_pickled = sps.ssp

        # loop again to get uncertainties
        for i, lbl in enumerate(result["theta_labels"]):

            ipostburn = int(-run_params["niter"] / 2)
            th_sample = result["chain"][:, ipostburn:, i].flatten()
            lnp = result["lnprobability"][:, ipostburn:].flatten()

            for ft in theta_fit:

                if lbl == "mass":
                    tt = np.log10(th_sample * sps_pickled.stellar_mass)
                else:
                    tt = th_sample

                if ft == "weightp":
                    w = np.exp(lnp) * th_sample
                else:
                    w = th_sample

                if ft == "maxp":
                    center = tt[np.argmax(lnp)]
                else:
                    center = np.sort(tt)[np.argmin(np.abs(w.cumsum() - w.sum() / 2))]

                theta_fit[ft][i] = center
                islo = np.argmin(
                    np.abs(w.cumsum() - w.sum() * (0.5 - cred_int / 100 / 2.0))
                )
                isup = np.argmin(
                    np.abs(w.cumsum() - w.sum() * (0.5 + cred_int / 100 / 2.0))
                )
                fit_info[ft][lbl] = (
                    center,
                    np.clip(np.abs(center - np.sort(tt)[islo]), 0, 999),
                    np.clip(np.abs(center - np.sort(tt)[isup]), 0, 999),
                )

        mphot_dict = {
            k: np.zeros(result["run_params"]["nwalkers"]) for k in obs["filternames"]
        }

        # loop over final state of all walkers to sample some spectra
        for i in range(result["run_params"]["nwalkers"]):
            theta = result["chain"][i, -1, :]
            mspec, mphot, _ = model.mean_model(theta, obs, sps=sps)
            specphot_obs.append((mspec, mphot))

            for l, k in enumerate(obs["filternames"]):
                mphot_dict[k][i] = mi2mg(mphot[l])

        if "zred" in fit_info["medpos"].keys():
            zplot = fit_info["medpos"]["zred"][0]
        else:
            zplot = phot["z"]

        # -----
        # plot some models
        colormap = plt.cm.viridis

        for i in range(result["run_params"]["nwalkers"]):
            mspec, mphot = specphot_obs[i]
            plt.plot(
                sps_pickled.wavelengths * (1 + zplot),
                mi2mg(mspec),
                lw=0.7,
                color=colormap(0.2),
                alpha=0.1,
                zorder=1,
            )
            plt.plot(
                obs["phot_wave"],
                mi2mg(mphot),
                label="",
                marker="o",
                alpha=0.1,
                ls="",
                lw=1,
                zorder=2,
                markersize=5,
                markerfacecolor="none",
                markeredgecolor=colormap(0.5),
                markeredgewidth=0.7,
            )

        nice_label = ""
        for k in model.free_params:
            nice_label += theta_lbl_latex[k].format(fit_info["medpos"][k]) + "  "
            if k == "dust2":
                nice_label += "\n"

        plt.plot(
            wspec * (1 + zplot),
            mi2mg(mspec_medpos) * 0,
            label=nice_label,
            lw=0.7,
            color=colormap(0.2),
            alpha=0.9,
        )
        plt.plot(
            wspec * (1 + zplot), mi2mg(mspec_medpos), lw=0.7, color="red", alpha=0.9
        )

        plt.legend(loc=4, fontsize=10, framealpha=0.7, facecolor="white")

        plt.xlim(obs["phot_wave"].min() / 1.1, obs["phot_wave"].max() * 1.2)
        plt.ylim(25, 15)
        plt.xscale("log")
        plt.show()

        cornerfig = reader.subcorner(
            result,
            start=0,
            thin=5,
            fig=plt.subplots(len(theta), len(theta), figsize=(10, 10))[0],
        )
        plt.show()
        return (wspec * (1 + zplot), mi2mg(mspec_medpos)), (obs['phot_wave'], obs['mags']), 

    #return pspec, pphot, pfrac
    return (wspec * (1 + zplot), mi2mg(pspec)), (obs['phot_wave'], obs['mags'], mi2mg(pphot)), theta_best



def build_obs(
    phot=None,
    filternames=None,
    mag_err_clip=0.05,
    standard=["galex", "ps", "sdss", "wise", "decam"],
    **extras
):
    """Build a dictionary of observational data.

    phot = table or dictionary with photometry. Should already be extinction corrected.
    Upper limits have flux=0 and "err" equal to the limit (in what units?).
    filternames is list of bands and must match name in sedpy.
    If None, names inferred from phot using standard root names.

    returns obs: A dictionary of observational data to use in the fit.
    """

    obs = {}
    if filternames is None:
        sel = lambda x: any([(st in x.lower()) for st in standard if "err" not in x])
        filternames = list(filter(sel, phot.keys()))
    flt_use = np.array([], dtype="S20")
    data_use, edata_use = np.array([]), np.array([])
    for k in filternames:
        if phot[k] != 0:
            flt_use = np.append(flt_use, k)
            data_use = np.append(data_use, phot[k])
            edata_use = np.append(edata_use, phot[k + "_err"])

    edata_use = np.array(edata_use)
    edata_use = np.clip(edata_use, mag_err_clip, 10)
    obs["islim"] = edata_use <= 0

    obs["filternames"] = filternames
    obs["filters"] = sedpy.observate.load_filters(filternames)
    obs["mags"] = data_use
    obs["mags_unc"] = edata_use
    obs["maggies"] = 10 ** (-0.4 * obs["mags"])  # fluxes in units of maggies (Jy/3631)
    obs["maggies_unc"] = np.log(10) / 2.5 * obs["mags_unc"] * obs["maggies"]

    # deal with upper limits
    # assume upper limits are 5-sigma (true for GALEX)
    nsigma = np.repeat(5, len(obs["filternames"]))
    iswise = [i for i, w in enumerate(filternames) if "wise" in w]
    nsigma[iswise] = 2  # WISE limits are 95%CL

    obs["maggies_unc"][obs["islim"]] = (
        10 ** (-0.4 * obs["mags"][obs["islim"]]) / nsigma[obs["islim"]]
    )
    obs["maggies"][obs["islim"]] = 0

    # Photometry mask, True = include in fit
    obs["phot_mask"] = np.array([True for f in obs["filters"]])

    # This is an array of effective wavelengths for each of the filters.
    # It is not necessary, but it can be useful for plotting so we store it here as a convenience
    obs["phot_wave"] = np.array([f.wave_effective for f in obs["filters"]])

    # We do not have a spectrum, so we set some required elements of the obs dictionary to None.
    obs["wavelength"] = None
    obs["spectrum"] = None
    obs["unc"] = None
    obs["mask"] = None

    # This function ensures all required keys are present in the obs dictionary,
    # adding default values if necessary
    obs = fix_obs(obs)

    return obs


def build_model(
    name="",
    results_dir="./results/",
    object_redshift=None,
    free_redshift=False,
    add_dust_emission=False,
    no_dust=False,
    model_template="parasfh",
    **extras
):

    """Instantiate and return a ProspectorParams model subclass.
    :param fixed_metallicity: (optional, default: None)
        If given, fix the model metallicity (:math:`log(Z/Z_sun)`) to the given value.
    :param add_dust_emission: (optional, default: False)
        If `True`, add dust emission and associated (fixed) parameters to the model.
    :returns model:
        An instance of prospect.models.SedModel
    """

    # Get (a copy of) one of the prepackaged model set dictionaries.
    # This is a dictionary of dictionaries, keyed by parameter name
#    model_params = TemplateLibrary["ssp"]
    model_params = TemplateLibrary["parametric_sfh"]
    model_params["sfh"]["init"] = 1  # this is 4 in template
    model_params["tburst"] = {
        "N": 1,
        "isfree": False,
        "init": 0,
        "units": """the age of the Universe when the burst occurs (Gyr)""",
    }

    # important: overwrite some default from the template
    model_params["dust_type"]["init"] = 2  # this is 0 in template
    model_params["imf_type"]["init"] = 1  # this seems to be 1 in template

    # Some initial values, should be adjusted probably
    model_params["dust2"]["init"] = 0.05
    model_params["logzsol"]["init"] = -0.63
    model_params["tage"]["init"] = 4.67
    model_params["mass"]["init"] = 10 ** (9)
    model_params["tage"]["prior"] = priors.TopHat(mini=0.1, maxi=12.5)
    model_params["dust2"]["prior"] = priors.TopHat(mini=0.0, maxi=1)
    model_params["mass"]["prior"] = priors.LogUniform(mini=1e7, maxi=1e12)

    if "tau" in model_params:
        model_params["tau"]["prior"] = priors.LogUniform(mini=0.1, maxi=1)

    # If we are going to be using emcee, it is useful to provide a
    # minimum scale for the cloud of walkers (the default is 0.1)
    model_params["mass"]["init_disp"] = 1e8
    model_params["mass"]["disp_floor"] = 1e7
    model_params["tage"]["init_disp"] = 3
    model_params["tage"]["disp_floor"] = 2
    model_params["logzsol"]["init_disp"] = 1
    model_params["logzsol"]["disp_floor"] = 0.5
    model_params["dust2"]["init_disp"] = 1
    model_params["dust2"]["disp_floor"] = 0.5

    if "tau" in model_params:
        model_params["tau"]["init_disp"] = 1
        model_params["tau"]["disp_floor"] = 1

    model_params["fagn"] = {
        "N": 1,
        "isfree": False,
        "init": 0.01,
        "units": "",
        "prior": priors.TopHat(mini=0.0, maxi=2),
    }

    # ensure correct redshift
    model_params["zred"]["isfree"] = free_redshift
    model_params["zred"]["init"] = object_redshift

    if add_dust_emission:
        # Add dust emission (with fixed dust SED parameters)
        # Since `model_params` is a dictionary of parameter specifications,
        # and `TemplateLibrary` returns dictionaries of parameter specifications,
        # we can just update `model_params` with the parameters described in the
        # pre-packaged `dust_emission` parameter set.
        model_params.update(TemplateLibrary["dust_emission"])

    # Now instantiate the model object using this dictionary of parameter specifications
    model = sedmodel.SedModel(model_params)

    return model


def mi2mg(maggies):
    return -2.5 * np.log10(maggies)

def read_h5(hfile, plot=True):
    # grab results (dictionary), the obs dictionary, and our corresponding models
    # When using parameter files set `dangerous=True`
    import prospect.io.read_results as reader

    result, obs, model = reader.results_from(hfile, dangerous=False)
    run_params = result['run_params']
    model = (build_model(**run_params),)[0]
    
    sps = CSPSpecBasis(zcontinuous=1,dust_type=2, imf_type=1, 
                       add_neb_emission=run_params["add_dust_emission"],
                       compute_vega_mags=False,
    )
    phot = run_params["phot"]
    
    # get the mean spectrum
    theta_medpos = np.median(result["chain"][:, -1, :], axis=0)
    mspec_medpos, mphot_medpos, _ = model.mean_model(theta_medpos, obs, sps=sps)

    # collect MCMC results
    specphot_obs = []

    cred_int = 68.1
    fit_info = {}
    fit_info["medpos"] = {}  # median of posterior (default?)
    fit_info["maxp"] = {}  # max of likelihood
    fit_info[
        "weightp"
    ] = {}  # median of samples from poserior wegithed by likelihood

    theta_fit = {}
    theta_fit["medpos"] = np.zeros(
        len(result["theta_labels"])
    )  # to store the median of the posterior
    theta_fit["maxp"] = np.zeros(
        len(result["theta_labels"])
    )  # to store the max likelihood estimate
    theta_fit["weightp"] = np.zeros(
        len(result["theta_labels"])
    )  # to store the weighted mean

    # run FSPS to get best-fit model and surviving mass
    theta = model.theta.copy()
    _, _, _ = model.mean_model(theta, obs, sps=sps)
    sps_pickled = sps.ssp

    # loop again to get uncertainties
    for i, lbl in enumerate(result["theta_labels"]):

        ipostburn = int(-run_params["niter"] / 2)
        th_sample = result["chain"][:, ipostburn:, i].flatten()
        lnp = result["lnprobability"][:, ipostburn:].flatten()

        for ft in theta_fit:

            if lbl == "mass":
                tt = np.log10(th_sample * sps_pickled.stellar_mass)
            else:
                tt = th_sample

            if ft == "weightp":
                w = np.exp(lnp) * th_sample
            else:
                w = th_sample

            if ft == "maxp":
                center = tt[np.argmax(lnp)]
            else:
                center = np.sort(tt)[np.argmin(np.abs(w.cumsum() - w.sum() / 2))]

            theta_fit[ft][i] = center
            islo = np.argmin(
                np.abs(w.cumsum() - w.sum() * (0.5 - cred_int / 100 / 2.0))
            )
            isup = np.argmin(
                np.abs(w.cumsum() - w.sum() * (0.5 + cred_int / 100 / 2.0))
            )
            fit_info[ft][lbl] = (
                center,
                np.clip(np.abs(center - np.sort(tt)[islo]), 0, 999),
                np.clip(np.abs(center - np.sort(tt)[isup]), 0, 999),
            )

    mphot_dict = {
        k: np.zeros(result["run_params"]["nwalkers"]) for k in obs["filternames"]
    }

    # loop over final state of all walkers to sample some spectra
    for i in range(result["run_params"]["nwalkers"]):
        theta = result["chain"][i, -1, :]
        mspec, mphot, _ = model.mean_model(theta, obs, sps=sps)
        specphot_obs.append((mspec, mphot))

        for l, k in enumerate(obs["filternames"]):
            mphot_dict[k][i] = mi2mg(mphot[l])
          
    if "zred" in fit_info["medpos"].keys():
        zplot = fit_info["medpos"]["zred"][0]
    else:
        zplot = phot["z"]
    
    wspec = sps.wavelengths

    # calculate SFR
    sfr = sfh.parametric_sfr(times=cosmo.lookback_time(phot["z"]).to('Gyr').value, sfh=1,
                             mass=10**fit_info["medpos"]["mass"][0],
                             tage=fit_info["medpos"]["tage"][0],
                             tau=fit_info["medpos"]["tau"][0])[0]
    fit_info["medpos"]["sfr"] = [sfr]  # TODO: add bounds
    
    if plot:
        # dict with latex labels for the plot legend
        theta_lbl_latex = {}
        downup = "${0[0]:0.1f}_{{{0[1]:0.1f}}}^{{{0[2]:0.1f}}}$"
        theta_lbl_latex["mass"] = r"log($M/M_\odot$)=" + downup
        theta_lbl_latex["zred"] = r"z=" + downup
        theta_lbl_latex["logzsol"] = r"$log(Z/Z_\odot$)=" + downup
        theta_lbl_latex["dust2"] = r"E(B-V)=" + downup
        theta_lbl_latex["tau"] = r"$\tau$=" + downup + " Gyr"
        theta_lbl_latex["tage"] = r"$t$=" + downup + " Gyr"
        theta_lbl_latex["fburst"] = r"$M_{{young}}/M$=" + downup
        theta_lbl_latex["fage_burst"] = r"$t_{{young}}/t$=" + downup
        theta_lbl_latex["fagn"] = r"$f_{{AGN}}$=" + downup
        
        # plot some models
        colormap = plt.cm.viridis

        for i in range(result["run_params"]["nwalkers"]):
            mspec, mphot = specphot_obs[i]
            plt.plot(
                sps_pickled.wavelengths * (1 + zplot),
                mi2mg(mspec),
                lw=0.7,
                color=colormap(0.2),
                alpha=0.1,
                zorder=1,
            )
            plt.plot(
                obs["phot_wave"],
                mi2mg(mphot),
                label="",
                marker="o",
                alpha=0.1,
                ls="",
                lw=1,
                zorder=2,
                markersize=5,
                markerfacecolor="none",
                markeredgecolor=colormap(0.5),
                markeredgewidth=0.7,
            )
            
        plt.plot(
            obs["phot_wave"],
            obs['mags'],
            label="Observations",
            marker="o",
            alpha=1,
            ls="",
            lw=1,
            zorder=2,
            markersize=5,
            markerfacecolor=colormap(0.8),
            markeredgecolor=colormap(0.8),
            markeredgewidth=0.7,
        )

        nice_label = ""
        for k in model.free_params:
            nice_label += theta_lbl_latex[k].format(fit_info["medpos"][k]) + "  "
            if k == "dust2":
                nice_label += "\n"

        plt.plot(
            wspec * (1 + zplot),
            mi2mg(mspec_medpos) * 0,
            label=nice_label,
            lw=0.7,
            color=colormap(0.2),
            alpha=0.9,
        )
        plt.plot(
            wspec * (1 + zplot), mi2mg(mspec_medpos), lw=0.7, color="red", alpha=0.9
        )

        plt.legend(loc=4, fontsize=10, framealpha=0.7, facecolor="white")

        plt.xlim(obs["phot_wave"].min() / 1.1, obs["phot_wave"].max() * 1.2)
        plt.ylim(
            [
                np.nanmax(obs["mags"][np.isfinite(obs["mags"])]) + 3.5,
                np.nanmin(obs["mags"][np.isfinite(obs["mags"])]) - 1,
            ]
        )
        plt.xscale("log")
        plt.show()

        cornerfig = reader.subcorner(
            result,
            start=0,
            thin=5,
            fig=plt.subplots(len(theta), len(theta), figsize=(10, 10))[0],
        )
        plt.show()
    
    return (wspec * (1 + zplot), mi2mg(mspec_medpos)), (obs['phot_wave'], obs['mags'], mi2mg(mphot)), fit_info['medpos']


def run_fit2(phot, hfile="results.h5", emcee=True, **params):
    """Just runs MCMC and saves a file
        NO ANALYSIS
    """

    # set default run_params, then overload
    run_params = {
        "object_redshift": phot["z"],
        "free_redshift": phot["free_redshift"],
        "zcontinuous": 1,
    }
    run_params["phot"] = phot
    if "name" in phot:
        run_params["name"] = str(phot["name"]).replace("_", "")
    else:
        run_params["name"] = ""
    run_params["model_template"] = "parasfh"  # default
    run_params["add_dust_emission"] = False  # include dust emission and fit WISE
    run_params["add_AGN_dust"] = False  # add AGN dust?
    run_params["fixed_metallicity"] = None  # fix metals?
    run_params["no_dust"] = False  # fix dust to zero
    run_params["verbose"] = True
    run_params["output_pickles"] = False
    run_params["dynesty"] = False
    run_params["optimize"] = True
    run_params["min_method"] = "lm"
    run_params["nmin"] = 10  # number of draws to start from
    for kk, vv in params:
        run_params[kk] = vv
    run_params["emcee"] = emcee

    # set it up
    obs = (build_obs(**run_params),)
    model = (build_model(**run_params),)
    sps = CSPSpecBasis(
        zcontinuous=1,
        dust_type=2,
        imf_type=1,
        add_neb_emission=run_params["add_dust_emission"],
        compute_vega_mags=False,
    )

    # not clear why these are tuples, but they should be dicts
    if isinstance(obs, tuple):
        obs = obs[0]
    if isinstance(model, tuple):
        model = model[0]

    if emcee:
        import prospect.io.read_results as reader
        from prospect.io import write_results as writer

        runit = False
        if os.path.exists(hfile):
            resp = input(f"Output file ({hfile}) already exists. \n 'y' to load those results or enter new file name.")
            if resp.lower() != 'y':
                hfile = resp
                print(f"Setting output file name to {hfile}.")
                runit = True
        else:
            runit = True

        run_params["optimize"] = False
        run_params["emcee"] = True
        run_params["nwalkers"] = 100  # Number of emcee walkers
        run_params["niter"] = 1000  # Number of iterations of the MCMC sampling
        run_params["nburn"] = [100]
        run_params["progress"] = False

        if runit:
            # Run emcee
            print("Starting MCMC")
            output = fit_model(obs, model, sps, lnprobfn=lnprobfn, **run_params)  # ,
            print("Done emcee in {0:0.1f}s".format(output["sampling"][1]))
            
            writer.write_hdf5(hfile, run_params, model, obs,
                  output["sampling"][0], output["optimization"][0],
                  tsample=output["sampling"][1],
                  toptimize=output["optimization"][1])
