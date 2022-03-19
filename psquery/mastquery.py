try:
    import mastcasjobs
except ImportError:
    print("mastcasjobs not available. Cannot use mastquery.")

# get the WSID and password if not already defined
import getpass
import os
import re

import numpy as np
from astropy import coordinates, units
from astropy.io import ascii
from astropy.table import Table

from . import get_coord

if not os.environ.get("CASJOBS_WSID"):
    os.environ["CASJOBS_WSID"] = input("Enter Casjobs WSID:")
if not os.environ.get("CASJOBS_PW"):
    os.environ["CASJOBS_PW"] = getpass.getpass("Enter Casjobs password:")


def cone_ps1strm(
    radec,
    radius=5/3600,
    selectcol=[
        "objID",
        "raMean",
        "decMean",
        "z_phot0",
        "z_photErr",
        "prob_Galaxy",
        "prob_Star",
        "prob_QSO",
    ],
):
    """cone search in pan-starrs strm classification and photo-z
    ra, dec in degrees, radius in degrees.
    Columns described in https://archive.stsci.edu/hlsp/ps1-strm
    """

    ra, dec = get_coord(radec, ret="radec")

    query = """SELECT ps1_strm.*, nearby.distance\nFROM fGetNearbyObjEq({0}, {1}, {2}) AS nearby\nINNER JOIN catalogRecordRowStore AS ps1_strm\nON ps1_strm.objID = nearby.objID""".format(
        ra, dec, radius*60
    )

    jobs = mastcasjobs.MastCasJobs(context="HLSP_PS1_STRM")
    tab = jobs.quick(query, task_name="python ps1strm cone search")

    return tab[selectcol]


def match_ps1strm(radec, radius, verbose=False):
    """
    Compare ra, dec location to the PS1 STRM association.
    ra, dec in degrees, radius in degrees.
    """

    ra, dec = get_coord(radec, ret="radec")
    co = get_coord(radec, ret="skycoord")

    tab = cone_ps1strm(ra, dec, radius)
    coo = coordinates.SkyCoord(
        tab["raMean"].tolist(), tab["decMean"].tolist(), unit=(units.deg, units.deg)
    )
    sep = co.separation(coo).to_value(units.arcsec)
    if len(sep):
        i = np.where(sep == sep.min())[0][0]
        if verbose:
            print("Source {0} separated by {1}".format(coo, sep))
        return (
            len(sep),
            sep[i],
            tab[i]["objID"],
            tab[i]["z_phot0"],
            tab[i]["z_photErr"],
            tab[i]["prob_Galaxy"],
            tab[i]["prob_Star"],
            tab[i]["prob_QSO"],
        )
    else:
        return None


def fixcolnames(tab):
    """
    Fix column names returned by the casjobs query
    """

    pat = re.compile(r"\[(?P<name>[^[]+)\]")
    for c in tab.colnames:
        m = pat.match(c)
        if not m:
            raise ValueError("Unable to parse column name '{}'".format(c))
        newname = m.group("name")
        tab.rename_column(c, newname)
    return tab


def cone_emline(
    ra,
    dec,
    radius=5/3600,
    selectcol=[
        "specObjID",
        "ra",
        "dec",
        "z",
        "zErr",
        "bpt",
        "Flux_Ha_6562",
        "Flux_NII_6583",
        "Flux_Hb_4861",
        "Flux_OIII_5006",
    ],
):
    """box search in emissionLinesPort table
    ra, dec in degrees, radius in degrees.
    Columns described in http://skyserver.sdss.org/dr16/en/help/browser/browser.aspx?cmd=description+emissionLinesPort+U#&&history=description+emissionLinesPort+U
    TODO: use selectcol in sql query
    """

    # dumb way
    #    query = """SELECT TOP 10 emline.*\nFROM emissionLinesPort AS emline\nWHERE ra > {0} and ra < {1}\nAND dec > {2} and dec < {3}""".format(ra-size/2, ra+size/2, dec-size/2, dec+size/2)
    colstr = ", "
    query = """SELECT TOP 10 G.specobjID, G.ra, G.dec, G.z, G.bpt, G.Flux_Ha_6562, G.Flux_NII_6583, G.Flux_Hb_4861, G.Flux_OIII_5006, N.distance\nFROM emissionLinesPort as G\nJOIN dbo.fGetNearbySpecObjEq({0}, {1}, {2}) AS N\nON G.specobjID = N.specobjID""".format(
        ra, dec, radius * 60
    )

    jobs = mastcasjobs.MastCasJobs(context="SDSSDR14")
    tab = jobs.quick(query, task_name="python emission line cone search")

    return tab


def cone_galaxymass(ra, dec, radius=5/3600, selectcol=[]):
    """Query stellarMassStarFormingPort table to get mass, star formation rate and other galaxy properties with BOSS redshifts.
    Based on Maraston et al (2006).
    ra, dec in degrees, size in degrees.
    """

    colstr = ", "
    query = """SELECT TOP 10 G.specobjID, G.ra, G.dec, G.z, G.zErr, G.logMass, G.SFR, G.Metallicity, G.age, N.distance\nFROM stellarMassStarFormingPort as G\nJOIN dbo.fGetNearbySpecObjEq({0}, {1}, {2}) AS N\nON G.specobjID = N.specobjID""".format(
        ra, dec, radius * 60
    )

    jobs = mastcasjobs.MastCasJobs(context="SDSSDR14")
    tab = jobs.quick(query, task_name="python galaxy mass port cone search")

    return tab


def cone_ps1_casjobs(radec, radius=5/3600, ndet=1, nr=1, query="MeanObject"):
    """cone search in ps1 via casjobs (similar to mastquery PS1STRM function)
    ra, dec in degrees, radius in degrees.
    ndet, nr define number of total and r band detections required.
    Gets primary detection from stacks for rkronRad.
    """

    ra, dec = get_coord(radec, ret="radec")

    if "select" in query:
        pass
    elif query == "ForcedGalaxyShape":
        query = f"""select o.objID, o.raMean, o.decMean, o.nDetections, o.ng, o.nr, o.ni, o.nz, o.ny, m.gGalMag, m.gGalMagErr, m.rGalMag, m.rGalMagErr, m.iGalMag, m.iGalMagErr, m.zGalMag, m.zGalMagErr, m.yGalMag, m.yGalMagErr\nfrom fGetNearbyObjEq({ra}, {dec}, {radius}) nb\ninner join ObjectThin o on o.objid = nb.objid and nDetections>{ndet} and nr>{nr}\ninner join ForcedGalaxyShape m on o.objid = m.objid\ninner join StackObjectAttributes d on o.objid = d.objid and d.primaryDetection = 1""".format(
            ra, dec, radius*60, ndet, nr
        )
    elif query == "MeanObject":
        query = f"""select o.objID, o.raMean, o.decMean, o.nDetections, o.ng, o.nr, o.ni, o.nz, o.ny, m.gMeanPSFMag, m.rMeanPSFMag, m.iMeanPSFMag, m.zMeanPSFMag, m.yMeanPSFMag, m.rMeanKronMag, d.rkronRad\nfrom fGetNearbyObjEq({ra}, {dec}, {radius}) nb\ninner join ObjectThin o on o.objid = nb.objid and nDetections>{ndet} and nr>{nr}\ninner join MeanObject m on o.objid = m.objid\ninner join StackObjectAttributes d on o.objid = d.objid and d.primaryDetection = 1""".format(
            ra, dec, radius*60, ndet, nr
        )
    elif query == "StackPetrosian":
        query = f"""select o.objID, o.raMean, o.decMean, o.nDetections, o.ng, o.nr, o.ni, o.nz, o.ny, m.gpetMag, m.gpetMagErr, m.rpetMag, m.rpetMagErr, m.ipetMag, m.ipetMagErr, m.zpetMag, m.zpetMagErr, m.ypetMag, m.ypetMagErr\nfrom fGetNearbyObjEq({ra}, {dec}, {radius}) nb\ninner join ObjectThin o on o.objid = nb.objid and nDetections>{ndet} and nr>{nr}\ninner join StackPetrosian m on o.objid = m.objid\ninner join StackObjectAttributes d on o.objid = d.objid and d.primaryDetection = 1""".format(
            ra, dec, radius*60, ndet, nr
        )

    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    tab = jobs.quick(query, task_name="python ps1 DR2 cone search")

    return tab


def cone_ps1_psc(radec, radius=5/3600):
    """cone search in ps1 PSC via casjobs (similar to mastquery PS1STRM function)
    ra, dec in degrees, radius in degres.
    See Tachibana & Miller (2018; https://iopscience.iop.org/article/10.1088/1538-3873/aae3d9).
    Optimal extended source has ps_score<0.83.
    """

    ra, dec = get_coord(radec, ret="radec")

    query = f"""select p.objID, p.ps_score, p.raMean, p.decMean\nfrom pointsource_magnitudes_view as p\ninner join fGetNearbyObjEq({ra}, {dec}, {radius}) as r on p.objid=r.objid and p.primaryDetection = 1""".format(
        ra, dec, radius*60
    )

    jobs = mastcasjobs.MastCasJobs(context="HLSP_PS1_PSC")
    tab = jobs.quick(query, task_name="python psc cone search")

    return tab
