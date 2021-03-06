import mastcasjobs
from astropy.io import ascii
from astropy.table import Table
from astropy import coordinates, units
import numpy as np
import os
import re

# get the WSID and password if not already defined
import getpass
if not os.environ.get('CASJOBS_WSID'):
    os.environ['CASJOBS_WSID'] = input('Enter Casjobs WSID:')
if not os.environ.get('CASJOBS_PW'):
    os.environ['CASJOBS_PW'] = getpass.getpass('Enter Casjobs password:')


def cone_ps1strm(ra, dec, radius=5, selectcol=['objID', 'raMean', 'decMean', 'z_phot0', 'z_photErr', 'prob_Galaxy', 'prob_Star', 'prob_QSO']):
    """ cone search in pan-starrs strm classification and photo-z
    ra, dec in degrees, radius in arcsec.
    Columns described in https://archive.stsci.edu/hlsp/ps1-strm
    """

    query = """SELECT ps1_strm.*, nearby.distance\nFROM fGetNearbyObjEq({0}, {1}, {2}/60.0) AS nearby\nINNER JOIN catalogRecordRowStore AS ps1_strm\nON ps1_strm.objID = nearby.objID""".format(ra, dec, radius)
  
    jobs = mastcasjobs.MastCasJobs(context="HLSP_PS1_STRM")
    tab = jobs.quick(query, task_name="python cone search")

    return tab[selectcol]


def match_ps1strm(ra, dec, radius, verbose=False):
    """ Compare ra, dec location to the PS1 STRM association.
    """

    co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))

    tab = cone_ps1strm(ra, dec, radius)
    coo = coordinates.SkyCoord(tab['raMean'].tolist(), tab['decMean'].tolist(), unit=(units.deg, units.deg))
    sep = co.separation(coo).to_value(units.arcsec)
    if len(sep):
        i = np.where(sep == sep.min())[0][0]
        if verbose:
            print("Source {0} separated by {1}".format(coo, sep))
        return len(sep), sep[i], tab[i]['objID'], tab[i]['z_phot0'], tab[i]['z_photErr'], tab[i]['prob_Galaxy'], tab[i]['prob_Star'], tab[i]['prob_QSO']
    else:
        return None


def fixcolnames(tab):
    """
    Fix column names returned by the casjobs query
    """

    pat = re.compile(r'\[(?P<name>[^[]+)\]')
    for c in tab.colnames:
        m = pat.match(c)
        if not m:
            raise ValueError("Unable to parse column name '{}'".format(c))
        newname = m.group('name')
        tab.rename_column(c,newname)
    return tab
