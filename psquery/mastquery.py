import mastcasjobs
from astropy.io import ascii
from astropy.table import Table
from astropy import coordinates, units

import sys
import os
import re
import numpy as np
import pylab
import json

try: # Python 3.x
    from urllib.parse import quote as urlencode
    from urllib.request import urlretrieve
except ImportError:  # Python 2.x
    from urllib import pathname2url as urlencode
    from urllib import urlretrieve

try: # Python 3.x
    import http.client as httplib 
except ImportError:  # Python 2.x
    import httplib

# get the WSID and password if not already defined
import getpass
if not os.environ.get('CASJOBS_WSID'):
    os.environ['CASJOBS_WSID'] = input('Enter Casjobs WSID:')
if not os.environ.get('CASJOBS_PW'):
    os.environ['CASJOBS_PW'] = getpass.getpass('Enter Casjobs password:')


def cone_ps1strm(ra, dec, radius=5/3600):
    """ cone search in pan-starrs dr2
    ra, dec in degrees, radius in arcsec.
    """

    query = """SELECT ps1_strm.*, nearby.distance\nFROM fGetNearbyObjEq({0}, {1}, {2}/60.0) AS nearby\nINNER JOIN catalogRecordRowStore AS ps1_strm\nON ps1_strm.objID = nearby.objID""".format(ra, dec, radius)
  
    jobs = mastcasjobs.MastCasJobs(context="HLSP_PS1_STRM")
    results = jobs.quick(query, task_name="python cone search")

    # ref coord
    co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))
    coo = coordinates.SkyCoord(results['raMean'].tolist(), results['decMean'].tolist(), unit=(units.deg, units.deg))
    sep = co.separation(coo).to_value(units.arcsec)
    print("Source {0} separated by {1}".format(coo, sep))
    return zip(sep, results['objID'], results['z_phot'], results['z_photErr'])


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
