import numpy as np

import requests 
from urllib.parse import quote as urlencode 
from urllib.request import urlretrieve 
import http.client as httplib
     
from astropy.io import ascii 
from astropy.table import Table 
from astropy import coordinates, units
import mastcasjobs

def query_radec(ra, dec, ndet=5, radius=5/3600, columns=['objID', 'raMean', 'decMean', 'gMeanKronMag', 'gMeanKronMagErr',
                                                         'rMeanKronMag', 'rMeanKronMagErr', 'iMeanKronMag', 'iMeanKronMagErr',
                                                         'zMeanKronMag', 'zMeanKronMagErr', 'yMeanKronMag', 'yMeanKronMagErr'],
                verbose=False, release='dr2'):
    """ cone search in pan-starrs dr2.
    radius is in degrees.
    Returns number of tuple with number of matches, separation in arcsec to nearest, and photometry of nearest.
    """

    if ndet:
        constraints = {'nDetections.gt': ndet}
    else:
        constraints = {}
    radius = radius

    co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))

    if isinstance(columns, str):
        columns = columns.split(",")
    elif isinstance(columns, list):
        pass
    else:
        print('columns must be list or comma-delimited string')
#    columns = [x.strip() for x in columnstr if x and not x.startswith('#')]
#    add "nDetections,ng,nr,ni,nz,ny"?

    results = ps1cone(ra, dec, radius, release=release, columns=columns, verbose=verbose, **constraints)
    lines = results.split('\n')
    if len(lines) == 3:
        if verbose:
            print('Found one')
        rao, deco = lines[1].split(',')[1:3]
        coo = coordinates.SkyCoord(float(rao), float(deco), unit=(units.deg, units.deg))
        sep = co.separation(coo).to_value(units.arcsec)
        if verbose:
            print("Source {0} separated by {1}".format(coo, sep))
        return 1, sep, lines[1]
    elif len(lines) > 3:
        if verbose:
            print('Found multiple:')
            print(lines)
#        coos = []
        line_min = ''
        sep_min = radius*3600
        for line in lines[1:]:
            if line:
                rao, deco = line.split(',')[1:3]
                coo = coordinates.SkyCoord(float(rao), float(deco), unit=(units.deg, units.deg))
                sep = co.separation(coo).to_value(units.arcsec)
                if sep < sep_min:
                    line_min = line
                    sep_min = sep
#                coos.append(coordinates.SkyCoord(float(rao), float(deco), unit=(units.deg, units.deg)))
#        return coos
        return len(lines)-2, sep_min, line_min
    else:
        if verbose:
            print('Nothing there.')
        return None

def cone_ps1_casjobs(ra, dec, radius=5, ndet=1, nr=1):
    """ cone search in ps1 via casjobs (similar to mastquery PS1STRM function)
    ra, dec in degrees, radius in arcsec.
    ndet, nr define number of total and r band detections required.
    Gets primary detection from stacks for rkronRad.
    """

    query = f"""select o.objID, o.raMean, o.decMean, o.nDetections, o.ng, o.nr, o.ni, o.nz, o.ny, m.gMeanPSFMag, m.rMeanPSFMag, m.iMeanPSFMag, m.zMeanPSFMag, m.yMeanPSFMag, m.rMeanKronMag, d.rkronRad\nfrom fGetNearbyObjEq({ra}, {dec}, {radius}/60.0) nb\ninner join ObjectThin o on o.objid = nb.objid and nDetections>{ndet} and nr>{nr}\ninner join MeanObject m on o.objid = m.objid\ninner join StackObjectAttributes d on o.objid = d.objid and d.primaryDetection = 1""".format(ra, dec, radius, ndet, nr)
  
    jobs = mastcasjobs.MastCasJobs(context="PanSTARRS_DR2")
    tab = jobs.quick(query, task_name="python ps1 DR2 cone search")

    return tab


def cone_ps1_psc(ra, dec, radius=5):
    """ cone search in ps1 PSC via casjobs (similar to mastquery PS1STRM function)
    ra, dec in degrees, radius in arcsec.
    See Tachibana & Miller (2018; https://iopscience.iop.org/article/10.1088/1538-3873/aae3d9).
    Optimal extended source has ps_score<0.83.
    """

    query = f"""select p.objID, p.ps_score, p.raMean, p.decMean\nfrom pointsource_magnitudes_view as p\ninner join fGetNearbyObjEq({ra}, {dec}, {radius}/60.) as r on p.objid=r.objid and p.primaryDetection = 1""".format(ra, dec, radius)
  
    jobs = mastcasjobs.MastCasJobs(context="HLSP_PS1_PSC")
    tab = jobs.quick(query, task_name="python psc cone search")

    return tab


def ps1cone(ra,dec,radius,table="mean",release="dr2",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a cone search of the PS1 catalog
    
    Parameters
    ----------
    ra (float): (degrees) J2000 Right Ascension
    dec (float): (degrees) J2000 Declination
    radius (float): (degrees) Search radius (<= 0.5 degrees)
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2)
    """
    
    data = kw.copy()
    data['ra'] = ra
    data['dec'] = dec
    data['radius'] = radius
    return ps1search(table=table,release=release,format=format,columns=columns,
                    baseurl=baseurl, verbose=verbose, **data)


def ps1search(table="mean",release="dr2",format="csv",columns=None,
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs", verbose=False,
           **kw):
    """Do a general search of the PS1 catalog (possibly without ra/dec/radius)
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    format: csv, votable, json
    columns: list of column names to include (None means use defaults)
    baseurl: base URL for the request
    verbose: print info about request
    **kw: other parameters (e.g., 'nDetections.min':2).  Note this is required!
    """
    
    data = kw.copy()
    if not data:
        raise ValueError("You must specify some parameters for search")
    checklegal(table,release)
    if format not in ("csv","votable","json"):
        raise ValueError("Bad value for format")
    url = "{baseurl}/{release}/{table}.{format}".format(**locals())
    if columns:
        # check that column values are legal
        # create a dictionary to speed this up
        dcols = {}
        for col in ps1metadata(table,release)['name']:
            dcols[col.lower()] = 1
        badcols = []
        for col in columns:
            if col.lower().strip() not in dcols:
                badcols.append(col)
        if badcols:
            raise ValueError('Some columns not found in table: {}'.format(', '.join(badcols)))
        # two different ways to specify a list of column values in the API
        # data['columns'] = columns
        data['columns'] = '[{}]'.format(','.join(columns))

# either get or post works
#    r = requests.post(url, data=data)
    r = requests.get(url, params=data)

    if verbose:
        print(r.url)
    r.raise_for_status()
    if format == "json":
        return r.json()
    else:
        return r.text


def checklegal(table,release):
    """Checks if this combination of table and release is acceptable
    
    Raises a VelueError exception if there is problem
    """
    
    releaselist = ("dr1", "dr2")
    if release not in ("dr1","dr2"):
        raise ValueError("Bad value for release (must be one of {})".format(', '.join(releaselist)))
    if release=="dr1":
        tablelist = ("mean", "stack")
    else:
        tablelist = ("mean", "stack", "detection")
    if table not in tablelist:
        raise ValueError("Bad value for table (for {} must be one of {})".format(release, ", ".join(tablelist)))


def ps1metadata(table="mean",release="dr2",
           baseurl="https://catalogs.mast.stsci.edu/api/v0.1/panstarrs"):
    """Return metadata for the specified catalog and table
    
    Parameters
    ----------
    table (string): mean, stack, or detection
    release (string): dr1 or dr2
    baseurl: base URL for the request
    
    Returns an astropy table with columns name, type, description
    """
    
    checklegal(table,release)
    url = "{baseurl}/{release}/{table}/metadata".format(**locals())
    r = requests.get(url)
    r.raise_for_status()
    v = r.json()
    # convert to astropy table
    tab = Table(rows=[(x['name'],x['type'],x['description']) for x in v],
               names=('name','type','description'))
    return tab


def parse_objinfoflag(objinfoflag):
    """ Get object flags.
    See reference info at https://outerspace.stsci.edu/display/PANSTARRS/PS1+Object+Flags and
    https://outerspace.stsci.edu/display/PANSTARRS/PS1+Database+object+and+detection+tables+for+bad+skycells
    """

    k = objinfoflag
    factors = [] 
    while True: 
        for j in 2**np.linspace(30, 1, 30, dtype=int): 
            if j <= k: 
                print(i, j, k) 
                factors.append(j) 
                k -= j 
                break 
        if k == 0: 
            break 

        return factors
