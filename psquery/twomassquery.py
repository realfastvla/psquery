import pandas as pd
import requests
from astropy import coordinates, units

from . import get_coord


def query_radec(radec, radius=5/3600, catname='fp_psc'):
    """Do a cone search of the 2MASS catalog
    
    Parameters
    ----------
    radec can be any format (parsed by get_coord)
    radius (float): (degree) Search radius
    catname: fp_psc (for 2MASS)
    """
    
    radius = str(radius*3600)
    cols = 'ra,dec,designation,ndet'
    outformat = '1'
    ra, dec = get_coord(radec)
    co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))
    coord = f'{co.ra.deg}+{co.dec.deg}'
    
    baseurl = 'https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=cone'
    url = baseurl + f'&catalog={catname}'
    url += f'&objstr={coord}'
    url += f'&radius={radius}&radunits=arcsec'
    url += f'&outfmt={outformat}'
    url += f'&selcols={cols}'
    
    print(url)
    r = requests.get(url)
    r.raise_for_status()
    rows = r.text.split('\n')
    
    for i, r in enumerate(rows):
        if r[0] == '|':
            break
    l = []
    for r in rows[i+4:]:
        if len(r):
            l.append(r.split())
    
    if len(l):
        col = []
        for c in rows[i].split():
            if c == '|':
                continue
            col.append(c[:-1])

        df = pd.DataFrame(l, columns=col)
        df['ra'] = [float(ra) for ra in df['ra'].tolist()]
        df['dec'] = [float(dec) for dec in df['dec'].tolist()]
        coords = coordinates.SkyCoord(df.ra.tolist()*units.deg, df.dec.tolist()*units.deg)
    else:
        print('Nothing found.')
        coords = None
    return coords
