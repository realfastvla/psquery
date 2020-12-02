from astropy.io import ascii
from urllib.request import urlopen

def cone_racs(ra, dec, radius=5., selectcol=['id', 'ra_deg_cont', 'dec_deg_cont', 'ra_err',
                                             'dec_err', 'flux_peak', 'flux_peak_err', 'flux_int',
                                             'flux_int_err', 'maj_axis', 'min_axis']):
    """ cone search of ASKAP/RACS.
    ra, dec in degrees, radius in arcsec.
    """

    base = 'http://casda.csiro.au/casda_vo_tools/tap/sync?request=doQuery&lang=ADQL&format=csv&query='
    adql = f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont),CIRCLE('ICRS',{ra},{dec},{radius})) and project_id = 23" 
    url = base+adql
    url = url.replace(' ','%20') 
    url = url.replace(chr(34),'%27') 
    data = urlopen(url).read()
    tab = ascii.read(data.decode('utf-8'))

    return tab[selectcol]
