from astroquery.utils.tap.core import TapPlus

def cone_racs(ra, dec, radius=5., selectcol=['id', 'ra_deg_cont', 'dec_deg_cont', 'ra_err',
                                             'dec_err', 'flux_peak', 'flux_peak_err', 'flux_int',
                                             'flux_int_err', 'maj_axis', 'min_axis']):
    """ cone search of ASKAP/RACS.
    ra, dec in degrees, radius in arcsec.
    """

#    base = 'http://casda.csiro.au/casda_vo_tools/tap/sync?request=doQuery&lang=ADQL&format=csv&query='
    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")

    query = f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont),CIRCLE('ICRS',{ra},{dec},{radius})) and project_id = 23" 

#    tab = ascii.read(data.decode('utf-8'))
    tab = casdatap.launch_job(query).get_results()
    
    return tab[selectcol]
