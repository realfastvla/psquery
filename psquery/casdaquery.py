from astroquery.utils.tap.core import TapPlus
from . import get_radec


def cone_racs(radec, radius=5/3600, selectcol=['id', 'ra_deg_cont', 'dec_deg_cont', 'ra_err',
                                             'dec_err', 'flux_peak', 'flux_peak_err', 'flux_int',
                                             'flux_int_err', 'maj_axis', 'min_axis']):
    """ cone search of ASKAP/RACS.
    ra, dec in any format (parsed by get_radec).
    radius in degrees.
    """

    ra, dec = get_radec(radec)
    
#    base = 'http://casda.csiro.au/casda_vo_tools/tap/sync?request=doQuery&lang=ADQL&format=csv&query='
    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")

    query = f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont), CIRCLE('ICRS', {ra}, {dec}, {radius})) and project_id = 23" 

#    tab = ascii.read(data.decode('utf-8'))
    tab = casdatap.launch_job(query).get_results()
    
    return tab[selectcol]
