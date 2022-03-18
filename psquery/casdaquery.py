try:
    from astroquery.utils.tap.core import TapPlus
except ImportError:
    print('astroquery not available. Cannot use casdaquery.')
from . import get_radec


def cone_racs(radec, radius=5/3600, selectcol=['id', 'ra_deg_cont', 'dec_deg_cont', 'ra_err',
                                               'dec_err', 'flux_peak', 'flux_peak_err', 'flux_int',
                                               'flux_int_err', 'maj_axis', 'min_axis'],
              catalog='initial'):
    """ cone search of ASKAP/RACS.
    ra, dec in any format (parsed by get_radec).
    radius in degrees.
    selectcol sets columns to return. None or empty list returns all columns.
    catalog can be 'initial' or 'hale'.
    """

    ra, dec = get_radec(radec)
    
#    base = 'http://casda.csiro.au/casda_vo_tools/tap/sync?request=doQuery&lang=ADQL&format=csv&query='
    casdatap = TapPlus(url="https://casda.csiro.au/casda_vo_tools/tap")

    if catalog == 'initial':
        query = f"SELECT * FROM casda.continuum_component where 1=CONTAINS(POINT('ICRS', ra_deg_cont, dec_deg_cont), CIRCLE('ICRS', {ra}, {dec}, {radius})) and project_id = 23"
    elif catalog == 'hale':
        query = f"SELECT * FROM AS110.racs_dr1_gaussians_galacticcut_v2021_08_v01 where 1=CONTAINS(POINT('ICRS', ra, dec),CIRCLE('ICRS',{ra},{dec},{radius}))"
    else:
        print('catalog not recognized')
        return
        
#    tab = ascii.read(data.decode('utf-8'))
    tab = casdatap.launch_job_async(query).get_results()

    if selectcol:
        return tab[selectcol]
    else:
        return tab
