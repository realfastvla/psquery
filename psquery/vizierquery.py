from astropy import coordinates, units
from astroquery.utils.tap.core import TapPlus

try:
    from astroquery import vizier
except ImportError:
    print("astroquery not available. Cannot use mwaquery.")
from . import get_coord

vz = vizier.Vizier()


def cone_gleam(radec, radius):
    """query gleam catalog
    radec can be any format (parsed by get_coord).
    radius in degrees.
    """

    ra, dec = get_coord(radec)
    co = coordinates.SkyCoord(ra, dec, unit=units.deg)
    res = vz.query_region(co, radius=radius * units.deg, catalog="VIII/100/gleamegc")

    return res.values()


def cone_xmm(radec, radius, selectcol=None, table='IX/65/xmm4d11s'):
    """ Use TAP interface to XMM LSS.
    radec can be any format (parsed by get_coord).
    radius in degrees.
    """

    vz = TapPlus(url='http://tapvizier.u-strasbg.fr/TAPVizieR/tap')
    ra, dec = get_coord(radec)
    query = f"""SELECT TOP 100 "{table}".Source,  "{table}"."4XMM",  "{table}".RA_ICRS,  "{table}".DE_ICRS,  "{table}".srcML,  "{table}".Flux8,  "{table}".e_Flux8, "{table}".HR1,  "{table}".HR2,  "{table}".HR3,  "{table}".HR4,  "{table}".ext,  "{table}".V,  "{table}".S,  "{table}".Nd
 FROM "{table}"
 WHERE 1=CONTAINS(POINT('ICRS',"{table}".RA_ICRS,"{table}".DE_ICRS), CIRCLE('ICRS', {ra}, {dec}, {radius}))"""

    tab = vz.launch_job_async(query).get_results()

    if len(tab):
        if selectcol:
            return tab[selectcol]
        else:
            return tab
    else:
        return None
