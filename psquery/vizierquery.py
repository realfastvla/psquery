from astropy import coordinates, units
from astroquery.utils.tap.core import TapPlus
import pyvo as vo

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

def cone_vlass(radec, radius):
    """
    """

    co = get_coord(radec, ret='skycoord')
    size = radius*units.deg
    scs_srv = vo.dal.SCSService("http://vizier.cds.unistra.fr/viz-bin/conesearch/J/ApJS/255/30/comp?")
    vlass_table = scs_srv.search(pos=co, radius=size).to_table()

    tab = vlass_table[['CompName', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000', 'Ftot', 'e_Ftot', 'Fpeak', 'e_Fpeak', 'DupFlag', 'QualFlag']]
    tab.rename_column('CompName', 'vlassname')
    tab.rename_column("RAJ2000", "ra")
    tab.rename_column("DEJ2000", "dec")
    sel = (tab["DupFlag"] < 2) & ((tab["QualFlag"] == 0) | (tab["QualFlag"] == 4))

    return tab[sel]
