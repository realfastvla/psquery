from astropy.io import ascii
from astropy import coordinates
from astropy import units as u

from . import get_coord


def parse_cat(filename):
    """ Parse a catalog from ascii format file and return astropy table.
    Enforces standard name for "ra" and "dec" in table, to simplify cone search.

    Some files floating around:
    - CIRADA_VLASS1QL_table1_components_v2.csv -- CIRADA public catalog
    - all_7sigma_e12_sources.csv
    - clean_e2_and_e1_table.csv
    - all_e2_gtr_1mJy.csv
    """

    pass


def cone_vlass(radec, radius=2, tab=None, filename=''):
    """ Run cone search for radec tuple on vlass catalogs.
    Assumes vlass columns in astropy table format.
    radec tuple for coordinate.
    radius in arcseconds.
    """

    if tab is None and filename is not None:
        tab = parse_cat(filename=filename)

    co = get_coord(radec, ret="skycoord")
    cat = coordinates.SkyCoord(tab['ra'], tab['dec'], unit=u.deg)

    sep = cat.separation(co).to_value(units.arcsec)
    sel = np.where(seps < radius*units.arcsec)[0]

    return cat[sel]

