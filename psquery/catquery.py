import numpy as np

from astropy import table
from astropy import coordinates
from astropy import units

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

    if "CIRADA" in filename:
        colnames = ['Component_name', 'RA', 'DEC', 'E_RA', 'E_DEC', 'Total_flux', 'E_Total_flux', 'Peak_flux', 'E_Peak_flux', 'Maj', 'E_Maj', 'Min', 'E_Min', 'PA', 'E_PA', 'Tile', 'NVSS_distance', 'Peak_to_ring', 'Duplicate_flag', 'Quality_flag']
    else:
        colnames = None

    tab = table.Table.read(filename, guess=False, format="csv", include_names=colnames)

    if "CIRADA" in filename:
        tab.rename_column('Component_name', 'vlassname')
        tab.rename_column("RA", "ra")
        tab.rename_column("DEC", "dec")
        sel = (tab['Duplicate_flag'] < 2) & ((tab['Quality_flag'] == 0) | (tab['Quality_flag'] == 4))
        return tab[sel]
    else:
        return tab


def cone_vlass(radec, radius=2, tab=None, filename='/Users/claw/data/vlass/CIRADA_VLASS1QL_table1_components_v2.csv'):
    """ Run cone search for radec tuple on vlass catalogs.
    Assumes vlass columns in astropy table format.
    radec tuple for coordinate.
    radius in arcseconds.
    """

    if tab is None:
        if filename is not None:
            tab = parse_cat(filename=filename)
        else:
            print("Must provide tab or filename")

    assert isinstance(tab, table.Table)

    co = get_coord(radec, ret="skycoord")
    cat = coordinates.SkyCoord(tab['ra'], tab['dec'], unit=units.deg)

    seps = cat.separation(co)
    sel = np.where(seps < radius*units.arcsec)[0]

    return tab[sel]

