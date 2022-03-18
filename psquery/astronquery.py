import numpy as np
from astropy import table, coordinates, units
try:
    import pyvo
except ImportError:
    print('pyvo not available. Cannot use astronquery.')
from . import get_coord

def cone_lotss(radec, radius=5/3600, selectcol=['ra', 'dec', 'peak_flux', 'e_peak_flux', 'total_flux'], getepoch=True):
    """ cone search of LoTSS.
    ra, dec (in any format parsed by get_radec).
    radius in degrees.
    selectcol sets columns to return. None or empty list returns all columns.
    catalog can be 'initial' or 'hale'.
    """

    ra, dec = get_coord(radec, ret='radec')
    tap = pyvo.dal.TAPService('https://vo.astron.nl/__system__/tap/run/tap')
    query = f"SELECT * FROM lotss_dr2.main_sources where 1=CONTAINS(POINT('ICRS', RA, DEC), CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    tab = tap.search(query).to_table()[selectcol + ['mosaic_id']]
    if not len(tab):
        return

#    co = coordinates.SkyCoord(float(ra), float(dec), unit=(units.deg, units.deg))
#    col = coordinates.SkyCoord(float(tab['ra']), float(tab['dec']), unit=(units.deg, units.deg))
#    sep = co.separation(col).to_value(units.arcsec)

    dos = []
    for mosaic_id in tab['mosaic_id']:
        query = f"SELECT dateobs FROM lotss_dr2.mosaics where mosaic_id='{mosaic_id}'"
        do = tap.search(query).to_table()['dateobs'][0]
        dos.append(do)

    tab2 = table.Table(data=[dos], names=['dateobs'])

    return table.hstack([tab, tab2])

def xmatch_lotss(radecs, radius=5/3600):
    """ Given list of (RA, Dec) tuples, list all LoTSS matches.
    Prints [RA, Dec, Fp, Fint, sep, epoch]
    radius in degrees.
    """

    for ra, dec in radecs:
        tab = cone_lotss((ra, dec), radius=radius)
        if tab is not None:
            if len(tab) == 1:
                ral, decl, fp, fi, do = tab['ra', 'dec', 'peak_flux', 'total_flux', 'dateobs'][0]
                co = coordinates.SkyCoord(float(ra), float(dec), unit=(units.deg, units.deg))
                col = coordinates.SkyCoord(float(ral), float(decl), unit=(units.deg, units.deg))
                sep = co.separation(col).to_value(units.arcsec)
                print(ral, decl, fp, fi, sep, do)
            else:
                print(f'Found {len(tab)} LoTSS counterparts. Skipping...')
                continue
        else:
            print()
