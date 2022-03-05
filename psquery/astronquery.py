from astropy import table
import pyvo
from . import get_radec


def cone_lotss(radec, radius=5/3600, selectcol=['ra', 'dec', 'peak_flux', 'e_peak_flux', 'total_flux'], getepoch=True):
    """ cone search of LoTSS.
    ra, dec in any format (parsed by get_radec).
    radius in degrees.
    selectcol sets columns to return. None or empty list returns all columns.
    catalog can be 'initial' or 'hale'.
    """

    ra, dec = get_radec(radec)
    tap = pyvo.dal.TAPService('https://vo.astron.nl/__system__/tap/run/tap')
    query = f"SELECT * FROM lotss_dr2.main_sources where 1=CONTAINS(POINT('ICRS', RA, DEC), CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    tab = tap.search(query).to_table()[selectcol + ['mosaic_id']]
    if not len(tab):
        return

    dos = []
    for mosaic_id in tab['mosaic_id']:
        query = f"SELECT dateobs FROM lotss_dr2.mosaics where mosaic_id='{mosaic_id}'"
        do = tap.search(query).to_table()['dateobs'][0]
        dos.append(do)
    tab2 = table.Table(rows=[dos], names=['dateobs'])

    return table.hstack([tab, tab2])
