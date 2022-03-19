import numpy as np
from astropy import coordinates, table, units

try:
    import pyvo
except ImportError:
    print("pyvo not available. Cannot use astronquery.")
from . import get_coord


def cone_lotss(
    radec,
    radius=5 / 3600,
    selectcol=["ra", "dec", "peak_flux", "e_peak_flux", "total_flux"],
    getepoch=True,
):
    """cone search of LoTSS.
    ra, dec (in any format parsed by get_radec).
    radius in degrees.
    selectcol sets columns to return. None or empty list returns all columns.
    catalog can be 'initial' or 'hale'.
    """

    ra, dec = get_coord(radec, ret="radec")
    tap = pyvo.dal.TAPService("https://vo.astron.nl/__system__/tap/run/tap")
    query = f"SELECT * FROM lotss_dr2.main_sources where 1=CONTAINS(POINT('ICRS', RA, DEC), CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    tab = tap.search(query).to_table()[selectcol + ["mosaic_id"]]
    if not len(tab):
        if False:   # need to get function to find coverage
            return -0.8 # 90% completeness in mJy
        else:
            return

    #    co = coordinates.SkyCoord(float(ra), float(dec), unit=(units.deg, units.deg))
    #    col = coordinates.SkyCoord(float(tab['ra']), float(tab['dec']), unit=(units.deg, units.deg))
    #    sep = co.separation(col).to_value(units.arcsec)

    dos = []
    for mosaic_id in tab["mosaic_id"]:
        query = f"SELECT dateobs FROM lotss_dr2.mosaics where mosaic_id='{mosaic_id}'"
        do = tap.search(query).to_table()["dateobs"][0]
        dos.append(do)

    tab2 = table.Table(data=[dos], names=["dateobs"])

    return table.hstack([tab, tab2])


def xmatch_lotss(radecs, radius=5 / 3600):
    """Given list of (RA, Dec) tuples, list all LoTSS matches.
    Prints [RA, Dec, Fp, Fint, sep, epoch]
    radius in degrees.
    """

    for ra, dec in radecs:
        tab = cone_lotss((ra, dec), radius=radius)
        if tab is not None:
            if len(tab) == 1:
                ral, decl, fp, fi, do = tab[
                    "ra", "dec", "peak_flux", "total_flux", "dateobs"
                ][0]
                co = coordinates.SkyCoord(
                    float(ra), float(dec), unit=(units.deg, units.deg)
                )
                col = coordinates.SkyCoord(
                    float(ral), float(decl), unit=(units.deg, units.deg)
                )
                sep = co.separation(col).to_value(units.arcsec)
                print(ral, decl, fp, fi, sep, do)
            else:
                print(f"Found {len(tab)} LoTSS counterparts. Skipping...")
                continue
        else:
            print()


def cone_tgss(
    radec,
    radius=5/3600,
    selectcol=["ra", "dec", "spk", "e_spk", "sint", "e_sint"],
):
    """cone search of TGSS
    ra, dec (in any format parsed by get_radec).
    radius in degrees.
    selectcol sets columns to return. None or empty list returns all columns.
    """

    ra, dec = get_coord(radec, ret="radec")
    tap = pyvo.dal.TAPService("https://vo.astron.nl/__system__/tap/run/tap")
    query = f"SELECT * FROM tgssadr.main where 1=CONTAINS(POINT('ICRS', RA, DEC), CIRCLE('ICRS', {ra}, {dec}, {radius}))"
    tab = tap.search(query).to_table()[selectcol]

    if not len(tab):
        if dec > -53:
            return -5*3.5 # 5 sigma in mJy
        else:
            return

    #    co = coordinates.SkyCoord(float(ra), float(dec), unit=(units.deg, units.deg))
    #    col = coordinates.SkyCoord(float(tab['ra']), float(tab['dec']), unit=(units.deg, units.deg))
    #    sep = co.separation(col).to_value(units.arcsec)

#    dos = []
#    for ra, dec in tab["ra", "dec"]:
#        query = f"SELECT dateObs FROM tgssadr.img_main where 1=CONTAINS(POINT('ICRS', centerAlpha, centerDelta), CIRCLE('ICRS', {ra}, {dec}, {radius}))"
#        do = tap.search(query).to_table()["dateObs"][0]
#        dos.append(do)
#
#    tab2 = table.Table(data=[dos], names=["dateObs"])
#    return table.hstack([tab, tab2])

    return tab

