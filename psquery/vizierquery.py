from astropy import coordinates, units

try:
    from astroquery import vizier
except ImportError:
    print('astroquery not available. Cannot use mwaquery.')
from . import get_coord

vz = vizier.Vizier()

def cone_gleam(radec, radius):
    """ query gleam catalog
    radec can be any format (parsed by get_coord).
    radius in degrees.
    """

    ra, dec = get_coord(radec)
    co = coordinates.SkyCoord(ra, dec, unit=units.deg)
    res = vz.query_region(co, radius=radius*units.deg, catalog='VIII/100/gleamegc')

    return res.values()
