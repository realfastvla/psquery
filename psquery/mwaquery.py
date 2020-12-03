from astropy import coordinates, units
from astroquery import vizier

vz = vizier.Vizier()

def cone_gleam(ra, dec, radius):
    """ query gleam catalog
    """

    co = coordinates.SkyCoord(ra, dec, unit=units.deg)
    res = vz.query_region(co, radius=1*units.deg, catalog='VIII/100/gleamegc')

    return res.values()
