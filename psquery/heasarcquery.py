from astropy import coordinates, units
from astroquery import heasarc

from . import get_coord

heq = heasarc.Heasarc()

def query_vlssr(radec, radius):
    """ 
    radec is parsed by get_coord
    returns dataframe
    """

    coord = get_coord(radec, ret='skycoord')
    try:
        tab = heq.query_region(coord, mission='vlssr', radius=0.5*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()


def query_first(radec, radius):
    """ 
    radec is parsed by get_coord
    returns dataframe
    """

    coord = get_coord(radec, ret='skycoord')

    try:
        tab = heq.query_region(source, mission='first', radius=0.5*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()
