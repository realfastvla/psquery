from astropy import coordinates, units
from astroquery import heasarc

from . import get_coord

heq = heasarc.Heasarc()


def cone_vlssr(radec, radius):
    """
    radec is parsed by get_coord
    radius in degrees.
    returns dataframe
    """

    coord = get_coord(radec, ret="skycoord")
    try:
        tab = heq.query_region(coord, mission="vlssr", radius=radius*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()


def cone_first(radec, radius):
    """
    radec is parsed by get_coord
    radius in degrees
    returns dataframe
    """

    coord = get_coord(radec, ret="skycoord")

    try:
        tab = heq.query_region(coord, mission="first", radius=radius*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()


def cone_nvss(radec, radius):
    """
    radec is parsed by get_coord
    radius in degrees
    returns dataframe
    """

    coord = get_coord(radec, ret="skycoord")

    try:
        tab = heq.query_region(coord, mission="nvss", radius=radius*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()


def cone_wenss(radec, radius):
    """
    radec is parsed by get_coord
    radius in degrees
    returns dataframe
    """

    coord = get_coord(radec, ret="skycoord")

    try:
        tab = heq.query_region(coord, mission="wenss", radius=radius*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()
