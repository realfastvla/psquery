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
        if coord.dec.value > -30:
            tab = -5*130 # 5 sigma limit in mJy
        else:
            tab = None

    return tab


def cone_first(radec, radius):
    """
    radec is parsed by get_coord
    radius in degrees
    returns dataframe
    Column for epoch is present.
    """

    coord = get_coord(radec, ret="skycoord")

    try:
        tab = heq.query_region(coord, mission="first", radius=radius*units.deg)
    except TypeError:
        print("No sources found")
        if coord.dec.value > -8 and coord.dec.value < 65 and coord.ra.value > 90 and coord.ra.value < 210:  # rough check
            tab = -5*0.15 # 5 sigma limit in mJy
        else:
            tab = None

    return tab


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
        if coord.dec.value > -40:
            tab = -5*0.45  # 5 sigma limit in mJy
        else:
            tab = None

    return tab


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
        if coord.dec.value > +30:
            tab = -18 # 5 sigma limit in mJy
        else:
            tab = None

    return tab
