from astroquery import heasarc
from astropy import coordinates, units

heq = heasarc.Heasarc()

def query_vlssr(source, radius):
    """ Source can be an radec string (no comma, hmsdms or deg) or source name as string.
    """

    try:
        tab = heq.query_region(source, mission='vlssr', radius=0.5*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()


def query_first(source, radius):
    """ Source can be an radec string (no comma, hmsdms or deg) or source name as string.
    """

    try:
        tab = heq.query_region(source, mission='first', radius=0.5*units.deg)
    except TypeError:
        print("No sources found")
        tab = None

    return tab.to_pandas()
