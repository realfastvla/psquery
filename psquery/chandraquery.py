# TODO: potential for image access via cadc API
#
#try:
#    from astroquery import cadc
#except ImportError:
#    print("astroquery not available. Cannot use irsaquery.")

import pyvo
from . import get_coord


def cone_chandra(radec, radius=5):
    """Query Chandra Source Catalog 2 with pyvo.
    radius in arcseconds.
    """

    url = "http://cda.cfa.harvard.edu/csc2scs/coneSearch"
    co = get_coord(radec, ret="skycoord")
    tab = pyvo.conesearch(url, co, radius=radius/3600).to_table()

    return tab
