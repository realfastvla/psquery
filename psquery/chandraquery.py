# TODO: potential for image access via cadc API
#
#try:
#    from astroquery import cadc
#except ImportError:
#    print("astroquery not available. Cannot use irsaquery.")

import pyvo
from . import get_coord


def cone_chandra(radec, radius=5, path="csc21_snapshot_scs"):
    """Query Chandra Source Catalog 2 with pyvo.
    radius in arcseconds.
    path is portion of VO cone search url (e.g., https://cda.cfa.harvard.edu/{path}/coneSearch)
    """

    oldpath = "csc2scs"
    url = f"https://cda.cfa.harvard.edu/{path}/coneSearch"
    co = get_coord(radec, ret="skycoord")
    tab = pyvo.conesearch(url, co, radius=radius/3600).to_table()

    return tab
