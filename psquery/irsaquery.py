try:
    from astroquery import irsa
except ImportError:
    print("astroquery not available. Cannot use irsaquery.")

from astropy import coordinates, units
import pandas as pd
import requests

from . import get_coord


def cone_wise(
    radec,
    radius=5 / 3600,
    selectcol=[
        "designation",
        "ra",
        "dec",
        "w1mpro",
        "w1sigmpro",
        "w2mpro",
        "w2sigmpro",
        "w3mpro",
        "w3sigmpro",
        "w4mpro",
        "w4sigmpro",
    ],
    catalog="allwise_p3as_psd",
):
    """cone search from wise cryogenic all-sky survey
    ra, dec in degrees, radius in arcsec.
    https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-dd?catalog=allwise_p3as_psd
    """

    ra, dec = get_coord(radec, ret="radec")
    co = coordinates.SkyCoord(ra, dec, unit="deg")
    tab = irsa.Irsa.query_region(co, catalog=catalog, selcols=",".join(selectcol))

    return tab


def cone_wise_pyvo(ra, dec, radius=5, table="allwise_p3as_psd"):
    """Query IRSA with pyvo
    TODO: use unwise? need to figure out which columns

    Default table is allwise.
    Cone search radius in arcseconds.
    """

    url = "https://irsa.ipac.caltech.edu/SCS?table=allwise_p3as_psd"
    co = coordinates.SkyCoord(ra, dec, unit="deg")
    tab = pyvo.conesearch(
        url, pos=(co.ra.deg, co.dec.deg), radius=rad / 3600
    ).to_table()

    return tab


def cone_twomass(radec, radius=5 / 3600, catname="fp_psc"):
    """Do a cone search of the 2MASS catalog

    Parameters
    ----------
    radec can be any format (parsed by get_coord)
    radius (float): (degree) Search radius
    catname: fp_psc (for 2MASS)
    """

    radius = str(radius * 3600)
    cols = "ra,dec,designation,ndet"
    outformat = "1"
    ra, dec = get_coord(radec)
    co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))
    coord = f"{co.ra.deg}+{co.dec.deg}"

    baseurl = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?spatial=cone"
    url = baseurl + f"&catalog={catname}"
    url += f"&objstr={coord}"
    url += f"&radius={radius}&radunits=arcsec"
    url += f"&outfmt={outformat}"
    url += f"&selcols={cols}"

    print(url)
    r = requests.get(url)
    r.raise_for_status()
    rows = r.text.split("\n")

    for i, r in enumerate(rows):
        if r[0] == "|":
            break
    l = []
    for r in rows[i + 4 :]:
        if len(r):
            l.append(r.split())

    if len(l):
        col = []
        for c in rows[i].split():
            if c == "|":
                continue
            col.append(c[:-1])

        df = pd.DataFrame(l, columns=col)
        df["ra"] = [float(ra) for ra in df["ra"].tolist()]
        df["dec"] = [float(dec) for dec in df["dec"].tolist()]
        coords = coordinates.SkyCoord(
            df.ra.tolist() * units.deg, df.dec.tolist() * units.deg
        )
    else:
        print("Nothing found.")
        coords = None
    return coords
