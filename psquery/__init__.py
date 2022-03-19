__all__ = [
    "astronquery",
    "casdaquery",
    "clutools",
    "heasarcquery",
    "irsaquery",
    "mastquery",
    "noaoquery",
    "psquery",
    "sed",
    "vizierquery",
    "get_coord",
]

from astropy import coordinates, units

def get_coord(*args, ret="radec"):
    if len(args) == 2:
        ra0, dec0 = args

    elif len(args) == 1:
        if isinstance(args[0], tuple):
            ra0, dec0 = args[0]
            if isinstance(ra0, str):
                if ":" in ra0:
                    co = coordinates.SkyCoord(
                        ra=ra0, dec=dec0, unit=(units.hourangle, units.deg)
                    )
                else:
                    co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)
            elif isinstance(ra0, float) or isinstance(ra0, int):
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)
        elif not isinstance(args[0], coordinates.SkyCoord):
            try:
                ra0, dec0 = args[0].split(",")
            except ValueError:
                ra0, dec0 = args[0].split(" ")

            if ":" in ra0:
                co = coordinates.SkyCoord(
                    ra=ra0, dec=dec0, unit=(units.hourangle, units.deg)
                )
            else:
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)
        else:
            co = args[0]
    else:
        print("could not parse args: ", args)

    if ret == "radec":
        return co.ra.value, co.dec.value
    elif ret == "skycoord":
        return co
    else:
        print("ret must be radec or skycoord")


def radio_survey_plot():
    """ Compile and visualize radio surveys in time, freq.
    """

    # survey, freq, area, start, stop
    # first, 1.4, 10000, 1994, 2005
    # nvss, 1.4, 33885, 1993, 1998
    # wenss, 0.325, ?, 1997
    # vlass, 3.0, 33885, 2017, 2022
    # vcss, 0.34, 33885, 2017?, 2022
    # tgss, 0.15, 37100, 2010, 2012
    # racs, 0.887, 34240, 2019, 2022?
    # lotss, 0.15, 5600, , 2022?
    # sumss (max dec=-30)
    pass

from psquery.astronquery import cone_lotss, cone_tgss
from psquery.casdaquery import cone_racs
from psquery.heasarcquery import cone_first, cone_nvss, cone_wenss, cone_vlssr
from psquery.irsaquery import cone_wise, cone_twomass # cone_pyvo should be renamed
from psquery.mastquery import cone_ps1strm, cone_ps1psc, cone_emline, cone_galaxymass
from psquery.noaoquery import cone_legacy
#from psquery.psquery import cone_ps1 # need to rename still
from psquery.vizierquery import cone_gleam

