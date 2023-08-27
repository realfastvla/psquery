__all__ = [
    "astronquery",
    "casdaquery",
    "catquery",
    "clutools",
    "heasarcquery",
    "irsaquery",
    "psquery",
    "sed",
    "vizierquery",
    "get_coord",
    "radio_survey_data",
# these are noisy
    "mastquery",
    "noaoquery",
    "chimequery"
]

from astropy import coordinates, units, io

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


def radio_survey_data():
    """ Compile and visualize radio surveys in time, freq.
    """

    return io.ascii.read('''
    name, freq, area, start, stop, notes
    sumss , , , , , south of -30 so not much overlap with others
    gb6, 5.0, 20000, 1986, 1987, 0<dec<75
    wenss, 0.325, 10300, 1991, 1997, north of +30
    nvss, 1.4, 33885, 1993, 1998, north of -40
    first, 1.4, 10000, 1994, 2005, north/south galactic caps
    tgss, 0.15, 37100, 2010, 2012, north of -53
    gleam, 0.15, 33885, 2013, 2015, south of +30
    vlass, 3.0, 33885, 2017, 2022, north of -40 (ongoing)
    vcss, 0.34, 33885, 2017, 2022, north of -40 (ongoing. start may be wrong)
    racs, 0.887, 34240, 2019, 2022, south of +30 (ongoing)
    lotss, 0.15, 5600, 2014, 2020, subset of north of 0
    ''')

from psquery.astronquery import cone_lotss, cone_tgss
from psquery.casdaquery import cone_racs
from psquery.heasarcquery import cone_first, cone_nvss, cone_wenss, cone_vlssr
from psquery.irsaquery import cone_wise, cone_twomass
from psquery.mastquery import cone_ps1strm, cone_ps1psc, cone_emline, cone_galaxymass, cone_galex
from psquery.noaoquery import cone_legacy
from psquery.psquery import cone_ps1
from psquery.vizierquery import cone_gleam, cone_xmm, cone_vlass
from psquery.chandraquery import cone_chandra
# TODO: add hsc queries with pyvo(?)
