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

from 
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
