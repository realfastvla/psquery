__all__ = ['psquery', 'irsaquery', 'mastquery', 'clutools', 'mwaquery', 'casdaquery', 'twomassquery', 'get_radec']

from astropy import coordinates, units

def get_radec(*args):
    if len(args) == 2:
        ra0, dec0 = args
        if isinstance(ra0, str):
            if ':' in ra0:
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=(units.hourangle, units.deg))
            else:
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)
        elif isinstance(ra0, float) or isinstance(ra0, int):
            co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)

    elif len(args) == 1:
        if not isinstance(args[0], coordinates.SkyCoord):
            try:
                ra0, dec0 = args[0].split(',')
            except ValueError:
                ra0, dec0 = args[0].split(' ')
            
            if ':' in ra0:
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=(units.hourangle, units.deg))
            else:
                co = coordinates.SkyCoord(ra=ra0, dec=dec0, unit=units.deg)
        else:
            co = args[0]
    else:
        print('could not parse args: ', args)

    return co.ra.value, co.dec.value
