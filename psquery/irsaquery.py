from astroquery import irsa
from astropy import coordinates, units
from . import get_radec

def cone_wise(radec, radius=5/3600, selectcol=['designation', 'ra', 'dec', 'w1mag', 'w1flg', 'w2mag', 'w2flg',
                                                 'w3mag', 'w3flg', 'w4mag', 'w4flg']):
    """ cone search from wise cryogenic all-sky survey
    ra, dec in degrees, radius in arcsec.
    Column definitions at https://wise2.ipac.caltech.edu/docs/release/allsky/expsup/sec2_2a.html
    """

    ra, dec = get_radec(radec)
    co = coordinates.SkyCoord(ra, dec, unit='deg')
    tab = irsa.Irsa.query_region(co, catalog='allsky_4band_p3as_psd')[selectcol]

    return tab


def cone_pyvo(ra, dec, radius=5, table='allwise_p3as_psd'):
    """ Query IRSA with pyvo
    TODO: use unwise? need to figure out which columns

    Default table is allwise.
    Cone search radius in arcseconds.
    """

    url = "https://irsa.ipac.caltech.edu/SCS?table=allwise_p3as_psd"
    co = coordinates.SkyCoord(ra, dec, unit='deg')
    tab = pyvo.conesearch(url, pos=(co.ra.deg, co.dec.deg), radius = rad/3600).to_table()

    return tab


def match_wise(ra, dec):
    """
    for j, (ra, dec, name) in enumerate(tab[['RA (deg)', 'DEC (deg)', 'Short name']]): 
     ...:     if j < j0: 
     ...:         continue 
     ...:     co = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg)) 
     ...:     res = irsaquery.cone_wise(ra, dec) 
     ...:     if res: 
     ...:         rao, deco = res['ra'], res['dec'] 
     ...:         coo = coordinates.SkyCoord(rao, deco, unit=(units.deg, units.deg)) 
     ...:         sep = co.separation(coo).to_value(units.arcsec) 
     ...:         i = np.where(sep == sep.min())[0][0] 
     ...: #        print(f'using {i}/{len(sep)}') 
     ...:         designation, w1mag, w2mag, w3mag, w4mag = res['designation'].tolist()[i].decode('utf-8'), res['w1mag'].tolist()[i], res['w2
     ...: mag'].tolist()[i], res['w3mag'].tolist()[i], res['w4mag'].tolist()[i] 
     ...:         print(f'{designation}\t{w1mag}\t{w2mag}\t{w3mag}\t{w4mag}') 
     ...:     else: 
     ...:         print() 
     ...:          
    """


    pass
