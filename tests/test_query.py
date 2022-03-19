import pytest
from astropy import coordinates
from psquery import get_coord, astronquery, casdaquery, heasarcquery, irsaquery, mastquery, noaoquery, psquery, vizierquery

radecn = '13:31:08.288,30:30:32.9'
#'19:58:17.74, +34:34:38.7'
radecn1 = '13:29:52.7, 47:11:43'
radecs = '00:22:25.425930, 00:14:56.161440'
          #(221.2578, -13.9480)
radius = 5/3600


def test_coord1():
    ra, dec = get_coord(radecn, ret="radec")
    assert isinstance(ra, float)
    assert isinstance(dec, float)

def test_coord2():
    co = get_coord(radecs, ret="skycoord")
    assert isinstance(co, coordinates.SkyCoord)

def test_lotss():
    tab = astronquery.cone_lotss(radecn, radius)
    if tab is not None:
        print(len(tab))

def test_tgss():
    tab = astronquery.cone_tgss(radecs, radius*10)
    if tab is not None:
        print(len(tab))

def test_racs():
    tab = casdaquery.cone_racs(radecs, radius)
    if tab is not None:
        print(len(tab))
    
def test_vlssr():
    tab = heasarcquery.cone_vlssr(radecn1, radius)
    if tab is not None:
        print(len(tab))
    
def test_first():
    tab = heasarcquery.cone_first(radecn, radius)
    if tab is not None:
        print(len(tab))
    
def test_nvss():
    tab = heasarcquery.cone_nvss(radecn, radius)
    if tab is not None:
        print(len(tab))
    
def test_wenss():
    tab = heasarcquery.cone_wenss(radecn, radius)
    if tab is not None:
        print(len(tab))
    
def test_wise():
    pass
#psquery/irsaquery.py:def cone_wise(
#psquery/irsaquery.py:def cone_pyvo(ra, dec, radius=5, table="allwise_p3as_psd"):

def test_ps1strm():
    tab = mastquery.cone_ps1strm(radecn, radius)
    if tab is not None:
        print(len(tab))

# others
#psquery/mastquery.py:def match_ps1strm(radec, radius, verbose=False):
#psquery/mastquery.py:def cone_ps1_psc(radec, radius=5):
#psquery/mastquery.py:def cone_emline(
#psquery/mastquery.py:def cone_galaxymass(ra, dec, radius=5, selectcol=[]):
#psquery/mastquery.py:def cone_ps1_casjobs(radec, radius=5, ndet=1, nr=1, query="MeanObject"):

def test_legacy():
    tab = noaoquery.cone_legacy(radecs, radius)
    if tab is not None:
        print(len(tab))

def test_ps1():
    tab = psquery.query_radec(radecn, radius)
    if tab is not None:
        print(len(tab))

def test_sed():
    pass
#psquery/sed.py:def extinct(ra, dec, phot):
#psquery/sed.py:def get_phot(ra, dec, radius, **kwargs):
#psquery/sed.py:def run_fit(phot, hfile="results.h5", emcee=False, plot=True, **params):

def test_twomass():
    tab = irsaquery.cone_twomass(radecn, radius, catname="fp_psc")
    if tab is not None:
        print(len(tab))
    
def test_gleam():
    tab = vizierquery.cone_gleam(radecs, radius)
    if tab is not None:
        print(len(tab))

