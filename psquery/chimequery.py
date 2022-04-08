import numpy as np
from astropy import coordinates, table, units
import pandas as pd

try:
    from cfod import catalog
    data = catalog.as_dataframe()
except ImportError:
    try:
        data = pd.read_csv("/home/ubuntu/chimefrbcat1.csv")
        print("read cfod catalog from disk")
    except:
        print("cfod catalog not available. Cannot use chimequery.")

from . import get_coord


def cone_catalog1(radec, radius=1, selectcols=["tns_name", "repeater_name", "ra", "dec", "bonsai_dm"]):
    """ Do cone search on CHIME/FRB catalog1.
    Returns DataFrame with selectcol columns plus a "separation" column (in degrees).
    """

    cat = coordinates.SkyCoord(data[["ra", "dec"]].to_numpy(), unit='deg')
    co = get_coord(radec, ret="skycoord")

    seps = cat.separation(co)
    sel = np.where(seps < radius*units.deg)[0]
    data_close = data.iloc[sel][selectcols]
    df_sep = pd.DataFrame(seps[np.where(seps < radius*units.deg)[0]], columns=['separation'], index=data_close.index)
    
    return pd.concat([data_close, df_sep], join="outer", axis=1)
