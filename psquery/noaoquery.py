from dl import queryClient as qc
from dl.helpers.utils import convert


def query_legacy(ra, dec, radius = 5):
    """ Query legacy catalog with (ra, dec) in degrees.
    Radius in degrees.
    Returns dataframe.
    """

    query1 = f"""
    SELECT ra,dec,type,ls_id,
    22.5-2.5*log10(flux_g) as mag_g,22.5-2.5*log10(flux_r) as mag_r,22.5-2.5*log10(flux_z) as mag_z,
    2.5*(flux_ivar_g)^(-1/2)/(flux_g*ln(10)) as mag_g_err,2.5*(flux_ivar_r)^(-1/2)/(flux_r*ln(10)) as mag_r_err,2.5*(flux_ivar_z)^(-1/2)/(flux_z*ln(10)) as mag_z_err
    FROM ls_dr9.tractor
    WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {radius/3600})
    """

    result1 = qc.query(sql=query1)
    df1 = convert(result1)
    return df1
