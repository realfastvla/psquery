import numpy as np
import pandas as pd
import h5py
from astropy import units as u
from astropy import coordinates

CLU_min_dist = 0*u.Mpc #change to 5 or something similar to exclude nearby resolved galaxies
CLU_max_dist = 200.0*u.Mpc


def compile_CLU_catalog(filename='CLU_20170106.hdf5', min_dist = CLU_min_dist, max_dist = CLU_max_dist):
    #Get CLU catalog
    f = h5py.File(filename, 'r')
    #f contains f['data'] and f['mask']. 
    #f['mask'] = True if the data has no value (i.e. 999999)
    #f['mask'] = False if the data has a value
    #f['data'] has a number of columns: 
    #They can be accessed as e.g. f['data']['ra']
    #COPY CLU (name,ra,dec,a,ratio_b2a,pa,z_NED,zerr_NED,dm_in,dm_range,source,ra_NED,dec_NED,type_NED,btc,B_r25,B_r25err,j_m_k20fe,j_msig_k20fe,h_m_k20fe,h_msig_k20fe,k_m_k20fe,k_msig_k20fe,FUV,FUVerr,NUV,NUVerr,ra_SDSS,dec_SDSS,petroMag_u,petroMag_g,petroMag_r,petroMag_i,petroMag_z,petroMagErr_u,petroMagErr_g,petroMagErr_r,petroMagErr_i,petroMagErr_z,modelMag_u,modelMag_g,modelMag_r,modelMag_i,modelMag_z,modelMagErr_u,modelMagErr_g,modelMagErr_r,modelMagErr_i,modelMagErr_z) FROM STDIN WITH DELIMITER ',';
    #Get original columns from CLU
    CLU_name = f['data']['name']
    CLU_ra = f['data']['ra']
    CLU_dec = f['data']['dec']
    CLU_z = f['data']['z']
    CLU_dm = f['data']['dm_kin']
    CLU_SDSS_rMag = f['data']['modelmag_r']
    CLU_a = f['data']['a']
    CLU_b2a = f['data']['CLU_ratio_b2a']
    CLU_pa = f['data']['CLU_pa']
#    CLU_petro_rMag = f['data']['petromag_r']
#    CLU_SDSS_rMagErr = f['data']['petromagerr_r']
#    CLU_petro_rMagErr = f['data']['petromagerr_r']
    #Calculate distance in Mpc to each CLU galaxy
    CLU_dist = dm_to_Mpc(CLU_dm)
    #Get columns filtered by CLU_min_dist and CLU_max_dist
    CLU_filter = np.where((min_dist<CLU_dist)&(CLU_dist<max_dist))
    CLU_name = CLU_name[CLU_filter]
    CLU_ra = CLU_ra[CLU_filter]
    CLU_dec = CLU_dec[CLU_filter]
    CLU_z = CLU_z[CLU_filter]
    CLU_dm = CLU_dm[CLU_filter]
    CLU_SDSS_rMag = CLU_SDSS_rMag[CLU_filter]
#    CLU_petro_rMag = CLU_petro_rMag[CLU_filter]
#    CLU_SDSS_rMagErr = CLU_SDSS_rMagErr[CLU_filter]
#    CLU_petro_rMagErr = CLU_petro_rMagErr[CLU_filter]
    CLU_dist = CLU_dist[CLU_filter]
    CLU_a = CLU_a[CLU_filter]
    CLU_b2a = CLU_b2a[CLU_filter]
    CLU_pa = CLU_pa[CLU_filter]
    d = {'CLU_name': CLU_name, 'CLU_ra': CLU_ra, 'CLU_dec': CLU_dec, 'CLU_dist_Mpc': CLU_dist, 'CLU_z': CLU_z,
         'CLU_SDSS_rMag': CLU_SDSS_rMag, 'CLU_a': CLU_a, 'CLU_b2a': CLU_b2a, 'CLU_pa': CLU_pa}#,'CLU_petro_rMag':CLU_petro_rMag,'CLU_SDSS_rMagErr':CLU_SDSS_rMagErr,'CLU_petro_rMagErr':CLU_petro_rMagErr}
    CLU_table = pd.DataFrame(d)
    return CLU_table


def table2cat(table):
    CLU_ra = table['CLU_ra']
    CLU_dec = table['CLU_dec']
    CLU_cat = coordinates.SkyCoord(ra = CLU_ra*u.degree,dec = CLU_dec*u.degree)
    return CLU_cat


def dm_to_Mpc(dm):
    '''Converts dm to distance in Mpc'''
    #return ((10**(1+(dm/5.0)))*u.pc).to(u.Mpc) #-----> Use astropy!!!
    return (coordinates.Distance(distmod = dm)).to(u.Mpc)


def select_circle(cat, ra0, dec0, radius):
    co0 = regions.PixCoord(ra0, dec0)
    cir = regions.circle.CirclePixelRegion(co0, radius)
    catpix = regions.PixCoord(cat.ra, cat.dec)
    return cat[cir.contains(catpix)]
