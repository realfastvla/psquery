# Taken from https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/-/blob/master/pdr3/hscReleaseQuery/hscReleaseQuery.py

import json
import urllib.request, urllib.error, urllib.parse
import time
import sys
import os
import os.path

from . import get_coord

API_URL = 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'

if not os.environ.get("HSC_SSP_CAS_PASSWORD"):
    print("need to set env HSC_SSP_CAS_PASSWORD")
    credential = None
else:
    credential = {'account_name': os.environ.get("HSC_SSP_CAS_ID"),
                  'password': os.environ.get("HSC_SSP_CAS_PASSWORD")}


def cone_hsc(radec, radius=5/3600):
    """Query Subaru HSC SSP catalog with (ra, dec) in degrees.
    Radius in degrees.
    Returns ?.
    """

    ra, dec = get_coord(radec, ret="radec")

    sql = f"""
    SELECT
      object_id
      , ra
      , dec
      , g_cmodel_mag
      , g_cmodel_magerr
      , g_extendedness_value
      , r_cmodel_mag
      , r_cmodel_magerr
      , r_extendedness_value
      , i_cmodel_mag
      , i_cmodel_magerr
      , i_extendedness_value
      , z_cmodel_mag
      , z_cmodel_magerr
      , z_extendedness_value
      , y_cmodel_mag
      , y_cmodel_magerr
      , y_extendedness_value
      , photoz_best
      , photoz_err68_min
      , photoz_err68_max
    FROM pdr3_wide.forced
    JOIN pdr3_wide.photoz_demp USING (object_id)
    WHERE isprimary
    AND conesearch(coord, {ra}, {dec}, {radius*3600});
    """

    job = None

    try:
        job = submitJob(credential, sql)
        blockUntilJobFinishes(credential, job['id'])
        download(credential, job['id'], sys.stdout.buffer)
        deleteJob(credential, job['id'])
    except urllib.error.HTTPError as e:
        if e.code == 401:
            print('invalid id or password.', file=sys.stderr)
        if e.code == 406:
            print(e.read(), file=sys.stderr)
        else:
            print(e, file=sys.stderr)
    except QueryError as e:
        print(e, file=sys.stderr)
    except KeyboardInterrupt:
        if job is not None:
            jobCancel(credential, job['id'])
        raise


def submitJob(credential, sql, out_format='csv', release_version='pdr3-citus-columnar'):
    """ Submits a SQL query with credentials.
    Assumes no email needed and no syntax check.
    release_version can be: pdr1 pdr2 pdr2-citus pdr3 pdr3-citus pdr3-citus-columnar
    """

    url = API_URL + 'submit'
    catalog_job = {
        'sql'                     : sql,
        'out_format'              : out_format,
        'include_metainfo_to_body': True,
        'release_version'         : release_version
    }
    postData = {'credential': credential, 'catalog_job': catalog_job, 'nomail': True, 'skip_syntax_check': True}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def httpJsonPost(url, data, version = '20190514.1'):
    data['clientVersion'] = version
    postData = json.dumps(data)
    return httpPost(url, postData, {'Content-type': 'application/json'})


def httpPost(url, postData, headers):
    req = urllib.request.Request(url, postData.encode('utf-8'), headers)
    res = urllib.request.urlopen(req)
    return res


def jobStatus(credential, job_id):
    url = API_URL + 'status'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    job = json.load(res)
    return job


def jobCancel(credential, job_id):
    url = API_URL + 'cancel'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


def blockUntilJobFinishes(credential, job_id):
    max_interval = 5 * 60 # sec.
    interval = 1
    while True:
        time.sleep(interval)
        job = jobStatus(credential, job_id)
        if job['status'] == 'error':
            raise QueryError('query error: ' + job['error'])
        if job['status'] == 'done':
            break
        interval *= 2
        if interval > max_interval:
            interval = max_interval


def download(credential, job_id, out):
    url = API_URL + 'download'
    postData = {'credential': credential, 'id': job_id}
    res = httpJsonPost(url, postData)
    bufSize = 64 * 1<<10 # 64k
    while True:
        buf = res.read(bufSize)
        out.write(buf)
        if len(buf) < bufSize:
            break


def deleteJob(credential, job_id):
    url = API_URL + 'delete'
    postData = {'credential': credential, 'id': job_id}
    httpJsonPost(url, postData)


class QueryError(Exception):
    pass
