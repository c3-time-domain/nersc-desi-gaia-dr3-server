# To run this, you have to be in a venv that's got the
#   pypi astro-datalab package installed :
#
# python -mvenv venv
# source venv/bin/activate
# pip install astro-datalab pytest
# decativate
#
# After that, you can source venv/bin/activate to get things going
# (It looks like you explicitly need pytest in the venv for it to
# find things installed in the venv.)
#
# (May also need to pip install requests in the venv if your
# venv creation didn't inherit it.  (venv is mysterious.))

import pytest
import requests

from dl import queryClient
import dl.helpers.utils

@pytest.fixture
def url_verify():
    # This is my dev server
    return ( "https://webap.ls4-nersc-gaia-dr3.development.svc.spin.nersc.org/gaiarect", False )
    # This is the production server for use with LS4
    # return( "https://ls4-gaia-dr3.lbl.gov/gaiarect", True )
    
@pytest.fixture
def sky_range():
    return ( 30., 30.15, -20.15, -20. )

# TODO : come up with a sky range that spans a healpix boundry or two

def test_server( url_verify, sky_range ):
    url, verify = url_verify
    ra0, ra1, dec0, dec1 = sky_range

    res = requests.post( f"{url}/{ra0}/{ra1}/{dec0}/{dec1}", verify=verify )
    assert res.status_code == 200

    data = res.json()
    assert set( data.keys() ) == { 'ra', 'dec', 'ra_error', 'dec_error',
                                   'phot_g_mean_mag', 'phot_g_mean_flux_over_error',
                                   'phot_bp_mean_mag', 'phot_bp_mean_flux_over_error',
                                   'phot_rp_mean_mag', 'phot_rp_mean_flux_over_error',
                                   'pm', 'pmra', 'pmdec',
                                   'classprob_dsc_combmod_star' }
    assert len( data['ra'] ) == 54
    assert min( data['ra'] ) == pytest.approx( 30.00731, abs=2e-5 )
    assert max( data['ra'] ) == pytest.approx( 30.14993, abs=2e-5 )
    assert min( data['dec'] ) == pytest.approx( -20.14340, abs=2e-5 )
    assert max( data['dec'] ) == pytest.approx( -20.00452, abs=2e-5 )
    assert min( data['phot_g_mean_mag'] ) == pytest.approx( 14.48, abs=0.01 )
    assert max( data['phot_g_mean_mag'] ) == pytest.approx( 21.60, abs=0.01 )

def test_minmaxmag( url_verify, sky_range ):
    url, verify = url_verify
    ra0, ra1, dec0, dec1 = sky_range

    res = requests.post( f"{url}/{ra0}/{ra1}/{dec0}/{dec1}/20.", verify=verify )
    assert res.status_code == 200

    data = res.json()
    assert len( data['ra'] ) == 27
    assert min( data['phot_g_mean_mag'] ) == pytest.approx( 14.48, abs=0.01 )
    assert max( data['phot_g_mean_mag'] ) == pytest.approx( 19.93, abs=0.01 )

    res = requests.post( f"{url}/{ra0}/{ra1}/{dec0}/{dec1}/20./18.", verify=verify )
    assert res.status_code == 200
    
    data = res.json()
    assert len( data['ra'] ) == 14
    assert min( data['phot_g_mean_mag'] ) == pytest.approx( 18.03, abs=0.01 )
    assert max( data['phot_g_mean_mag'] ) == pytest.approx( 19.93, abs=0.01 )

    
    
def test_server_rawrap( url_verify ):
    url, verify = url_verify
    ra0 = 359.9
    ra1 = 0.1
    dec0 = 0.
    dec1 = 0.1

    res = requests.post( f"{url}/{ra0}/{ra1}/{dec0}/{dec1}", verify=verify )
    assert res.status_code == 200
    
    data = res.json()
    assert len( data['ra'] ) == 44
    highra = [ r for r in data['ra'] if r > 180. ]
    assert min( highra ) == pytest.approx( 359.90093, abs=2e-5 )
    lowra = [ r for r in data['ra'] if r < 180. ]
    assert max( lowra ) == pytest.approx( 0.09951, abs=2e-5 )

    
def test_failure_modes( url_verify ) :
    url, verify = url_verify

    res = requests.post( f"{url}/this/is/not/float", verify=verify )
    assert res.status_code == 500
    assert res.text == "Error converting ra/dec values to float"

    res = requests.post( f"{url}/20./21./0./0.1", verify=verify )
    assert res.status_code == 500
    assert res.text == "Error, can't handle Δra or Δdec > 55 arcmin"
    
    res = requests.post( f"{url}/20./20.1/-1./0.5", verify=verify )
    assert res.status_code == 500
    assert res.text == "Error, can't handle Δra or Δdec > 55 arcmin"

    res = requests.post( f"{url}/359.5/0.5/0./0.2", verify=verify )
    assert res.status_code == 500
    assert res.text == "Error, can't handle Δra or Δdec > 55 arcmin"
    

def test_server_vs_noirlab( url_verify, sky_range ):
    url, verify = url_verify
    ra0, ra1, dec0, dec1 = sky_range

    res = requests.post( f"{url}/{ra0}/{ra1}/{dec0}/{dec1}", verify=verify )
    assert res.status_code == 200

    data = res.json()

    q = ( f"SELECT ra, dec, ra_error, dec_error, pm, pmra, pmdec, "
          f"       phot_g_mean_mag, phot_g_mean_flux_over_error, "
          f"       phot_bp_mean_mag, phot_bp_mean_flux_over_error, "
          f"       phot_rp_mean_mag, phot_rp_mean_flux_over_error, "
          f"       classprob_dsc_combmod_star "
          f"FROM gaia_dr3.gaia_source "
          f"WHERE ra>={ra0} AND ra<={ra1} AND dec>={dec0} AND dec<={dec1} "
         )
    qresult = queryClient.query( sql=q )

    dldata = dl.helpers.utils.convert( qresult, 'structarray' )

    assert len( dldata['ra'] ) == len( data['ra'] )
    assert min( dldata['ra'] ) == pytest.approx( min( data['ra'] ), abs=2e-6 )
    assert max( dldata['ra'] ) == pytest.approx( max( data['ra'] ), abs=2e-6 )
    assert min( dldata['dec'] ) == pytest.approx( min( data['dec'] ), abs=2e-6 )
    assert max( dldata['dec'] ) == pytest.approx( max( data['dec'] ), abs=2e-6 )

    
