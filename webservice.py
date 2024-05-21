import sys
import io
import re
import pathlib
import math

import numpy
from astropy.table import Table
import healpy

import flask

app = flask.Flask( __name__, instance_relative_config=True )
# app.logger.setLevel( logging.INFO )
app.logger.setLevel( logging.DEBUG )


@app.route( "/gaiarect/<float:ra0>/<float:ra1>/<float:dec0>/<float:dec1>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<float:ra0>/<float:ra1>/<float:dec0>/<float:dec1>/<float:maxmag>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<float:ra0>/<float:ra1>/<float:dec0>/<float:dec1>/<float:maxmag>/<float:minmag>",
            methods=['GET','POST'], strict_slashes=False )
def gaiarect( ra0, ra1, dec0, dec1, maxmag=None, minmag=None ):
    if dec0 > dec1:
        tmp = dec0
        dec0 = dec1
        dec1 = tmp
    ddec = dec1 - dec0
    if ra0 > ra1:
        tmp = ra0
        ra0 = ra1
        ra1 = tmp
    dra = ra1 - ra0
    # Try to detect ra around 0
    cyclic = False
    if ( ra1 - ra0 ) > 180.:
        tmp = ra0
        ra0 = ra1
        ra1 = ra0
        cyclic = True

    # Punt if too big of a patch was requested
    # (This dra is bad for dec close to the poles....)
    pixsize = healpy.nside2resol( 32 ) * 180. / math.pi
    if ( dra * math.cos( (dec0+dec1)/2. * math.pi / 180. ) > pixsize/2. ) or ( ddec > pixsize/2. ):
        app.logger.error( f"Error: ({ra0:.4f}:{ra1:.4f} , {dec0:.4f}:{dec1:.4f}) is too big a rectangle" )
        return f"Error, can't handle Δra or Δdec > {pixsize/2.*206265:.0f} arcsec", 500

    # Pick spots around the edge in hopes that we get all overlapping pixels
    # The area cut above means that no pixel will be entirely inside the rectangle
    ras = numpy.array( [ ra0, ra0+dra/4., ra0+dra/2., ra0+3.,*dra/4., ra1,
                         ra0, ra1,
                         ra0, ra1,
                         ra0, ra1,
                         ra0, ra0+dra/4., ra0+dra/2., ra0+3.,*dra/4., ra1, ] )
    decs = numpy.array( [ dec0, dec0, dec0, dec0, dec0,
                          dec0+ddec/4., dec0+ddec/4.,
                          dec0+ddec/2., dec0+ddec/2.,
                          dec0+3.*ddec/4., dec0+3.*ddec/4.,
                          dec1, dec1, dec1, dec1, dec1 ] )

    hps = set( healpy.ang2pix( 32, ras, decs, nest=True, lonlat=True ) )
    app.logger.debug( f"For ({ra0:.4f}:{ra1:.4f} , {dec0:.4f}:{dec1:.4f}), reading files for healpix: {hps}" )
               
    # Make the returns the same as what you'd get from NOIRLab Data Lab
    retval = {
        'ra': [],
        'dec': [],
        'ra_error': [],
        'dec_error': [],
        'phot_g_mean_mag': [],
        'phot_g_mean_flux_over_error': [],
        'phot_bp_mean_mag': [],
        'phot_bp_mean_flux_over_error': [],
        'phot_rp_mean_mag': [],
        'phot_rp_mean_flux_over_error': [],
        'pm': [],
        'pmra': [],
        'pmdec': [],
        'classprob_dsc_combmod_star': []
    }
    
    datadir = pathlib.Path( "/data" )
    for hp in hps:
        t = Table.read( datadir / f"healpix-{hp:05d}.fits" )
        t = t[ ( t['DEC'] >= dec0 ) & ( t['DEC'] <= dec1 ) ]
        if cyclic:
            t = t[ ( t['RA'] >= ra0 ) | ( t['RA'] <= ra1 ) ]
        else:
            t = t[ ( t['RA'] >= ra0 ) & ( t['RA'] <= ra1 ) ]
        for kw in retval.keys():
            retval[ kw ].extend( list( t[ kw.upper() ] ) )

    return retval
