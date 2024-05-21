# I have no real idea why matplotlib is getting imported, but something
#   below does, and upon import it is trying to create a config directory,
#   which I consider to be highly dysfunctional, but what can you do.
import pathlib
mplconfigdir = pathlib.Path( '/tmp/matplotlib' )
mplconfigdir.mkdir( exist_ok=True )
import os
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'

import sys
import io
import re
import math
import logging

import numpy
from astropy.table import Table
import healpy

import flask

app = flask.Flask( __name__, instance_relative_config=True )
# app.logger.setLevel( logging.INFO )
app.logger.setLevel( logging.DEBUG )


# TODO : Implement minmag, maxmag

# @app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>",
#             methods=['GET','POST'], strict_slashes=False )
# @app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>/<string:minmag>",
#             methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>",
            methods=['GET','POST'], strict_slashes=False )
def gaiarect( ra0, ra1, dec0, dec1, maxmag=None, minmag=None ):
    try:
        ra0 = float(ra0)
        ra1 = float(ra1)
        dec0 = float(dec0)
        dec1 = float(dec1)
    except Exception as ex:
        app.logger.error( ex )
        return f"Error converting ra/dec values to float", 500

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

    # Pick spots around the edge in hopes that we get all overlapping pixels
    # The area cut above means that no pixel will be entirely inside the rectangle
    ras = numpy.array( [ ra0, ra0+dra/4., ra0+dra/2., ra0+3.*dra/4., ra1,
                         ra0, ra1,
                         ra0, ra1,
                         ra0, ra1,
                         ra0, ra0+dra/4., ra0+dra/2., ra0+3.*dra/4., ra1 ] )
    decs = numpy.array( [ dec0, dec0, dec0, dec0, dec0,
                          dec0+ddec/4., dec0+ddec/4.,
                          dec0+ddec/2., dec0+ddec/2.,
                          dec0+3.*ddec/4., dec0+3.*ddec/4.,
                          dec1, dec1, dec1, dec1, dec1 ] )

    # Try to detect ra around 0
    cyclic = False
    if ( ra1 - ra0 ) > 180.:
        app.logger.debug( "Detected ra around 0" )
        cyclic = True
        tmp = ra0
        ra0 = ra1
        ra1 = tmp
        dra = (360.-ra0) + ra1
        # These aren't evenly spaced, but whatevs
        ras = numpy.array( [ ra0, ra0+(360.-ra0)/2., 0., ra1/2., ra1,
                             ra0, ra1,
                             ra0, ra1,
                             ra0, ra1,
                             ra0, ra0+(360.-ra0)/2., 0., ra1/2., ra1 ] )
        decs = numpy.array( [ dec0, dec0, dec0, dec0, dec0,
                              dec0+ddec/4., dec0+ddec/4.,
                              dec0+ddec/2., dec0+ddec/2.,
                              dec0+3.*ddec/4., dec0+3.*ddec/4.,
                              dec1, dec1, dec1, dec1, dec1 ] )

    app.logger.debug( f"ra0={ra0}, ra1={ra1}, dec0={dec0}, dec1={dec1}, "
                      f"dra={dra}, ddec={ddec}, ras={ras}, decs={decs}" )
    app.logger.debug( f"ras={ras}, decs={decs}" )

    # Punt if too big of a patch was requested
    # (This dra is bad for dec close to the poles....)
    pixsize = healpy.nside2resol( 32 ) * 180. / math.pi
    if ( dra * math.cos( (dec0+dec1)/2. * math.pi / 180. ) > pixsize/2. ) or ( ddec > pixsize/2. ):
        app.logger.error( f"Error: ({ra0:.4f}:{ra1:.4f} , {dec0:.4f}:{dec1:.4f}) is too big a rectangle" )
        return f"Error, can't handle Δra or Δdec > {60.*pixsize/2.:.0f} arcmin", 500

    hps = set( healpy.ang2pix( 32, ras, decs, nest=True, lonlat=True ) )
    app.logger.debug( f"For ({ra0:.4f}:{ra1:.4f} , {dec0:.4f}:{dec1:.4f}), reading files for healpix: {hps}" )

    # Make the keywords of the returns the same as what you'd get from NOIRLab Data Lab
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
        app.logger.debug( f"{len(t)} stars from healpix {hp}" )
        for kw in retval.keys():
            # Gotta convert to floats because the json
            #  encoder doesn't know how to handle numpy float32
            t[ kw.upper() ] = t[ kw.upper() ].astype( float )
            retval[ kw ].extend( list( t[ kw.upper() ] ) )

    app.logger.debug( f"Returning {len(retval['ra'])} stars" )
    return retval

@app.route( "/hello", methods=['GET','POST'], strict_slashes=False )
def hello():
    return "Hello, world"
