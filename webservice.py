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
app.logger.setLevel( logging.INFO )
# app.logger.setLevel( logging.DEBUG )


@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>/<string:minmag>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>",
            methods=['GET','POST'], strict_slashes=False )
def gaiarect( ra0, ra1, dec0, dec1, maxmag=None, minmag=None ):
    try:
        ra0 = float(ra0)
        ra1 = float(ra1)
        dec0 = float(dec0)
        dec1 = float(dec1)
        maxmag = None if maxmag is None else float(maxmag)
        minmag = None if minmag is None else float(minmag)
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

    if ( dec0 < -90. ) or ( dec1 > 90. ) or ( ra0 < 0. ) or ( ra1 >= 360. ):
        return f"Error, invalid coordinates; δ must be [-90,90], α must be [0,360)", 500
    
    if ( dec1 > 89.9 ) or ( dec1 < -89.9 ):
        return f"Error, currently can't handle poles (|δ|>89.9)", 500
    
    # Try to detect ra around 0
    cyclic = False
    if ( ra1 - ra0 ) > 180.:
        cyclic = True
        tmp = ra0
        ra0 = ra1
        ra1 = tmp + 360.
    dra = ra1 - ra0

    # To make sure that we hit all of the possible overlapping healpix, we need very fine
    #   sampling at the edges (in case it's a small overlap), and than sampling that's
    #   roughly the size of a healpix in the middle.  Because healpix will in general
    #   be tilted relative to ra/dec lines, I'm going to use pixsize/2 as the sampling
    #   spacing internally.  Externally, pixsize/8., though we should probably do even
    #   better than that.  (Thought required.)
    # (Is there a better way to figure out which healpix are overlapped by a rectangle?)

    ras = []
    decs = []

    pixsize = healpy.nside2resol( 32 ) * 180. / math.pi
    ndecedge = max( 2, int(math.ceil( ( dec1 - dec0 ) / ( pixsize / 8. ) )) )

    for deci in range(ndecedge+1):
        dec = dec0 + deci * ( dec1 - dec0 ) / ndecedge
        edgeraonly = False
        if ( deci == 0 ) or ( deci == ndecedge-1 ):
            pixfracra = 8
        elif deci % 4 == 0:
            pixfracra = 2
        else:
            egeraonly = True

        if edgeraonly:
            decs.extend( [ dec, dec ] )
            ras.extend( [ ra0, ra1 ] )
        else:
            nra = max( 2, int(math.ceil( ( ra1 - ra0 ) /
                                         ( pixsize / pixfracra / math.cos( dec * math.pi / 180. ) )
                                        )) )
            decs.extend( [ dec for i in range(nra+1) ] )
            ras.extend( [ ra0 + i * ( ra1 - ra0 ) / nra for i in range(nra+1) ] )

    ras = numpy.array( ras )
    ras[ ras >= 360. ] -= 360.
    decs = numpy.array( decs )
                                        
    app.logger.debug( f"ra0={ra0}, ra1={ra1}, dec0={dec0}, dec1={dec1}, "
                      f"len(ras)={len(ras)}, len(decs)={len(decs)}" )
    app.logger.debug( f"ras={ras}, decs={decs}" )

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
            t = t[ ( t['RA'] >= ra0 ) | ( t['RA'] <= ( ra1 - 360. ) ) ]
        else:
            t = t[ ( t['RA'] >= ra0 ) & ( t['RA'] <= ra1 ) ]
        if maxmag is not None:
            t = t[ t['PHOT_G_MEAN_MAG'] <= maxmag ]
        if minmag is not None:
            t = t[ t['PHOT_G_MEAN_MAG'] >= minmag ]
        app.logger.debug( f"{len(t)} stars from healpix {hp}" )
        for kw in retval.keys():
            # Gotta convert to floats because the json
            #  encoder doesn't know how to handle numpy float32
            t[ kw.upper() ] = t[ kw.upper() ].astype( float )
            retval[ kw ].extend( list( t[ kw.upper() ] ) )

    app.logger.debug( f"Returning {len(retval['ra'])} stars" )
    return retval

@app.route( "/", methods=['GET','POST'], strict_slashes=False )
def root():
    return "Hit /gaiarect/ra0/ra1/dec0/dec1 , optionally adding /maxmag/minmag"
