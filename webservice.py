# I have no real idea why matplotlib is getting imported, but something
#   below does, and upon import it is trying to create a config directory,
#   which I consider to be highly dysfunctional, but what can you do.
import pathlib
mplconfigdir = pathlib.Path( '/tmp/matplotlib' )
mplconfigdir.mkdir( exist_ok=True )
import os
os.environ['MPLCONFIGDIR'] = '/tmp/matplotlib'

import sys
import time
# import math
import logging
import textwrap

import psycopg
from psycopg import sql
# import numpy
# from astropy.table import Table
# import healpy

import flask

_default_log_level = 'WARNING'

app = flask.Flask( __name__, instance_relative_config=True )
app.logger.setLevel( os.getenv( 'GAIA_DR3_SERVER_LOG_LEVEL', _default_log_level ) )
_hndlrs = app.logger.handlers
for hndlr in _hndlrs:
    app.logger.removeHandler( hndlr )
_hndlr = logging.StreamHandler( sys.stderr )
_formatter = logging.Formatter( '[%(asctime)s - %(levelname)s] - %(message)s',
                                datefmt='%Y-%m-%d %H:%M:%S' )
_hndlr.setFormatter( _formatter )
app.logger.addHandler( _hndlr )


_dbconfig = {}
for fname, param in zip( [ 'pghost', 'pgport', 'pgdbname', 'pguser', 'pgpass' ],
                         [ 'host',   'port',   'dbname',   'user',   'password' ] ):
    with open( f"/secrets/{fname}" ) as ifp:
        _dbconfig[param] = ifp.readline().strip()
_dbconfig['client_encoding'] = 'UTF8'


@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>/<string:maxmag>/<string:minmag>",
            methods=['GET','POST'], strict_slashes=False )
@app.route( "/gaiarect/<string:ra0>/<string:ra1>/<string:dec0>/<string:dec1>",
            methods=['GET','POST'], strict_slashes=False )
def gaiarect( ra0, ra1, dec0, dec1, maxmag=None, minmag=None ):
    global _dbconfig, app

    t0 = time.monotonic()
    pid = os.getpid()
    app.logger.info( f"PID {pid} got request to {flask.request.base_url}" )

    try:

        # Validate inputs
        try:
            ra0 = float(ra0)
            ra1 = float(ra1)
            dec0 = float(dec0)
            dec1 = float(dec1)
            maxmag = None if maxmag is None else float(maxmag)
            minmag = None if minmag is None else float(minmag)
        except Exception as ex:
            app.logger.error( ex )
            return "Error converting ra/dec values to float", 500

        # Make sure they are sorted right
        if dec0 > dec1:
            tmp = dec0
            dec0 = dec1
            dec1 = tmp
        if ra0 > ra1:
            tmp = ra0
            ra0 = ra1
            ra1 = tmp

        # ...more validating inputs
        if ( dec0 < -90. ) or ( dec1 > 90. ) or ( ra0 < 0. ) or ( ra1 >= 360. ):
            app.logger.error( f"Error, invalid coordinates ({ra0}, {ra1}, {dec0}, {dec1}); "
                              f"δ must be [-90,90], α must be [0,360)" )
            return "Error, invalid coordinates; δ must be [-90,90], α must be [0,360)", 500

        if ( dec1 > 89.9 ) or ( dec0 < -89.9 ):
            app.logger.error( f"Error, invalid dec ({dec0}, {dec1}): "
                              f"coordinates; δ must be [-89.9,89.9], α must be [0,360)" )
            return "Error, currently can't handle poles (|δ|>89.9)", 500

        # Try to detect ra around 0.  Assume that the ra range is <180°, and if it's
        #   greater than that, it means that we're around 0.  In that case, swap
        #   the two ranges so that ra0 > ra1, e.g. ra0=359° and ra1=1°.
        # cyclic = False
        if ( ra1 - ra0 ) > 180.:
            # cyclic = True
            tmp = ra0
            ra0 = ra1
            ra1 = tmp + 360.

        app.logger.debug( f"PID {pid} has validated input" )

        columns = [ 'source_id', 'ra', 'dec', 'ra_error', 'dec_error',
                    'phot_g_mean_mag', 'phot_g_mean_flux_over_error',
                    'phot_bp_mean_mag', 'phot_bp_mean_flux_over_error',
                    'phot_rp_mean_mag', 'phot_rp_mean_flux_over_error',
                    'pm', 'pmra', 'pmdec',
                    'classprob_dsc_combmod_quasar', 'classprob_dsc_combmod_galaxy', 'classprob_dsc_combmod_star' ]
        nls = { 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'pm', 'classprob_dsc_combmod_quasar' }
        nl = '\n'
        colq = sql.SQL( ", ".join( f'{nl if c in nls else ""}{{{c}}}' for c in columns )
                       ).format( **{ c: sql.Identifier(c) for c in columns } )
        q = sql.SQL( textwrap.dedent(
            """
            SELECT {colq}
            FROM gaia_dr3
            WHERE q3c_poly_query( ra, dec, ARRAY[ {ra0}, {dec0}, {ra0}, {dec1}, {ra1}, {dec1}, {ra1}, {dec0} ] )
            """
        ) ).format( colq=colq, ra0=ra0, dec0=dec0, ra1=ra1, dec1=dec1 )
        if minmag is not None:
            q += sql.SQL( "  AND phot_g_mean_mag>={minmag}\n" ).format( minmag=minmag )
        if maxmag is not None:
            q += sql.SQL( "  AND phot_g_mean_mag<={maxmag}\n" ).format( maxmag=maxmag )
        app.logger.debug( f"Sending query: {q.as_string()}" )

        with psycopg.connect( **_dbconfig ) as con:
            cursor = con.cursor()
            cursor.execute( q )
            rows = cursor.fetchall()

        rval = { c: [] for c in columns }
        for row in rows:
            for col, val in zip( columns, row ):
                rval[col].append( val )

        app.logger.info( f"PID {pid} returning {len(rval['ra'])} stars after {time.monotonic()-t0:.2f}s" )
        return rval

    except Exception as ex:
        app.logger.error( f"Exception after {time.monotonic()-t0:.2f}s : {ex}" )
        raise


@app.route( "/", methods=['GET','POST'], strict_slashes=False )
def root():
    return "Hit /gaiarect/ra0/ra1/dec0/dec1 , optionally adding /maxmag/minmag"



######################################################################
# Old code, from when we read desi healpix files

        # # To make sure that we hit all of the possible overlapping healpix, we need very fine
        # #   sampling at the edges (in case it's a small overlap), and than sampling that's
        # #   roughly the size of a healpix in the middle.  Because healpix will in general
        # #   be tilted relative to ra/dec lines, I'm going to use pixsize/2 as the sampling
        # #   spacing internally.  Externally, pixsize/8., though we should probably do even
        # #   better than that.  (Thought required.)
        # # (Is there a better way to figure out which healpix are overlapped by a rectangle?)

        # ras = []
        # decs = []

        # pixsize = healpy.nside2resol( 32 ) * 180. / math.pi
        # ndecedge = max( 2, int(math.ceil( ( dec1 - dec0 ) / ( pixsize / 8. ) )) )

        # for deci in range(ndecedge+1):
        #     dec = dec0 + deci * ( dec1 - dec0 ) / ndecedge
        #     edgeraonly = False
        #     if ( deci == 0 ) or ( deci == ndecedge-1 ):
        #         pixfracra = 8
        #     elif deci % 4 == 0:
        #         pixfracra = 2
        #     else:
        #         edgeraonly = True

        #     if edgeraonly:
        #         decs.extend( [ dec, dec ] )
        #         ras.extend( [ ra0, ra1 ] )
        #     else:
        #         nra = max( 2, int(math.ceil( ( ra1 - ra0 ) /
        #                                      ( pixsize / pixfracra / math.cos( dec * math.pi / 180. ) )
        #                                     )) )
        #         decs.extend( [ dec for i in range(nra+1) ] )
        #         ras.extend( [ ra0 + i * ( ra1 - ra0 ) / nra for i in range(nra+1) ] )

        # ras = numpy.array( ras )
        # ras[ ras >= 360. ] -= 360.
        # decs = numpy.array( decs )

        # app.logger.debug( f"PID {pid} has made the ra/dec grid\n"
        #                   f"        ra0={ra0}, ra1={ra1}, dec0={dec0}, dec1={dec1}, "
        #                   f"len(ras)={len(ras)}, len(decs)={len(decs)}\n"
        #                   f"        ras={ras}, decs={decs}" )


        # hps = set( healpy.ang2pix( 32, ras, decs, nest=True, lonlat=True ) )
        # # ... I hate that the string value of numpy.int64(6812) is "np.int64(6812)"
        # app.logger.debug( f"PID {pid} for ({ra0:.4f}:{ra1:.4f} , {dec0:.4f}:{dec1:.4f}), "
        #                   f"reading files for healpix: {[int(h) for h in hps]}" )

        # # Make the keywords of the returns the same as what you'd get from NOIRLab Data Lab
        # retval = {
        #     'source_id': [],
        #     'ra': [],
        #     'dec': [],
        #     'ra_error': [],
        #     'dec_error': [],
        #     'phot_g_mean_mag': [],
        #     'phot_g_mean_flux_over_error': [],
        #     'phot_bp_mean_mag': [],
        #     'phot_bp_mean_flux_over_error': [],
        #     'phot_rp_mean_mag': [],
        #     'phot_rp_mean_flux_over_error': [],
        #     'pm': [],
        #     'pmra': [],
        #     'pmdec': [],
        #     'classprob_dsc_combmod_star': [],
        #     'classprob_dsc_combmod_quasar': [],
        #     'classprob_dsc_combmod_galaxy': []
        # }

        # #

        # datadir = pathlib.Path( "/data" )
        # for hp in hps:
        #     app.logger.debug( f"PID {pid} reading healpix-{hp:05d}.fits..." )
        #     t = Table.read( datadir / f"healpix-{hp:05d}.fits" )
        #     app.logger.debug( f"...PID {pid} read healpix-{hp:05d}.fits." )
        #     t = t[ ( t['DEC'] >= dec0 ) & ( t['DEC'] <= dec1 ) ]
        #     if cyclic:
        #         t = t[ ( t['RA'] >= ra0 ) | ( t['RA'] <= ( ra1 - 360. ) ) ]
        #     else:
        #         t = t[ ( t['RA'] >= ra0 ) & ( t['RA'] <= ra1 ) ]
        #     if maxmag is not None:
        #         t = t[ t['PHOT_G_MEAN_MAG'] <= maxmag ]
        #     if minmag is not None:
        #         t = t[ t['PHOT_G_MEAN_MAG'] >= minmag ]
        #         app.logger.debug( f"{len(t)} stars from healpix {hp}" )
        #     for kw in retval.keys():
        #         # Gotta convert to floats because the json
        #         #  encoder doesn't know how to handle numpy float32
        #         t[ kw.upper() ] = t[ kw.upper() ].astype( float )
        #         retval[ kw ].extend( list( t[ kw.upper() ] ) )
        #     app.logger.debug( f"PID {pid} done with healpix-{hp:05d}.fits" )
