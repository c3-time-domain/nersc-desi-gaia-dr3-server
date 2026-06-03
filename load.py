import sys
import os
import pathlib

import psycopg
from astropy.table import Table


cols = [ 'source_id', 'ra', 'dec', 'ra_error', 'dec_error', 'phot_g_mean_mag', 'phot_g_mean_flux_over_error',
         'phot_bp_mean_mag', 'phot_bp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_rp_mean_flux_over_error',
         'pm', 'pmra', 'pmdec', 'classprob_dsc_combmod_quasar', 'classprob_dsc_combmod_galaxy',
         'classprob_dsc_combmod_star' ]

con = psycopg.connect( host='postgres', dbname='gaia_dr3_somecols', user='postgres', password=os.getenv('PGPASSWORD') )
cursor = con.cursor()

n = 0
for fname in pathlib.Path( '/gaia-dr3/fits' ).glob( 'GaiaSource*.fits' ):
    sys.stderr.write( f"Reading {fname.name}\n" )
    data = Table.read( fname )
    with cursor.copy( f"COPY gaia_dr3({",".join(cols)}) FROM STDIN" ) as copier:
        for row in data:
            datarow = [ row[c.upper()] for c in cols ]
            copier.write_row( datarow )

    con.commit()

    n += 1
    sys.stderr.write( f"Finished {n} files\n" )

con.close()


