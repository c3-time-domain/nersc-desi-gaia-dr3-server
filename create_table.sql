CREATE TABLE gaia_dr3(
   source_id bigint PRIMARY KEY,
   ra double precision,
   dec double precision,
   ra_error double precision,
   dec_error double precision,
   phot_g_mean_mag real,
   phot_g_mean_flux_over_error real,
   phot_bp_mean_mag real,
   phot_bp_mean_flux_over_error real,
   phot_rp_mean_mag real,
   phot_rp_mean_flux_over_error real,
   pm  real,
   pmra real,
   pmdec real,
   classprob_dsc_combmod_quasar real,
   classprob_dsc_combmod_galaxy real,
   classprob_dsc_combmod_star real
);
CREATE INDEX ix_gaia_dr3_q3c ON gaia_dr3(q3c_ang2ipix(ra, dec));
