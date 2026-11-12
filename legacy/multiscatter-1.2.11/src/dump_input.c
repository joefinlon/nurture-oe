/* dump_input.c -- Save input data to a file when a problem is found,
   so that it can be diagnosed

   Copyright (C) 2007  University of Reading

   Author: Robin Hogan <r.j.hogan@reading.ac.uk>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/


#include <math.h>

#include "ms.h"

/* Return 1 if any of the n elements of x contains "not a number", 0
   otherwise.  This is useful to check the inputs or outputs of the
   algorithm. */
int
ms_isnan(int n, const ms_real* x) {
  int i;
  for (i = 0; i < n; i++) {
    if (isnan(x[i])) {
      return 1;
    }
  }
  return 0;
}

/* Write an "input" file from the variables provided; this is useful
   when the multiscatter algorithm is embedded within a larger
   retrieval algorithm, a particular profile has produced NaN outputs
   and it is necessary to know what the inputs were. */
int
ms_dump_input(
    const char* filename,
    /* Input data */
    int n,                    /* Number of input gates */
    ms_config *config,        /* Configuration information */
    ms_instrument instrument, /* Structure containing instrument variables */
    ms_surface surface,       /* Surface scattering variables */
    const ms_real *range,     /* Height of each range gate, metres */
    const ms_real *radius,    /* Cloud/aerosol equivalent radius, microns */
    const ms_real *ext,       /* Total extinction coefficient, m-1 */
    const ms_real *ssa,       /* Total single-scatter albedo */
    const ms_real *g,         /* Total asymmetry factor */
    const ms_real *ext_bscat_ratio, /* Cloud/aerosol ext./bscat. ratio, sr */
    const ms_real *ext_air,   /* Air ext. coefft., m-1 (NULL for vacuum) */
    const ms_real *ssa_air,   /* Air single-scatter albedo */
    const ms_real *droplet_fraction, /* Fraction of extinction from droplets */
    const ms_real *pristine_ice_fraction,/* Fraction of ext from pristine ice */
    /* Adjoint inputs */
    const ms_real *bscat_AD,
    const ms_real *bscat_air_AD)
{
  FILE* file;
  int i;

  if (filename) {
    file = fopen(filename, "w");
  }
  else {
    file = fopen("nandump.in", "w");
  }
  if (!file) {
    return MS_FAILURE;
  }

  fprintf(file, "%d %g %g %g", n, instrument.wavelength, 
	  instrument.altitude, instrument.rho_transmitter);
  for (i = 0; i < instrument.nfov; i++) {
    fprintf(file, " %g", instrument.rho_receiver[i]);
  }
  fprintf(file, "\n");

  for (i = 0; i < n; i++) {
    fprintf(file, "%g %g %g %g", range[i], ext[i], radius[i],
	    ext_bscat_ratio[i]);
    if (ext_air) {
      fprintf(file, " %g", ext_air[i]);
    }
    else {
      fprintf(file, " 0");
    }
    if (ssa && g) {
      fprintf(file, " %g %g", ssa[i], g[i]);
    }
    else {
      fprintf(file, " 0 0");
    }
    if (ssa_air) {
      fprintf(file, " %g", ssa_air[i]);
    }
    else {
      fprintf(file, " 0");
    }
    if (droplet_fraction && pristine_ice_fraction) {
      fprintf(file, " %g %g", droplet_fraction[i],
	      pristine_ice_fraction[i]);
    }
    else {
      fprintf(file, " 0 0");
    }
    if (bscat_AD) {
      int ifov;
      for (ifov = 0; ifov < instrument.nfov; ifov++) {
	fprintf(file, " %g", bscat_AD[i + n*ifov]);
      }
    }
    if (bscat_air_AD) {
      int ifov;
      for (ifov = 0; ifov < instrument.nfov; ifov++) {
	fprintf(file, " %g", bscat_air_AD[i + n*ifov]);
      }
    }
    fprintf(file, "\n");
  }
  fclose(file);
  return MS_SUCCESS;
}
