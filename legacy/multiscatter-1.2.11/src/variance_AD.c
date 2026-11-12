/* variance_AD.c -- Adjoint of the lateral variance of outgoing photons

   Copyright (C) 2011  University of Reading

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

/* Adjoint of the calculation of the lateral variance of the outgoing
   unscattered and forward-scattered photon distribution */
int
ms_variance_AD(int n, ms_real wavelength, ms_real rho_laser,
	       ms_real inst_altitude, const ms_real *range, 
	       const ms_real *radius, const ms_real *ext, 
	       ms_real *variance_out,
	       ms_real *variance_AD,
	       ms_real *radius_AD, ms_real *ext_AD)
{
  ms_real angle_var = rho_laser*rho_laser;
  ms_real covar = fabs(range[0]-inst_altitude)*angle_var;

  variance_out[0]
    = (range[0]-inst_altitude)*(range[0]-inst_altitude)*angle_var;

  int i;
  for (i = 0; i < n-1; i++) {
    ms_real drange = ms_get_drange(i, n, range);
    ms_real theta = wavelength/(MS_PI*radius[i]);
    variance_out[i+1] = variance_out[i]
      + angle_var*drange*drange + covar*drange;
    covar += angle_var*drange;
    angle_var += 0.5*ext[i]*theta*theta*drange;
  }

  /* ADJOINT CALCULATION */
  ms_real angle_var_AD = 0.0;
  ms_real covar_AD = 0.0;

  for (i = n-2; i >= 0; i--) {
    ms_real drange = ms_get_drange(i, n, range);
    ms_real theta = wavelength/(MS_PI*radius[i]);
    ms_real theta_AD = 0.0;

    /* Need to be careful with the order here so that the *_AD
       variables can be reused */
    ext_AD[i] += angle_var_AD*0.5*theta*theta*drange;
    theta_AD += theta*drange*ext[i]*angle_var_AD;
    radius_AD[i] -= theta_AD*theta/radius[i];

    angle_var_AD += covar_AD*drange + variance_AD[i+1]*drange*drange;
    covar_AD += variance_AD[i+1]*drange;

    variance_AD[i] += variance_AD[i+1];
    variance_AD[i+1] = 0.0;
  }
  return MS_SUCCESS;
}

