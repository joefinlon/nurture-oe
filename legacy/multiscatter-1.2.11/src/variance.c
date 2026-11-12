/* variance.c -- Fast calculation of the lateral variance of outgoing photons

   Copyright (C) 2006  University of Reading

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

/* Calculate the lateral variance of the outgoing unscattered and
   forward-scattered photon distribution */
int
ms_variance(int n, ms_real wavelength, ms_real rho_laser,
	    ms_real inst_altitude, const ms_real *range, 
	    const ms_real *radius, const ms_real *ext,
	    ms_real *variance_out)
{
  ms_real angle_var = rho_laser*rho_laser;
  ms_real pos_var = variance_out[0]
    = (range[0]-inst_altitude)*(range[0]-inst_altitude)*angle_var;
  ms_real covar = fabs(range[0]-inst_altitude)*angle_var;

  int i;
  for (i = 0; i < n-1; i++) {
    ms_real drange = ms_get_drange(i, n, range);
    ms_real theta = wavelength/(MS_PI*radius[i]);
    pos_var += angle_var*drange*drange + covar*drange;
    covar += angle_var*drange;
    angle_var += 0.5*ext[i]*theta*theta*drange;
    variance_out[i+1] = pos_var;
  }
  return MS_SUCCESS;
}

