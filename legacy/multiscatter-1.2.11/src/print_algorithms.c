/* print_algorithms.c -- Print the multiscatter algorithms used

   Copyright (C) 2004-2010  University of Reading

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

#include <stdio.h>
#include <math.h>

#include "multiscatter.h"

/* Print the characteristics of the algorithms that will be run to the
   specified file (e.g. stderr) */
void
ms_print_algorithms(ms_config config,
		    ms_instrument instrument,
		    const ms_real* range, 
		    int use_isotropic_pp,
		    int separate_bscat_air,
		    int output_adjoint,
		    int output_jacobian,
		    int calc_jacobian,
		    FILE* file)
{
  int i;
  fprintf(file, "Instrument characteristics:\n");
  fprintf(file, "   Wavelength: %g m\n", instrument.wavelength);
  fprintf(file, "   Gaussian transmitter pattern 1/e half-width: %g radians\n",
	  instrument.rho_transmitter);
  fprintf(file, "   Number of fields-of-view: %d\n", instrument.nfov);
  if (instrument.receiver_type == MS_TOP_HAT) {
    fprintf(file, "   Top-hat receiver pattern half-width(s) in radians:");
    }
  else {
    fprintf(file, "   Gaussian receiver pattern 1/e half-width(s) in radians:");
  }
  for (i = 0; i < instrument.nfov; i++) {
    fprintf(file, " %g", instrument.rho_receiver[i]);
  }
  if (config.options & MS_ANNULAR_DETECTORS) {
    if (instrument.nfov > 2) {
      fprintf(file, "\n   Detectors 2-%d are annular", instrument.nfov);
    }
    else {
      fprintf(file, "\n   Detector 2 is annular");
    }
  }

  fprintf(file, "\n   Altitude: %g m\n", instrument.altitude);
  fprintf(file, "   Distance to first gate: %g m\n", 
	  fabs(instrument.altitude-range[0]));
  fprintf(file, "   Receiver full-width footprint at first gate: %g m\n",
	  fabs(instrument.altitude-range[0])*instrument.rho_receiver[0]*2.0);
  if (!output_jacobian) {
    fprintf(file, "The output from the following algorithm(s) are summed:\n");
  }
  else if (calc_jacobian == 1) {
    fprintf(file, "The Jacobian is calculated from the following algorithm(s):\n");
    }
  else {
    fprintf(file, "The Jacobian with respect to extinction is calculated from the following algorithm(s):\n");
  }
  
  switch (config.small_angle_algorithm) {
  case MS_SINGLE_AND_SMALL_ANGLE_NONE:
    break;
  case MS_SINGLE_SCATTERING:
    fprintf(file, "   Single scattering\n");
    break;
  case MS_SMALL_ANGLE_PVC_EXPLICIT:
    fprintf(file, "   Small-angle scattering using Eloranta-like explicit method\n");
    fprintf(file, "      taken to %d orders of scattering\n", config.max_scattering_order);
    break;
  case MS_SMALL_ANGLE_PVC_FAST:
    fprintf(file, "   Small-angle scattering using O(N) photon variance-covariance (PVC) method\n");
    break;
  case MS_SMALL_ANGLE_PVC_ORIGINAL:
    fprintf(file, "   Small-angle scattering using O(N^2) photon variance-covariance (PVC) method\n");
    break;
  default:
    fprintf(file, "   UNKNOWN SMALL-ANGLE ALGORITHM WITH CODE %d\n",
	    config.small_angle_algorithm);
  }

  if (config.small_angle_algorithm > MS_SINGLE_SCATTERING) {
    if (use_isotropic_pp 
	|| (config.small_angle_algorithm == MS_SMALL_ANGLE_PVC_EXPLICIT)) {
      fprintf(file, "      Phase function assumed to be isotropic near 180 degrees\n");
    }
    else {
      fprintf(file, "      Using anisotropic phase function near 180 degrees according to specified\n"
	      "      droplet and pristine ice fractions\n");
    }
  }

  switch (config.wide_angle_algorithm) {
  case MS_WIDE_ANGLE_NONE:
    break;
  case MS_WIDE_ANGLE_TDTS_FORWARD_LOBE:
    fprintf(file, "   Wide-angle scattering using time-dependent two-stream (TDTS) method assuming a forward lobe\n");
    break;
  case MS_WIDE_ANGLE_TDTS_NO_FORWARD_LOBE:
    fprintf(file, "   Wide-angle scattering using time-dependent two-stream (TDTS) method assuming no forward lobe\n");
    break;
  default:
    fprintf(file, "   UNKNOWN WIDE-ANGLE ALGORITHM WITH CODE %d\n",
	    config.wide_angle_algorithm);
  }

  if (config.wide_angle_algorithm != MS_WIDE_ANGLE_NONE) {
    fprintf(file, "      Coherent backscatter enhancement: %g\n",
	    config.coherent_backscatter_enhancement);
  }

  if (config.small_angle_algorithm == MS_SINGLE_AND_SMALL_ANGLE_NONE
      && config.wide_angle_algorithm == MS_WIDE_ANGLE_NONE) {
    fprintf(file, "   None\n");
  }
  
  if (separate_bscat_air) {
    fprintf(file, "The backscatter due to particulates and air will be reported separately (e.g. for HSRL or Raman)\n");
  }
  else {
    fprintf(file, "The backscatter due to particulates and air will be summed\n");
  }
  
  if (output_adjoint) {
    fprintf(file, "The adjoint will also be calculated\n");
  }
}
