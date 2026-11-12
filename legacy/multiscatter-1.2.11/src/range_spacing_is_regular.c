/* range_spacing_is_regular.c -- Check even spacing of range gates

   Copyright (C) 2004-2007  University of Reading 

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
#include <stdio.h>

#include "ms.h"

/* Check that range gate spacing is monotonic and that the spacings
   are the same to within the specified tolerance (e.g. 0.05 for 5% or
   0.01 for 1%; return 1 if regular and 0 if irregular. */
int
ms_range_spacing_is_regular(int n, const ms_real* range, ms_real tolerance)
{
  ms_real max_drange = -1.0e37;
  ms_real min_drange = 1.0e37;
  int i;
  for (i = 0; i < n-1; i++) {
    ms_real drange = range[i+1] - range[i];
    if (drange > max_drange) {
      max_drange = drange;
    }
    if (drange < min_drange) {
      min_drange = drange;
    }
  }

  if (max_drange * min_drange < 0.0) {
    /* Range is not monotonic */
    return 0;
  }

  if (tolerance >= 1.0 || tolerance <= 0.0) {
    fprintf(stderr, "Warning: tolerance argument to ms_check_range_spacing should be greater than 0 and less than 1\n");
  }

  if (fabs(max_drange) > fabs(min_drange)*(1.0+tolerance)
      || fabs(min_drange) > fabs(max_drange)*(1.0+tolerance)) {
    return 0;
  }

  return 1;
}
