/* config.c -- Configuration manipulation for multiscatter algorithms

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

#include "ms.h"

void
ms_set_default_config(ms_config* config)
{
  static const ms_config local_config = MS_DEFAULT_CONFIG;
  *config = local_config;
}

void 
ms_set_options(ms_config* config, int new_options)
{
  config->options = new_options;
}
void
ms_add_options(ms_config* config, int new_options)
{
  config->options |= new_options;
}


int
ms_set_max_scattering_order(ms_config* config, int norder)
{
  if (norder > 1) {
    config->max_scattering_order = norder;
    return MS_SUCCESS;
  }
  else {
    return MS_FAILURE;
  }
}


void
ms_set_first_wide_angle_gate(ms_config* config, int first_gate)
{
  if (config) {
    config->first_wide_angle_gate = first_gate;
  }
}
