"""lidar_ms: lidar multiple-scattering operator for optimal estimation (appendix C).

Solvers `single`, `pvc` (Hogan 2008), `pvc+tdts` (Hogan & Battaglia 2008) and an offline Monte
Carlo reference `mc`. Outputs are apparent backscatter by backscatterer class and polarization.
Depends on NumPy and SciPy only; never imports `nurture_oe` [DECISION D29].
"""

__version__ = "0.1.0"

from . import constants  # noqa: E402,F401
from .schema import (Instrument, McConfig, NearBackscatter, OpticalProfile, OpticalTangent,  # noqa: E402,F401
                     PhaseTables, RangeGrid, Returns, ReturnsTangent, SolverConfig)
from .forward import simulate  # noqa: E402,F401
