"""Fit constants, thresholds and flag bits for `lidar_ms` (appendix C §C.2, §C.3.3).

Every number here is either a published fit constant (H08 = Hogan 2008; HB08 = Hogan & Battaglia
2008), a threshold fixed by appendix C, or a value reproduced from `multiscatter` 1.2.11 for the
`compat` switch. None is an instrument specification.
"""

import numpy as np

# --- general -----------------------------------------------------------------------------------
PI = np.pi
SPEED_OF_LIGHT = 2.99792458e8          # m s-1

# --- mixing (§C.2.2) ---------------------------------------------------------------------------
THETA_MAX_DEFAULT = 0.1                # rad, lobe-forming limit (multiscatter default)
RAYLEIGH_BSCAT_EXT = 3.0 / (8.0 * PI)  # sr-1, molecular backscatter-to-extinction ratio (C code; compat)
COMPAT_AIR_BSCAT_MAX_WAVELENGTH = 1e-6  # m, 1.2.11 sets air backscatter to zero above this (Dep1)

# --- pvc (§C.2.4) ------------------------------------------------------------------------------
FORWARD_FRACTION_THRESHOLD = 1e-12     # step 2: T_a <= 1e-12 T_r  =>  M = 0
LOGMEAN_TOL = 1e-12                    # step 6: |Y- - Y+| <= 1e-12 Y-  =>  beta Y-
COMPAT_ONE_THIRD = 0.333333            # 1.2.11 writes 1/3 as 0.333333 in fast_small_angle.c

# H08 Eq. 36, droplet near-backscatter fit
DROPLET_U1 = 0.3
DROPLET_U2 = 0.5
DROPLET_U0 = 1.0 - DROPLET_U1 - DROPLET_U2   # 0.2
DROPLET_V1_SQ = 16.0                   # V1 = 4
DROPLET_V2_SQ = 0.16                   # V2 = 0.4

# H08 Eq. 37, pristine-ice (Yang et al. 2000) near-backscatter fit
YANG_W = 0.89
YANG_GAMMA0 = 0.038                    # rad

# --- closures (§C.2.6) -------------------------------------------------------------------------
C1_DELTA_ACC_MAX = 0.45                # validity of Hu (2007) closed form
C1_MIN_FOOTPRINT = 3.0                 # m, Hu et al. (2006) reduced accuracy below this footprint
LOBE_MIXTURE_RATIO = 3.0               # flag bit 4
TABLE_TAIL_WEIGHT = 1e-3               # flag bit 5
DELTA_TABLE_MISMATCH = 0.02            # NearBackscatter warning
EXTINGUISHED_TU = 1e-12                # flag bit 6

# --- tdts (§C.2.5) -----------------------------------------------------------------------------
MU1 = 0.5
RANGE_SPACING_TOL = 0.05

# --- derivatives (§C.2.8) ----------------------------------------------------------------------
COMPLEX_STEP = 1e-30

# --- flag bits (§C.3.3) ------------------------------------------------------------------------
FLAG_BITS = {
    "C1_DELTA_RANGE": 0,
    "C1_SMALL_FOOTPRINT": 1,
    "C1_NEGATIVE": 2,
    "C1_UNDER_OTHER_LOBE": 3,
    "LOBE_MIXTURE": 4,
    "TABLE_EXTRAPOLATED": 5,
    "EXTINGUISHED": 6,
    "WIDE_HEURISTIC": 7,
}
FLAG_MEANINGS = tuple(FLAG_BITS)       # ordered by bit
FLAG_MASKS = tuple(1 << b for b in FLAG_BITS.values())


def flag(name: str) -> int:
    """Mask of the named flag bit."""
    return 1 << FLAG_BITS[name]
