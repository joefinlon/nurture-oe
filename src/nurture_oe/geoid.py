"""Geoid undulation along track (EGM2008 / EGM96) via PROJ geoid grids.

Used to check ``[DECISION D23]``: HALO altitudes are MSL, CloudCube ``alt`` is above the
WGS84 ellipsoid, so ``alt_HALO - alt_CC`` should equal ``-N`` (MSL = ellipsoidal - N).
Needs ``pyproj`` and the ``proj-data`` grids ``us_nga_egm08_25.tif`` / ``us_nga_egm96_15.tif``
(conda-forge). Network grid fetching is not used; a missing grid raises.
"""

import numpy as np

GRIDS = {"EGM2008": "us_nga_egm08_25.tif", "EGM96": "us_nga_egm96_15.tif"}


def geoid_undulation(lat_deg, lon_deg, model="EGM2008"):
    """Geoid height N (m above the WGS84 ellipsoid) at each point.

    Parameters
    ----------
    lat_deg, lon_deg : array-like
    model : {"EGM2008", "EGM96"}

    Returns
    -------
    ndarray
        N in metres; NaN where the input is not finite.
    """
    from pyproj import Transformer

    grid = GRIDS[model]
    t = Transformer.from_pipeline(
        "+proj=pipeline +step +proj=unitconvert +xy_in=deg +xy_out=rad "
        f"+step +proj=vgridshift +grids={grid} +multiplier=1 "
        "+step +proj=unitconvert +xy_in=rad +xy_out=deg")
    lat = np.asarray(lat_deg, float)
    lon = np.asarray(lon_deg, float)
    ok = np.isfinite(lat) & np.isfinite(lon)
    n = np.full(lat.shape, np.nan)
    _, _, h = t.transform(lon[ok], lat[ok], np.zeros(ok.sum()), errcheck=True)
    n[ok] = h
    return n
