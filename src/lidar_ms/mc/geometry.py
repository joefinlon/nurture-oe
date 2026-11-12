"""Vector algebra, scattering plane and Stokes rotation for solver `mc` (appendix C §C.2.7).

Vectors are stored as structure of arrays: shape (3, n), one column per photon packet. Angles of
rotation are never formed explicitly: the kernels carry (cos 2psi, sin 2psi).

Conventions (fixed by tests P1-P3): psi is the angle from the old reference axis e to the new one,
counted right-handed about the propagation direction; a Stokes vector (I, Q, U, V) referred to e
describes linear polarization at angle chi from e (right-handed about the direction) as
(1, cos 2chi, sin 2chi, 0); L(psi) refers it to the rotated axis.
"""
import numpy as np

DEGENERATE = 1e-12          # |u x u'| below this: scattering plane undefined (§C.2.7)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a, b):
    return np.stack([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])


def norm(a):
    return np.sqrt(dot(a, a))


def normalize(a):
    return a / norm(a)


def rotation_matrix(psi):
    """L(psi), shape (..., 4, 4) (§C.2.7)."""
    psi = np.asarray(psi, dtype=float)
    c, s = np.cos(2 * psi), np.sin(2 * psi)
    L = np.zeros(psi.shape + (4, 4))
    L[..., 0, 0] = L[..., 3, 3] = 1.0
    L[..., 1, 1] = L[..., 2, 2] = c
    L[..., 1, 2] = s
    L[..., 2, 1] = -s
    return L


def rotate_stokes(S, c2, s2):
    """Apply L(psi) given cos 2psi, sin 2psi; S has shape (4, n)."""
    q = c2 * S[1] + s2 * S[2]
    u = -s2 * S[1] + c2 * S[2]
    return np.stack([S[0], q, u, S[3]])


def double_angle(c, s):
    """(cos 2psi, sin 2psi) from an unnormalized (cos psi, sin psi)."""
    h = c * c + s * s
    ok = h > 0
    h = np.where(ok, h, 1.0)
    return np.where(ok, (c * c - s * s) / h, 1.0), np.where(ok, 2 * c * s / h, 0.0)


def scattering_plane(u, e, v):
    """Geometry of a scattering from direction u (reference axis e) into direction v.

    Returns
    -------
    cos_t, sin_t : (n,)  cosine and sine of the scattering angle (sin_t = |u x v| >= 0)
    c2, s2 : (n,)        cos 2psi, sin 2psi: rotation from e into the scattering plane
    l_out : (3, n)       reference axis after scattering (in the plane, perpendicular to v)

    Notes
    -----
    n = u x v / |u x v|, l_in = n x u, l_out = n x v,
    psi = atan2((e x l_in) . u, e . l_in). If |u x v| < 1e-12, l_in = e (n = u x e) and psi = 0.
    """
    nvec = cross(u, v)
    sin_t = norm(nvec)
    cos_t = dot(u, v)
    deg = sin_t < DEGENERATE
    n_deg = cross(u, e)
    nvec = np.where(deg, n_deg, nvec / np.where(deg, 1.0, sin_t))
    l_in = cross(nvec, u)
    c = dot(e, l_in)
    s = dot(cross(e, l_in), u)
    c = np.where(deg, 1.0, c)
    s = np.where(deg, 0.0, s)
    c2, s2 = double_angle(c, s)
    l_out = cross(nvec, v)
    return cos_t, sin_t, c2, s2, l_out


def analyzer_rotation(l_out, d, p):
    """cos 2psi2, sin 2psi2 of the rotation from l_out to the analyzer axis p about d (§C.2.7)."""
    c = dot(l_out, p)
    s = dot(cross(l_out, p), d)
    return double_angle(c, s)


def perpendicular_axis(a, u):
    """Unit vector along a - (a . u) u; a is (3,) or (3, n)."""
    a = np.asarray(a, dtype=float)
    if a.ndim == 1:
        a = a[:, None]
    p = a - dot(a, u) * u
    return p / norm(p)


def new_direction(u, e, cos_t, phi):
    """Direction at polar angle acos(cos_t) about u and azimuth phi counted from e towards u x e."""
    sin_t = np.sqrt(np.maximum(0.0, (1.0 - cos_t) * (1.0 + cos_t)))
    w = np.cos(phi) * e + np.sin(phi) * cross(u, e)
    v = cos_t * u + sin_t * w
    return v / norm(v)


def orthonormal_pair(d):
    """Two unit vectors a, b with (a, b, d) right-handed orthonormal; d is (3, n)."""
    ref = np.where(np.abs(d[0]) < 0.9, 1.0, 0.0)
    helper = np.stack([ref, 1.0 - ref, np.zeros_like(ref)])
    a = normalize(helper - dot(helper, d) * d)
    return a, cross(d, a)
