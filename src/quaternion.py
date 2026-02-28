"""
Quaternion arithmetic for binary polyhedral groups in SU(2).

Convention: q = (w, x, y, z) with w = scalar, (x,y,z) = vector.
SU(2) correspondence: U = w*I + i*(x*σ1 + y*σ2 + z*σ3).
Unit quaternions |q| = 1 ↔ SU(2).
"""
import numpy as np


def qmul(a, b):
    """Hamilton product of quaternion arrays. Shape (..., 4)."""
    w1, x1, y1, z1 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    w2, x2, y2, z2 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ], axis=-1)


def qconj(a):
    """Quaternion conjugate. For unit quaternions, qconj = inverse."""
    return a * np.array([1, -1, -1, -1])


def qinv(q):
    """Inverse of a unit quaternion: (w, -x, -y, -z)."""
    r = q.copy()
    r[..., 1:] *= -1
    return r


def qnorm(a):
    """Quaternion norm."""
    return np.sqrt(np.sum(a**2, axis=-1, keepdims=True))


def normalize(a):
    """Project to unit quaternions."""
    return a / (qnorm(a) + 1e-30)


ROUND_DIGITS = 12


def qkey(q):
    """Canonical hashable key for a quaternion (rounded tuple)."""
    return tuple(round(float(x), ROUND_DIGITS) for x in q)


def quat_to_su2(q):
    """Convert unit quaternion to SU(2) matrix.
    q = (w, x, y, z) → U = w·I + i·(x·σ₁ + y·σ₂ + z·σ₃).
    """
    w, x, y, z = q
    return np.array([
        [w + 1j*z,  1j*x - y],
        [1j*x + y,  w - 1j*z]
    ])


def quat_to_so3(q):
    """
    Convert unit quaternion to SO(3) rotation matrix.
    q = (w, x, y, z) → 3×3 matrix R.
    Key property: q and -q give the SAME R.
    """
    w, x, y, z = q[..., 0], q[..., 1], q[..., 2], q[..., 3]
    R = np.zeros(q.shape[:-1] + (3, 3))
    R[..., 0, 0] = 1 - 2*(y**2 + z**2)
    R[..., 0, 1] = 2*(x*y - w*z)
    R[..., 0, 2] = 2*(x*z + w*y)
    R[..., 1, 0] = 2*(x*y + w*z)
    R[..., 1, 1] = 1 - 2*(x**2 + z**2)
    R[..., 1, 2] = 2*(y*z - w*x)
    R[..., 2, 0] = 2*(x*z - w*y)
    R[..., 2, 1] = 2*(y*z + w*x)
    R[..., 2, 2] = 1 - 2*(x**2 + y**2)
    return R
