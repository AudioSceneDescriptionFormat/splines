"""Quaternions and unit-quaternion splines."""
import math as _math

import numpy as _np
from . import _check_param


class Quaternion:
    """A very simple quaternion class.

    This is the base class for the more relevant `UnitQuaternion` class.

    """

    __slots__ = '_scalar', '_vector'

    def __new__(cls, scalar, vector):
        obj = super().__new__(cls)
        obj._scalar = scalar
        x, y, z = vector
        obj._vector = x, y, z
        return obj

    @property
    def scalar(self):
        return self._scalar

    @property
    def vector(self):
        return self._vector

    def __repr__(self):
        name = type(self).__name__
        return f'{name}(scalar={self.scalar}, vector={self.vector})'

    def __eq__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        return self.xyzw == other.xyzw

    def __mul__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        if isinstance(self, UnitQuaternion) and isinstance(other, UnitQuaternion):
            result_type = UnitQuaternion
        else:
            result_type = Quaternion
        a1 = self.scalar
        b1, c1, d1 = self.vector
        a2 = other.scalar
        b2, c2, d2 = other.vector
        return Quaternion.__new__(
            result_type,
            # TODO: properly derive this
            a1*a2 - b1*b2 - c1*c2 - d1*d2,
            (
                a1*b2 + b1*a2 + c1*d2 - d1*c2,
                a1*c2 - b1*d2 + c1*a2 + d1*b2,
                a1*d2 + b1*c2 - c1*b2 + d1*a2,
            )
        )

    def __neg__(self):
        x, y, z = self.vector
        return Quaternion.__new__(type(self), -self.scalar, (-x, -y, -z))

    def conjugate(self):
        x, y, z = self.vector
        return Quaternion.__new__(type(self), self.scalar, (-x, -y, -z))

    def normalize(self):
        norm = self.norm
        x, y, z, w = self.xyzw
        return UnitQuaternion.from_unit_xyzw(
            (x / norm, y / norm, z / norm, w / norm))

    def dot(self, other):
        """Dot product of two quaternions.

        This is the 4-dimensional dot product, yielding a scalar result.
        This operation is commutative.

        Note that this is different from the quaternion multiplication
        (``q1 * q2``), which produces another quaternion
        (and is non-commutative).

        """
        # NB: math.prod() is available since Python 3.8
        prod = _np.multiply.reduce
        return sum(map(prod, zip(self.xyzw, other.xyzw)))

    @property
    def norm(self):
        x, y, z, w = self.xyzw
        return _math.sqrt(x**2 + y**2 + z**2 + w**2)

    @property
    def xyzw(self):
        x, y, z = self.vector
        return x, y, z, self.scalar

    @property
    def wxyz(self):
        x, y, z = self.vector
        return self.scalar, x, y, z


class UnitQuaternion(Quaternion):
    """Unit quaternion."""

    __slots__ = ()

    def __new__(cls, *args, **kwargs):
        raise TypeError('Use UnitQuaternion.from_*() to create a UnitQuaternion')

    @classmethod
    def from_axis_angle(cls, axis, angle):
        """Create a unit quaternion from a rotation axis and angle.

        :param axis: Three-component rotation axis. This will be normalized.
        :param angle: Rotation angle in radians.

        """
        x, y, z = axis
        norm = _math.sqrt(x**2 + y**2 + z**2)
        s = angle / (2 * norm)
        return cls.exp_map((s * x, s * y, s * z))

    @classmethod
    def from_unit_xyzw(cls, xyzw):
        """Create a unit quaternion from another unit quaternion.

        :param xyzw: Components of a unit quaternion (scalar last).
            This will *not* be normalized, it must already have unit length.

        """
        x, y, z, w = xyzw
        return super().__new__(cls, w, (x, y, z))

    def __pow__(self, alpha):
        if not _np.isscalar(alpha):
            return _np.array([self.__pow__(a) for a in alpha])
        if self.scalar == 1:
            return super().__new__(UnitQuaternion, self.scalar, self.vector)
        return UnitQuaternion.from_axis_angle(self.axis, alpha * self.angle)

    inverse = Quaternion.conjugate
    """Multiplicative inverse.

    For unit quaternions, this is the same as `conjugate()`.

    """

    @classmethod
    def exp_map(cls, value):
        """Exponential map from :math:`R^3` to quaternions.

        This is the inverse operation to `log_map()`.

        :param value: Element of the tangent space at the quaternion identity.
        :param type: 3-tuple
        :returns: Corresponding unit quaternion.

        """
        x, y, z = value
        norm = _math.sqrt(x**2 + y**2 + z**2)
        if norm == 0:
            zero, one = norm, norm + 1  # to get appropriate numeric types
            return super().__new__(cls, one, (zero, zero, zero))
        s = _math.sin(norm) / norm
        return super().__new__(cls, _math.cos(norm), (s * x, s * y, s * z))

    def log_map(self):
        """Logarithmic map from quaternions to :math:`R^3`.

        :returns: Corresponding vector in the tangent space at the
            quaternion identity.
        :rtype: 3-tuple

        """
        length = self.angle / 2
        if self.scalar == 1:
            zero = 0 * self.scalar  # to get appropriate numeric type
            return zero, zero, zero
        x, y, z = self.axis
        return x * length, y * length, z * length

    @property
    def axis(self):
        assert self.scalar <= 1
        # NB: This is the same as sqrt(x**2 + y**2 + z**2)
        norm = _math.sqrt(1 - self.scalar**2)
        x, y, z = self.vector
        return x / norm, y / norm, z / norm

    @property
    def angle(self):
        return 2 * _math.acos(self.scalar)

    def rotate_vector(self, v):
        rotated = self * Quaternion(0, v) * self.inverse()
        return rotated.vector


def slerp(one, two, t):
    """Spherical linear interpolation.

    :param one: Start quaternion.
    :param two: End quaternion.
    :param t: Parameter value(s) between 0 and 1.

    """
    return (two * one.inverse())**t * one


def canonicalized(quaternions):
    """Iterator adapter to ensure minimal angles between quaternions."""
    p = UnitQuaternion.from_unit_xyzw((0, 0, 0, 1))
    for q in quaternions:
        if p.dot(q) < 0:
            q = -q
        yield q
        p = q


class PiecewiseSlerp:
    """Piecewise Slerp, see __init__()."""

    def __init__(self, rotations, *, grid=None, closed=False):
        """Piecewise Slerp.

        :param rotations: Sequence of quaternions to be interpolated.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            Must have the same length as *rotations*, except when *closed*
            is ``True``, where it must be one element longer.
        :type grid: optional
        :param closed: If ``True``, the first quaternion is repeated at
            the end.
        :type closed: optional

        """
        rotations = list(rotations)
        if len(rotations) < 2:
            raise ValueError('At least two rotations are required')
        if closed:
            rotations += rotations[:1]
        rotations = list(canonicalized(rotations))
        if grid is None:
            grid = range(len(rotations))
        if len(rotations) != len(grid):
            raise ValueError(
                'Number of grid values must be same as '
                'rotations (one more for closed curves)')
        self.rotations = rotations
        self.grid = list(grid)

    def evaluate(self, t, n=0):
        if n != 0:
            raise NotImplementedError('Derivatives are not implemented yet')
        if not _np.isscalar(t):
            return _np.array([self.evaluate(t) for t in t])
        idx = _check_param('t', t, self.grid)
        t0, t1 = self.grid[idx:idx + 2]
        return slerp(
            self.rotations[idx],
            self.rotations[idx + 1],
            (t - t0) / (t1 - t0))


class DeCasteljau:
    """De Casteljau's algorithm, see __init__()."""

    def __init__(self, segments, grid=None):
        """De Casteljau's algorithm using `slerp()`.

        See `the corresponding notebook`__ for details.

        __ ../rotation/de-casteljau.ipynb

        :param segments: Sequence of segments,
            each one consisting of multiple control quaternions.
            Different segments can have different numbers of control points.
        :param grid: Sequence of parameter values corresponding to
            segment boundaries.  Must be strictly increasing.
        :type grid: optional

        """
        if len(segments) < 1:
            raise ValueError('There must be at least one segment')
        if grid is None:
            # uniform parameterization
            grid = range(len(segments) + 1)
        self.segments = segments
        self.grid = list(grid)

    def evaluate(self, t, n=0):
        """Get value or angular velocity at given parameter value(s).

        :param t: Parameter value(s).
        :param n: Use ``0`` for calculating the value (a quaternion),
            ``1`` for the angular velocity (a 3-element vector).
        :type n: {0, 1}, optional

        """
        if not _np.isscalar(t):
            return _np.array([self.evaluate(t, n) for t in t])
        t, segment = self._select_segment_and_normalize_t(t)
        if n == 0:
            return slerp(*_reduce(segment, t), t)
        elif n == 1:
            one, two = _reduce(segment, t)
            x, y, z = (two * one.inverse()).log_map()
            degree = len(segment) - 1
            # NB: twice the angle
            return x * 2 * degree, y * 2 * degree, z * 2 * degree
        else:
            raise ValueError('Unsupported n: {!r}'.format(n))

    def _select_segment_and_normalize_t(self, t):
        idx = _check_param('t', t, self.grid)
        t0, t1 = self.grid[idx:idx + 2]
        t = (t - t0) / (t1 - t0)
        return t, self.segments[idx]


def _reduce(segment, t):
    """Obtain two quaternions for the last step of the algorithm.

    Recursively applies `slerp()` to neighboring control quaternions
    until only two are left.

    """
    if len(segment) < 2:
        raise ValueError('Segment must have at least two quaternions')
    if len(segment) == 2:
        return segment
    return _reduce([
        slerp(one, two, t)
        for one, two in zip(segment, segment[1:])
    ], t)
