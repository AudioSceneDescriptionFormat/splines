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

    def __add__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        s_x, s_y, s_z = self.vector
        o_x, o_y, o_z = other.vector
        return Quaternion(
            scalar=self.scalar + other.scalar,
            vector=(s_x + o_x, s_y + o_y, s_z + o_z))

    def __sub__(self, other):
        if not isinstance(other, Quaternion):
            return NotImplemented
        s_x, s_y, s_z = self.vector
        o_x, o_y, o_z = other.vector
        return Quaternion(
            scalar=self.scalar - other.scalar,
            vector=(s_x - o_x, s_y - o_y, s_z - o_z))

    def __mul__(self, other):
        if isinstance(self, UnitQuaternion) and isinstance(other, UnitQuaternion):
            result_type = UnitQuaternion
        elif isinstance(other, Quaternion):
            result_type = Quaternion
        else:
            x, y, z = self.vector
            return Quaternion(
                scalar=self.scalar * other,
                vector=(x * other, y * other, z * other))
        a1 = self.scalar
        b1, c1, d1 = self.vector
        a2 = other.scalar
        b2, c2, d2 = other.vector
        return Quaternion.__new__(
            result_type,
            a1*a2 - b1*b2 - c1*c2 - d1*d2,
            (
                a1*b2 + b1*a2 + c1*d2 - d1*c2,
                a1*c2 - b1*d2 + c1*a2 + d1*b2,
                a1*d2 + b1*c2 - c1*b2 + d1*a2,
            )
        )

    def __rmul__(self, other):
        if not isinstance(other, Quaternion):
            return self.__mul__(other)
        assert False

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

    def __new__(cls):
        return super().__new__(cls, 1.0, (0.0, 0.0, 0.0))

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
        if self.scalar >= 1:
            return super().__new__(UnitQuaternion, self.scalar, self.vector)
        return UnitQuaternion.exp_map(alpha * self.log_map())

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

        """
        if self.scalar >= 1:
            return _np.zeros_like(self.vector)
        length = self.angle / 2
        return self.axis * length

    @property
    def axis(self):
        assert self.scalar <= 1
        # NB: This is the same as sqrt(x**2 + y**2 + z**2)
        norm = _math.sqrt(1 - self.scalar**2)
        return _np.true_divide(self.vector, norm)

    @property
    def angle(self):
        return 2 * _math.acos(self.scalar)

    def rotation_to(self, other):
        """Rotation required to rotate *self* into *other*.

        See :ref:`/rotation/quaternions.ipynb#Relative-Rotation-(Global-Frame-of-Reference)`.

        """
        return other * self.inverse()

    def rotate_vector(self, v):
        rotated = self * Quaternion(0, v) * self.inverse()
        return rotated.vector


def slerp(one, two, t):
    """Spherical linear interpolation.

    See :ref:`/rotation/slerp.ipynb`.

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

    def __init__(self, quaternions, *, grid=None, closed=False):
        """Piecewise Slerp.

        See :ref:`/rotation/slerp.ipynb#Piecewise-Slerp`.

        :param quaternions: Sequence of rotations to be interpolated.
            The quaternions will be `canonicalized()`.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            Must have the same length as *quaternions*, except when *closed*
            is ``True``, where it must be one element longer.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param closed: If ``True``, the first quaternion is repeated at
            the end.
        :type closed: optional

        """
        quaternions = _check_quaternions(quaternions, closed=closed)
        if grid is None:
            grid = range(len(quaternions))
        if len(quaternions) != len(grid):
            raise ValueError(
                'Number of grid values must be same as '
                'quaternions (one more for closed curves)')
        self.quaternions = quaternions
        self.grid = list(grid)

    def evaluate(self, t, n=0):
        if n != 0:
            raise NotImplementedError('Derivatives are not implemented yet')
        if not _np.isscalar(t):
            return _np.array([self.evaluate(t) for t in t])
        idx = _check_param('t', t, self.grid)
        t0, t1 = self.grid[idx:idx + 2]
        return slerp(
            self.quaternions[idx],
            self.quaternions[idx + 1],
            (t - t0) / (t1 - t0))


class DeCasteljau:
    """Spline using De Casteljau's algorithm, see __init__()."""

    def __init__(self, segments, grid=None):
        """Spline using De Casteljau's algorithm with `slerp()`.

        See `the corresponding notebook`__ for details.

        __ ../rotation/de-casteljau.ipynb

        :param segments: Sequence of segments,
            each one consisting of multiple control quaternions.
            Different segments can have different numbers of control points.
        :param grid: Sequence of parameter values corresponding to
            segment boundaries.  Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional

        """
        if len(segments) < 1:
            raise ValueError('There must be at least one segment')
        if grid is None:
            # uniform parameterization
            grid = range(len(segments) + 1)
        else:
            if len(segments) + 1 != len(grid):
                raise ValueError(
                    'Number of grid values must be one more '
                    'than number of segments')
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
        segment, t, delta_t = self._select_segment_and_normalize_t(t)
        if n == 0:
            return slerp(*_reduce_de_casteljau(segment, t), t)
        elif n == 1:
            one, two = _reduce_de_casteljau(segment, t)
            tangent = (two * one.inverse()).log_map()
            degree = len(segment) - 1
            # NB: twice the angle
            return tangent * 2 * degree / delta_t
        else:
            raise ValueError('Unsupported n: {!r}'.format(n))

    def _select_segment_and_normalize_t(self, t):
        idx = _check_param('t', t, self.grid)
        t0, t1 = self.grid[idx:idx + 2]
        delta_t = t1 - t0
        t = (t - t0) / delta_t
        return self.segments[idx], t, delta_t


def _reduce_de_casteljau(segment, t):
    """Obtain two quaternions for the last step of De Castelau's algorithm.

    Repeatedly applies `slerp()` to neighboring control quaternions
    until only two are left.

    """
    if len(segment) < 2:
        raise ValueError('Segment must have at least two quaternions')
    while len(segment) > 2:
        segment = [
            slerp(one, two, t)
            for one, two in zip(segment, segment[1:])]
    return segment


class KochanekBartels(DeCasteljau):
    """Kochanek--Bartels-like rotation spline, see __init__()."""

    @staticmethod
    def _calculate_control_quaternions(quaternions, times, tcb):
        q_1, q0, q1 = quaternions
        t_1, t0, t1 = times
        T, C, B = tcb
        a = (1 - T) * (1 + C) * (1 + B)
        b = (1 - T) * (1 - C) * (1 - B)
        c = (1 - T) * (1 - C) * (1 + B)
        d = (1 - T) * (1 + C) * (1 - B)

        q_in = q0 * q_1.inverse()
        q_out = q1 * q0.inverse()
        v_in = q_in.log_map() / (t0 - t_1)
        v_out = q_out.log_map() / (t1 - t0)

        def v0(weight_in, weight_out):
            return (
                weight_in * (t1 - t0) * v_in +
                weight_out * (t0 - t_1) * v_out
            ) / (t1 - t_1)

        return [
            UnitQuaternion.exp_map(-v0(c, d) * (t0 - t_1) / 3) * q0,
            UnitQuaternion.exp_map(v0(a, b) * (t1 - t0) / 3) * q0,
        ]

    def __init__(self, quaternions, grid=None, *, tcb=(0, 0, 0), alpha=None,
                 endconditions='natural'):
        """Kochanek--Bartels-like rotation spline.

        :param quaternions: Sequence of rotations to be interpolated.
            The quaternions will be `canonicalized()`.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param tcb: Sequence of *tension*, *continuity* and *bias* triples.
            TCB values can only be given for the interior quaternions.
            If only two quaternions are given, TCB values are ignored.
        :type tcb: optional
        :param alpha: TODO
        :type alpha: optional
        :param endconditions: Start/end conditions. Can be ``'closed'``,
            ``'natural'`` or pair of tangent vectors (a.k.a. "clamped").

            TODO: clamped

            If ``'closed'``, the first rotation is re-used as last rotation
            and an additional *grid* time has to be specified.
        :type endconditions: optional

        """
        closed = endconditions == 'closed'
        if closed:
            tcb_slots = len(quaternions)
        else:
            tcb_slots = len(quaternions) - 2
        quaternions = _check_quaternions(quaternions, closed=closed)
        grid = _check_grid(grid, alpha, quaternions)
        tcb = _np.asarray(tcb)
        if tcb.ndim == 1 and len(tcb) == 3:
            tcb = _np.tile(tcb, (tcb_slots, 1))
        elif len(tcb) != tcb_slots:
            raise ValueError(
                'There must be two more quaternions than TCB values '
                '(except for closed curves)')
        if closed:
            tcb = _np.row_stack([tcb, tcb[0]])
        start, end, zip_quaternions, zip_grid = _check_endconditions(
            endconditions, quaternions, grid)

        control_points = []
        for qs, ts, tcb in zip(zip_quaternions, zip_grid, tcb):
            q_before, q_after = self._calculate_control_quaternions(
                qs, ts, tcb)
            control_points.extend([q_before, qs[1], qs[1], q_after])
        if closed:
            assert len(grid) * 4 == len(control_points)
            control_points = control_points[2:-2]
        elif not control_points:
            # two quaternions -> spherical linear interpolation
            assert len(quaternions) == 2
            assert len(grid) == 2
            assert not closed
            assert not tcb
            q0, q1 = quaternions
            offset = (q1 * q0.inverse())**(1 / 3)  # "cubic" spline, degree 3
            control_points = [q0, offset * q0, offset.inverse() * q1, q1]
        else:
            control_points.insert(0, _end_control_quaternion(
                start,
                [quaternions[0], control_points[0], control_points[1]]))
            control_points.insert(0, quaternions[0])
            control_points.append(_end_control_quaternion(
                end,
                [quaternions[-1], control_points[-1], control_points[-2]]))
            control_points.append(quaternions[-1])
        segments = list(zip(*[iter(control_points)] * 4))
        DeCasteljau.__init__(self, segments, grid)


class CatmullRom(KochanekBartels):
    """Catmull--Rom-like rotation spline, see __init__()."""

    def __init__(self, quaternions, grid=None, *, alpha=None,
                 endconditions='natural'):
        """Catmull--Rom-like rotation spline.

        This is just `KochanekBartels` without TCB values.

        """
        super().__init__(
            quaternions, grid, tcb=(0, 0, 0),
            alpha=alpha, endconditions=endconditions)


def _check_quaternions(quaternions, *, closed):
    """Canonicalize and if closed, append first rotation at the end."""
    quaternions = list(quaternions)
    if len(quaternions) < 2:
        raise ValueError('At least two quaternions are required')
    if closed:
        quaternions = quaternions + quaternions[:1]
    return list(canonicalized(quaternions))


def _check_grid(grid, alpha, quaternions):
    if grid is None:
        if alpha is None:
            # NB: This is the same as alpha=0, except the type is int
            return range(len(quaternions))
        grid = [0]
        for a, b in zip(quaternions, quaternions[1:]):
            delta = (2 * _math.acos(a.dot(b)))**alpha
            if delta == 0:
                raise ValueError(
                    'Repeated quaternions are not possible with alpha != 0')
            grid.append(grid[-1] + delta)
    else:
        if alpha is not None:
            raise TypeError('Only one of {grid, alpha} is allowed')
        if len(quaternions) != len(grid):
            raise ValueError('Number of grid values must be same as '
                             'quaternions (one more for closed curves)')
        # TODO: check if grid values are increasing?
    return grid


def _check_endconditions(endconditions, quaternions, grid):
    if endconditions == 'closed':
        # NB: the first quaternion has already been appended to the end in
        #     _check_quaternions()
        prefix = quaternions[-2]
        if prefix.dot(quaternions[0]) < 0:
            prefix = -prefix
        suffix = quaternions[1]
        if quaternions[-1].dot(suffix) < 0:
            suffix = -suffix
        quaternions = [prefix] + quaternions + [suffix]
        grid = [
            grid[0] - (grid[-1] - grid[-2]),
            *grid,
            grid[-1] + (grid[1] - grid[0]),
        ]
        start = end = None
    elif isinstance(endconditions, str):
        start = end = endconditions
    else:
        try:
            start, end = endconditions
        except (TypeError, ValueError):
            raise TypeError('endconditions must be a string or a pair')
    assert len(quaternions) == len(grid)
    triples = [zip(arg, arg[1:], arg[2:]) for arg in (quaternions, grid)]
    return (start, end, *triples)


def _end_control_quaternion(condition, quaternions):
    if condition == 'natural':
        return _natural_control_quaternion(quaternions)
    elif _np.shape(condition) == _np.shape(quaternions[0]):
        #tangent = condition
        raise NotImplementedError('TODO')
    raise ValueError(
        f'{condition!r} is not a valid start/end condition')


def _natural_control_quaternion(quaternions):
    """Return second control quaternion given the other three."""
    outer, inner_control, inner = quaternions
    return (
        (inner_control * inner.inverse()) *
        (inner * outer.inverse())
    )**(1 / 2) * outer


class BarryGoldman:
    """Rotation spline using Barry--Goldman algorithm, see __init__()."""

    def __init__(self, quaternions, grid=None, *, alpha=None):
        """Rotation spline using Barry--Goldman algorithm.

        Always closed (for now).

        """
        # TODO: what happens when exactly 2 quaternions are given?
        self.quaternions = _check_quaternions(quaternions, closed=True)
        self.grid = list(_check_grid(grid, alpha, self.quaternions))

    def evaluate(self, t):
        if not _np.isscalar(t):
            return _np.array([self.evaluate(t) for t in t])
        idx = _check_param('t', t, self.grid)
        q0, q1 = self.quaternions[idx:idx + 2]
        t0, t1 = self.grid[idx:idx + 2]
        if idx == 0:
            q_1 = self.quaternions[-2]
            if q_1.dot(q0) < 0:
                q_1 = -q_1
            t_1 = t0 - (self.grid[-1] - self.grid[-2])
        else:
            q_1 = self.quaternions[idx - 1]
            t_1 = self.grid[idx - 1]
        if idx + 2 == len(self.quaternions):
            q2 = self.quaternions[1]
            if q1.dot(q2) < 0:
                q2 = -q2
            assert len(self.quaternions) == len(self.grid)
            t2 = t1 + (self.grid[1] - self.grid[0])
        else:
            q2 = self.quaternions[idx + 2]
            t2 = self.grid[idx + 2]
        slerp_0_1 = slerp(q0, q1, (t - t0) / (t1 - t0))
        return slerp(
            slerp(
                slerp(q_1, q0, (t - t_1) / (t0 - t_1)),
                slerp_0_1,
                (t - t_1) / (t1 - t_1)),
            slerp(
                slerp_0_1,
                slerp(q1, q2, (t - t1) / (t2 - t1)),
                (t - t0) / (t2 - t0)),
            (t - t0) / (t1 - t0))
