"""Piecewise polynomial curves (in Euclidean space).

.. rubric:: Submodules

.. autosummary::

    quaternion

"""
from bisect import bisect_right as _bisect_right, bisect_left as _bisect_left
from itertools import accumulate as _accumulate
from math import factorial as _factorial

import numpy as _np


__version__ = '0.1.0'


class Monomial:
    """Piecewise polynomial curve, see __init__()."""

    def __init__(self, segments, grid):
        r"""Piecewise polynomial curve using monomial basis.

        Arbitrary degree, arbitrary dimension.

        .. math::

            \boldsymbol{p}_i(t) = \sum_{k=0}^n
                \boldsymbol{a}_k \left(\frac{t - t_i}{t_{i+1} - t_i}\right)^k
                \text{ for } t_i \leq t < t_{i+1}

        Similar to https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PPoly.html,
        which states:

            "High-order polynomials in the power basis can be numerically
            unstable.  Precision problems can start to appear for orders
            larger than 20-30."


        :param segments: Sequence of polynomial segments.
            Each segment contains coefficients for the monomial basis
            (in order of decreasing degree).
            Different segments can have different polynomial degree.
        :param grid: Sequence of parameter values corresponding to
            segment boundaries.  Must be strictly increasing.

        """
        self.segments = [_np.array(coefficients, copy=True)
                         for coefficients in segments]
        if grid is None:
            grid = range(len(segments) + 1)
        self.grid = list(grid)

    def evaluate(self, t, n=0):
        """Get value (or *n*-th derivative) at given parameter value(s)."""
        if not _np.isscalar(t):
            return _np.array([self.evaluate(time, n) for time in t])

        idx = _check_param('t', t, self.grid)

        t0, t1 = self.grid[idx:idx + 2]
        t = (t - t0) / (t1 - t0)
        coefficients = self.segments[idx][:-n or None]
        powers = _np.arange(len(coefficients))[::-1]
        product = _np.multiply.reduce
        weights = product([powers + 1 + i for i in range(n)]) / (t1 - t0)**n
        return t**powers * weights @ coefficients


class Bernstein:
    """Piecewise Bézier curve, see __init__()."""

    @staticmethod
    def basis(degree, t):
        r"""Bernstein basis polynomials of given *degree*, evaluated at *t*.

        Returns a list of values corresponding to
        :math:`i = 0, \ldots, n`, given the degree :math:`n`,
        using the formula

        .. math::

            b_{i,n}(t) = {n \choose i} t^i \left( 1 - t \right)^{n - i},

        with the *binomial coefficient*
        :math:`{n \choose i} = \frac{n!}{i!(n - i)!}`.

        """
        return [
            _comb(degree, i) * t**i * (1 - t)**(degree - i)
            for i in range(degree + 1)]

    def __init__(self, segments, grid=None):
        """Piecewise Bézier curve using Bernstein basis.

        :param segments: Sequence of segments,
            each one consisting of multiple Bézier control points.
            Different segments can have different numbers of control points
            (and therefore different polynomial degrees).
        :param grid: Sequence of parameter values corresponding to
            segment boundaries.  Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional

        """
        self.segments = [_np.array(control_points, copy=True)
                         for control_points in segments]
        if grid is None:
            grid = range(len(segments) + 1)
        self.grid = list(grid)

    def evaluate(self, t, n=0):
        """Get value at the given parameter value(s)."""
        if n != 0:
            raise NotImplementedError('Derivatives are not implemented yet')
        if not _np.isscalar(t):
            return _np.array([self.evaluate(time, n) for time in t])
        idx = _check_param('t', t, self.grid)
        t0, t1 = self.grid[idx:idx + 2]
        t = (t - t0) / (t1 - t0)
        control_points = self.segments[idx]
        degree = len(control_points) - 1
        return sum(
            a * b
            for a, b in zip(control_points, self.basis(degree, t)))


def _check_param(name, param, grid):
    if param < grid[0]:
        raise ValueError(f'{name} too small: {param}')
    elif param < grid[-1]:
        idx = _bisect_right(grid, param) - 1
    elif param == grid[-1]:
        idx = len(grid) - 2
    else:
        raise ValueError(f'{name} too big: {param}')
    return idx


def _comb(n, k):
    # NB: Python 3.8 has math.comb()
    return _factorial(n) // _factorial(k) // _factorial(n - k)


def _check_vertices(vertices, *, closed):
    """For closed curves, append first vertex at the end."""
    if len(vertices) < 2:
        raise ValueError('At least two vertices are required')
    if closed:
        vertices = _np.concatenate([vertices, vertices[:1]])
    return vertices


def _check_grid(grid, alpha, vertices):
    if grid is None:
        if alpha is None:
            # NB: This is the same as alpha=0, except the type is int
            return range(len(vertices))
        vertices = _np.asarray(vertices)
        grid = [0]
        for x0, x1 in zip(vertices, vertices[1:]):
            delta = _np.linalg.norm(x1 - x0)**alpha
            if delta == 0:
                raise ValueError(
                    'Repeated vertices are not possible with alpha != 0')
            grid.append(grid[-1] + delta)
    else:
        if alpha is not None:
            raise TypeError('Only one of {grid, alpha} is allowed')
        if len(vertices) != len(grid):
            raise ValueError('Number of grid values must be same as '
                             'vertices (one more for closed curves)')
        # TODO: check if grid values are increasing?
    return grid


def _check_endconditions(endconditions, vertices, grid):
    if endconditions == 'closed':
        second_vertex = vertices[1:2]
        vertices = _np.concatenate([vertices, second_vertex])
        first_interval = grid[1] - grid[0]
        grid = list(grid) + [grid[-1] + first_interval]
        start = end = None
    elif isinstance(endconditions, str):
        start = end = endconditions
    else:
        try:
            start, end = endconditions
        except (TypeError, ValueError):
            raise TypeError('endconditions must be a string or a pair')
    triples = [zip(arg, arg[1:], arg[2:]) for arg in (vertices, grid)]
    return (start, end, *triples)


def _check_tangents(tangents, vertices, grid, start, end, *, closed):
    if closed:
        # Move last (outgoing) tangent to the beginning:
        tangents = tangents[-1:] + tangents[:-1]
    elif not tangents:
        # straight line
        assert len(vertices) == 2
        assert len(grid) == 2
        vertices = _np.asarray(vertices)
        tangents = [(vertices[1] - vertices[0]) / (grid[1] - grid[0])] * 2
    else:
        tangents.insert(0, _end_tangent(
            start, vertices[:2], grid[:2], tangents[0]))
        tangents.append(_end_tangent(
            end, vertices[-2:], grid[-2:], tangents[-1]))
    return tangents


def _end_tangent(condition, vertices, times, other_tangent):
    if condition == 'natural':
        tangent = _natural_tangent(vertices, times, other_tangent)
    elif _np.shape(condition) == _np.shape(vertices[0]):
        tangent = condition
    else:
        raise ValueError(
            f'{condition!r} is not a valid start/end condition')
    return tangent


def _natural_tangent(vertices, times, tangent):
    """Calculate tangent for "natural" end condition.

    Given 2 points and one tangent, this returns the tangent for the
    other side that results from the second derivative being zero.

    See :ref:`end-conditions-natural.ipynb`.

    """
    x0, x1 = _np.asarray(vertices)
    t0, t1 = times
    delta = t1 - t0
    return (3 * x1 - 3 * x0 - delta * tangent) / (2 * delta)


class CubicHermite(Monomial):
    """Cubic Hermite curve, see __init__()."""

    matrix = _np.array([
        [2, -2, 1, 1],
        [-3, 3, -2, -1],
        [0, 0, 1, 0],
        [1, 0, 0, 0]])

    def __init__(self, vertices, tangents, grid=None):
        """Cubic Hermite curve.

        See :ref:`/euclidean/hermite.ipynb`.

        :param vertices: Sequence of vertices.
        :param tangents: Sequence of tangent vectors
            (two per segment, outgoing and incoming).
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional

        """
        if len(vertices) < 2:
            raise ValueError('At least 2 vertices are needed')
        if len(tangents) != 2 * (len(vertices) - 1):
            raise ValueError('Exactly 2 tangents per segment are needed')
        if grid is None:
            grid = range(len(vertices))
        if len(vertices) != len(grid):
            raise ValueError('As many grid times as vertices are needed')
        tangents = _np.asarray(tangents)
        segments = [
            self.matrix @ [x0, x1, (t1 - t0) * v0, (t1 - t0) * v1]
            for (x0, x1), (v0, v1), (t0, t1) in zip(
                zip(vertices, vertices[1:]),
                zip(tangents[::2], tangents[1::2]),
                zip(grid, grid[1:]))]
        Monomial.__init__(self, segments, grid)


class CatmullRom(CubicHermite):
    """Catmull--Rom spline, see __init__()."""

    # NB: Catmull-Rom could be implemented as special case of Kochanek-Bartels,
    #     but here we chose not to.

    # NB: We could use the basis matrix for Catmull-Rom splines, but
    #     this wouldn't work if only 3 vertices are given by the user.
    #     Since we have to handle this special case anyway, we use the same
    #     method for everything.  Apart from reducing the amount of code, this
    #     also allows us to define derived classes that overwrite
    #     _calculate_tangent().

    @staticmethod
    def _calculate_tangent(points, times):
        x_1, x0, x1 = _np.asarray(points)
        t_1, t0, t1 = times
        delta_1 = t0 - t_1
        delta0 = t1 - t0
        return ((delta0**2 * (x0 - x_1) + delta_1**2 * (x1 - x0)) /
                (delta0 * delta_1 * (delta0 + delta_1)))

    def __init__(self, vertices, grid=None, *, alpha=None,
                 endconditions='natural'):
        """Catmull--Rom spline.

        This class implements one specific member of the family of
        splines described in :cite:`catmull1974splines`,
        which is commonly known as *Catmull--Rom spline*:
        The cubic spline that can be constructed by linear Lagrange
        interpolation (and extrapolation) followed by quadratic B-spline
        blending, or equivalently, quadratic Lagrange interpolation
        followed by linear B-spline blending.

        The implementation used in this class, however, does nothing of
        that sort.  It simply calculates the appropriate tangent vectors
        at the control points and instantiates a `CubicHermite` spline.

        See :ref:`/euclidean/catmull-rom.ipynb`.

        :param vertices: Sequence of vertices.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param alpha: TODO
        :type alpha: optional
        :param endconditions: Start/end conditions. Can be ``'closed'``,
            ``'natural'`` or pair of tangent vectors (a.k.a. "clamped").
            If ``'closed'``, the first vertex is re-used as last vertex
            and an additional *grid* time has to be specified.
        :type endconditions: optional

        """
        closed = endconditions == 'closed'
        vertices = _check_vertices(vertices, closed=closed)
        grid = _check_grid(grid, alpha, vertices)
        start, end, zip_vertices, zip_grid = _check_endconditions(
            endconditions, vertices, grid)
        tangents = [
            self._calculate_tangent(points, times)
            for points, times in zip(zip_vertices, zip_grid)]
        # Duplicate tangents (incoming and outgoing are the same):
        tangents = [x for tangent in tangents for x in (tangent, tangent)]
        tangents = _check_tangents(
            tangents, vertices, grid, start, end, closed=closed)
        CubicHermite.__init__(self, vertices, tangents, grid)


class KochanekBartels(CubicHermite):
    """Kochanek--Bartels spline, see __init__()."""

    @staticmethod
    def _calculate_tangents(points, times, tcb):
        x_1, x0, x1 = _np.asarray(points)
        t_1, t0, t1 = times
        T, C, B = tcb
        a = (1 - T) * (1 + C) * (1 + B)
        b = (1 - T) * (1 - C) * (1 - B)
        c = (1 - T) * (1 - C) * (1 + B)
        d = (1 - T) * (1 + C) * (1 - B)
        incoming = (
            c * (t1 - t0)**2 * (x0 - x_1) + d * (t0 - t_1)**2 * (x1 - x0)
        ) / (
            (t1 - t0) * (t0 - t_1) * (t1 - t_1)
        )
        outgoing = (
            a * (t1 - t0)**2 * (x0 - x_1) + b * (t0 - t_1)**2 * (x1 - x0)
        ) / (
            (t1 - t0) * (t0 - t_1) * (t1 - t_1)
        )
        return incoming, outgoing

    def __init__(self, vertices, grid=None, *, tcb=(0, 0, 0), alpha=None,
                 endconditions='natural'):
        """Kochanek--Bartels spline.

        See :ref:`/euclidean/kochanek-bartels.ipynb`.

        :param vertices: Sequence of vertices.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param tcb: Sequence of *tension*, *continuity* and *bias* triples.
            TCB values can only be given for the interior vertices.
        :type tcb: optional
        :param alpha: TODO
        :type alpha: optional
        :param endconditions: Start/end conditions. Can be ``'closed'``,
            ``'natural'`` or pair of tangent vectors (a.k.a. "clamped").
            If ``'closed'``, the first vertex is re-used as last vertex
            and an additional *grid* time has to be specified.
        :type endconditions: optional

        """
        closed = endconditions == 'closed'
        if closed:
            tcb_slots = len(vertices)
        else:
            tcb_slots = len(vertices) - 2
        vertices = _check_vertices(vertices, closed=closed)
        grid = _check_grid(grid, alpha, vertices)
        tcb = _np.asarray(tcb)
        if tcb.ndim == 1 and len(tcb) == 3:
            tcb = _np.tile(tcb, (tcb_slots, 1))
        elif len(tcb) != tcb_slots:
            raise ValueError('There must be two more vertices than TCB values '
                             '(except for closed curves)')
        if closed:
            # Move first TCB value to the end:
            tcb = _np.roll(tcb, -1, axis=0)
        start, end, zip_vertices, zip_grid = _check_endconditions(
            endconditions, vertices, grid)
        tangents = [
            tangent
            for points, times, tcb in zip(zip_vertices, zip_grid, tcb)
            for tangent in self._calculate_tangents(points, times, tcb)]
        tangents = _check_tangents(
            tangents, vertices, grid, start, end, closed=closed)
        CubicHermite.__init__(self, vertices, tangents, grid)


class Natural(CubicHermite):
    """Natural spline, see __init__()."""

    def __init__(self, vertices, grid=None, *, alpha=None,
                 endconditions='natural'):
        """Natural spline.

        See :ref:`/euclidean/natural.ipynb`.

        :param vertices: Sequence of vertices.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param alpha: TODO
        :type alpha: optional
        :param endconditions: Start/end conditions. Can be ``'closed'``,
            ``'natural'`` or pair of tangent vectors (a.k.a. "clamped").
            If ``'closed'``, the first vertex is re-used as last vertex
            and an additional *grid* time has to be specified.
        :type endconditions: optional

        """
        N = len(vertices)
        A = _np.zeros((N, N))
        b = _np.zeros_like(vertices)
        closed = endconditions == 'closed'
        vertices = _check_vertices(vertices, closed=closed)
        vertices = _np.asarray(vertices)
        grid = _check_grid(grid, alpha, vertices)
        delta = _np.diff(grid)

        for i in range(0, N) if closed else range(1, N - 1):
            A[i, i - 1] = 1 / delta[i - 1]
            A[i, i] = (2 / delta[i - 1] + 2 / delta[i])
            A[i, (i + 1) % N] = 1 / delta[i]
            b[i] = 3 * (
                (vertices[i] - vertices[(i - 1) % N]) / delta[i - 1]**2 +
                (vertices[i + 1] - vertices[i]) / delta[i]**2)

        if closed:
            # Nothing to do here, the first and last row have already
            # been populated in the for-loop above.
            pass
        else:
            if isinstance(endconditions, str):
                start = end = endconditions
            else:
                try:
                    start, end = endconditions
                except (TypeError, ValueError):
                    raise TypeError('endconditions must be a string or a pair')
            if start == 'natural':
                A[0, 0:2] = 2 * delta[0], delta[1]
                b[0] = 3 * (vertices[1] - vertices[0])
            elif _np.shape(start) == _np.shape(vertices[0]):
                A[0, 0] = 1
                b[0] = start
            else:
                raise ValueError(f'{start!r} is not a valid start condition')
            if end == 'natural':
                A[N - 1, N - 2:] = delta[N - 2], 2 * delta[N - 2]
                b[-1] = 3 * (vertices[N - 1] - vertices[N - 2])
            elif _np.shape(end) == _np.shape(vertices[0]):
                A[N - 1, N - 1] = 1
                b[N - 1] = end
            else:
                raise ValueError(f'{end!r} is not a valid end condition')

        tangents = _np.linalg.solve(A, b)

        if closed:
            tangents = _np.concatenate([tangents, tangents[:1]])
        # Duplicate inner tangents (incoming and outgoing are the same):
        tangents = _np.concatenate(
            [tangents[i:i + 2] for i in range(len(tangents) - 1)])
        CubicHermite.__init__(self, vertices, tangents, grid)


def _monotone_end_condition(inner_slope, secant_slope):
    """

    Return the "outer" (i.e. first or last) slope given the "inner"
    (i.e. second or penultimate) slope and the secant slope of the
    first or last segment.

    """
    # NB: This is a very ad-hoc algorithm meant to minimize the change in slope
    # within the first/last curve segment.  Especially, this should avoid a
    # change from negative to positive acceleration (and vice versa).
    # There might be a better method available!?!
    if secant_slope < 0:
        return -_monotone_end_condition(-inner_slope, -secant_slope)
    assert 0 <= inner_slope <= 3 * secant_slope
    if inner_slope <= secant_slope:
        return 3 * secant_slope - 2 * inner_slope
    else:
        return (3 * secant_slope - inner_slope) / 2


class PiecewiseMonotoneCubic(CatmullRom):
    """Piecewise monotone cubic curve, see __init__()."""

    def __init__(self, values, grid=None, slopes=None, *, alpha=None,
                 closed=False):
        """Piecewise monotone cubic curve.

        See :ref:`/euclidean/piecewise-monotone.ipynb`.

        This only works for one-dimensional values.

        For undefined slopes, ``_calculate_tangent()`` is called on
        the base class.

        :param values: Sequence of values to be interpolated.
        :param grid: Sequence of parameter values.
            Must be strictly increasing.
            If not specified, a uniform grid is used (0, 1, 2, 3, ...).
        :type grid: optional
        :param slopes: Sequence of slopes or ``None`` if slope should be
            computed from neighboring values.  An error is raised
            if a segment would become non-monotone with a given slope.
        :type slopes: optional

        """
        if len(values) < 2:
            raise ValueError('At least two values are required')
        if closed:
            values = _np.concatenate([values, values[:1]])
        grid = _check_grid(grid, alpha, values)
        if slopes is None:
            slopes = (None,) * len(values)
        elif closed:
            slopes = *slopes, slopes[0]
        if len(values) != len(slopes):
            raise ValueError('Same number of values and slopes is required')

        # TODO: check strictly increasing times?

        if closed:
            second_value = values[1:2]
            values = _np.concatenate([values, second_value])
            first_interval = grid[1] - grid[0]
            grid = list(grid) + [grid[-1] + first_interval]

        def fix_slope(slope, left, right):
            """Manipulate the slope to preserve monotonicity.

            See Dougherty et al. (1989), eq. (4.2)

            """
            if left * right <= 0:
                return 0
            elif right > 0:
                return min(max(0, slope), 3 * min(abs(left), abs(right)))
            else:
                return max(min(0, slope), -3 * min(abs(left), abs(right)))

        final_slopes = []
        for xs, ts, slope in zip(zip(values, values[1:], values[2:]),
                                 zip(grid, grid[1:], grid[2:]),
                                 slopes[1:]):
            x_1, x0, x1 = xs
            t_1, t0, t1 = ts
            left = (x0 - x_1) / (t0 - t_1)
            right = (x1 - x0) / (t1 - t0)
            if slope is None:
                # NB: This has to be defined on the parent class:
                slope = self._calculate_tangent(xs, ts)
                slope = fix_slope(slope, left, right)
            else:
                if slope != fix_slope(slope, left, right):
                    raise ValueError(f'Slope too steep: {slope}')
            final_slopes.append(slope)  # incoming
            final_slopes.append(slope)  # outgoing

        if closed:
            # Move last outgoing slope to front:
            final_slopes = final_slopes[-1:] + final_slopes[:-1]
            values = values[:-1]
            grid = grid[:-1]
        elif not final_slopes:
            secant_slope = (values[1] - values[0]) / (grid[1] - grid[0])
            one, two = slopes

            def check_slope(slope):
                if slope != fix_slope(slope, secant_slope, secant_slope):
                    raise ValueError(f'Slope too steep or wrong sign: {slope}')

            if one is None:
                if two is None:
                    final_slopes = [secant_slope] * 2
                else:
                    check_slope(two)
                    final_slopes = [
                        _monotone_end_condition(two, secant_slope),
                        two]
            else:
                if two is None:
                    check_slope(one)
                    final_slopes = [
                        one,
                        _monotone_end_condition(one, secant_slope)]
                else:
                    check_slope(one)
                    check_slope(two)
                    final_slopes = [one, two]
        else:

            def end_slope(outer, inner, secant_slope):
                if outer is None:
                    outer = _monotone_end_condition(inner, secant_slope)
                else:
                    if outer != fix_slope(outer, secant_slope, secant_slope):
                        raise ValueError(
                            f'Slope too steep or wrong sign: {outer}')
                return outer

            final_slopes.insert(
                0, end_slope(slopes[0], final_slopes[0],
                             (values[1] - values[0]) / (grid[1] - grid[0])))
            final_slopes.append(
                end_slope(slopes[-1], final_slopes[-1],
                          (values[-1] - values[-2]) / (grid[-1] - grid[-2])))

        CubicHermite.__init__(self, values, final_slopes, grid)


class MonotoneCubic(PiecewiseMonotoneCubic):
    """Monotone cubic curve, see __init__()."""

    def __init__(self, values, *args, **kwargs):
        """Monotone cubic curve.

        This takes the same arguments as `PiecewiseMonotoneCubic`
        (except ``closed``), but it raises an error if the given values
        are not montone.

        See :ref:`/euclidean/piecewise-monotone.ipynb#Monotone-Interpolation`.

        """
        if 'closed' in kwargs:
            raise TypeError('The "closed" argument is not allowed')
        PiecewiseMonotoneCubic.__init__(self, values, *args, **kwargs)
        diffs = _np.diff(values)
        if not (all(diffs >= 0) or all(diffs <= 0)):
            raise ValueError('Only monotone values are allowed')

    # TODO: rename to something with "solve"?
    def get_time(self, value):
        """Get the time instance for the given value.

        If the solution is not unique (i.e. there is a plateau),
        ``None`` is returned.

        """
        if not _np.isscalar(value):
            return _np.array([self.get_time(v) for v in value])

        values = self.evaluate(self.grid)
        if values[0] <= value <= values[-1]:
            # Increasing values

            def get_index(values, value):
                return _bisect_right(values, value) - 1

        elif values[-1] <= value <= values[0]:
            # Decreasing values

            def get_index(values, value):
                return len(values) - _bisect_left(values[::-1], value) - 1

        else:
            raise ValueError(f'value outside allowed range: {value}')
        # First, check for exact matches to find plateaus
        matches, = _np.nonzero(values == value)
        if len(matches) > 1:
            return None
        if len(matches) == 1:
            return self.grid[matches[0]]

        idx = get_index(values, value)
        coeffs = self.segments[idx]
        # Solve for p - value = 0
        roots = (_np.poly1d(coeffs) - value).roots
        # Segment is only defined for t in [0, 1]
        roots = roots[_np.isreal(roots) & (roots >= 0) & (roots <= 1)]
        assert len(roots) == 1 and _np.isreal(roots)
        time, = roots.real
        t0, t1 = self.grid[idx:idx + 2]
        return time * (t1 - t0) + t0


class ConstantSpeedAdapter:
    """Re-parameterize a spline to have constant speed, see __init__()."""

    def __init__(self, curve):
        """Re-parameterize a spline to have constant speed.

        For splines in Euclidean space this amounts to arc-length
        parameterization.

        However, this class is implemented in a way that also allows using
        rotation splines which will be re-parameterized to have constant
        angular speed.

        The parameter *s* represents the cumulative arc-length or the
        cumulative rotation angle, respectively.

        """
        self.curve = curve
        lengths = (
            self._integrated_speed(i, t0, t1)
            for i, (t0, t1) in enumerate(zip(curve.grid, curve.grid[1:])))

        # NB: "initial" argument to itertools.accumulate since Python 3.8
        #self.grid = list(_accumulate(lengths, initial=0))
        self.grid = [0] + list(_accumulate(lengths))

    def _integrated_speed(self, idx, t0, t1):
        """Integral over the speed on a curve segment.

        *t0* and *t1* must be within the given segment.

        """
        if not 0 <= idx < len(self.curve.grid) - 1:
            raise ValueError(f'invalid idx: {idx}')
        if not self.curve.grid[idx] <= t0 <= t1 <= self.curve.grid[idx + 1]:
            raise ValueError('Invalid t0 or t1')

        def speed(t):
            return _np.linalg.norm(self.curve.evaluate(t, 1), axis=-1)

        from scipy import integrate
        value, abserr = integrate.quad(speed, t0, t1)
        return value

    def _s2t(self, s):
        """Convert integrated speed to time value."""
        idx = _check_param('s', s, self.grid)
        s -= self.grid[idx]
        t0 = self.curve.grid[idx]
        t1 = self.curve.grid[idx + 1]

        def length(t):
            return self._integrated_speed(idx, t0, t) - s

        from scipy.optimize import bisect
        return bisect(length, t0, t1)

    def evaluate(self, s):
        if not _np.isscalar(s):
            return _np.array([self.evaluate(s) for s in s])
        return self.curve.evaluate(self._s2t(s))


class NewGridAdapter:
    """Re-parameterize a spline with new grid values, see __init__()."""

    def __init__(self, curve, new_grid=1):
        """Re-parameterize a spline with new grid values.

        :param curve: A spline.
        :param new_grid: If a single number is given, the new parameter
            will range from 0 to that number.  Otherwise, a sequence
            of numbers has to be given, one for each grid value.
            Instead of a value, ``None`` can be specified to choose a
            value automatically.
            The first and last value cannot be ``None``.
        :type new_grid: optional

        """
        if _np.isscalar(new_grid):
            new_grid = [0] + [None] * (len(curve.grid) - 2) + [new_grid]
        if len(new_grid) != len(curve.grid):
            raise ValueError('new_grid must have same length as curve.grid')
        if new_grid[0] is None or new_grid[-1] is None:
            raise TypeError('first/last element of new_grid cannot be None')
        old_values, new_values = [], []
        for old, new in zip(curve.grid, new_grid):
            # TODO: allow NaN?
            if new is None:
                continue
            new_values.append(new)
            old_values.append(old)
        self._new2old = MonotoneCubic1D(old_values, grid=new_values)
        self.grid = []
        for old, new in zip(curve.grid, new_grid):
            if new is None:
                new = self._new2old.get_time(old)
            self.grid.append(new)
        self.curve = curve

    def evaluate(self, u):
        if not _np.isscalar(u):
            return _np.array([self.evaluate(u) for u in u])
        idx = _check_param('u', u, self.grid)
        return self.curve.evaluate(self._new2old.evaluate(u))
