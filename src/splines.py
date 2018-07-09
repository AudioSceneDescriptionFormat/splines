"""Piecewise Polynomial Curves.

"""
from bisect import bisect_right as _bisect_right
import numpy as _np

__version__ = '0.0'


class PiecewiseCurve:
    r"""

    Arbitrary degree, arbitrary dimension.

    Uses monomial basis.

    .. math::

        \boldsymbol{p}_i(t) = \sum_{k=0}^n
            \boldsymbol{a}_k \left(\frac{t - t_i}{t_{i+1} - t_i}\right)^k
            \text{ for } t_i \leq t \lt t_{i+1}

    Similar to https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PPoly.html,
    which states:

        "High-order polynomials in the power basis can be numerically
        unstable.  Precision problems can start to appear for orders
        larger than 20-30."

    """

    def __init__(self, segments, grid):
        """

        *grid* must be strictly increasing.

        Theoretically, different segments could have different
        polynomial degree.

        """
        self.segments = [_np.array(coefficients, copy=True)
                         for coefficients in segments]
        self.grid = list(grid)

    def __call__(self, t, n=0):
        """Get value (or n-th derivative) at given time(s)."""
        if not _np.isscalar(t):
            return _np.array([self(time, n) for time in t])

        if t < self.grid[0]:
            raise ValueError(f't too small: {t}')
        elif t < self.grid[-1]:
            idx = _bisect_right(self.grid, t) - 1
        elif t == self.grid[-1]:
            idx = len(self.segments) - 1
        else:
            raise ValueError(f't too big: {t}')

        t0, t1 = self.grid[idx:idx + 2]
        t = (t - t0) / (t1 - t0)
        coefficients = self.segments[idx][:-n or None]
        powers = _np.arange(len(coefficients))[::-1]
        product = _np.multiply.reduce
        weights = product([powers + 1 + i for i in range(n)]) / (t1 - t0)**n
        return t**powers * weights @ coefficients


def _check_vertices(vertices, *, closed):
    """For closed curves, append first vertex at the end."""
    if len(vertices) < 3:
        raise ValueError('At least three vertices are required')
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
            grid.append(grid[-1] + _np.linalg.norm(x1 - x0)**alpha)
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


def _end_tangent(condition, vertices, times, other_tangent):
    if condition == 'natural':
        tangent = _natural_tangent(vertices, times, other_tangent)
    elif _np.shape(condition) == _np.shape(vertices[0]):
        tangent = condition
    else:
        raise ValueError(
            f'{condition!r} is not a valid start/end condition')
    return tangent


def _cubic_hermite_matrix(delta):
    return _np.array([[ 2, -2,      delta,  delta ],
                      [-3,  3, -2 * delta, -delta ],
                      [ 0,  0,      delta,      0 ],
                      [ 1,  0,          0,      0.]])


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


class CubicHermite(PiecewiseCurve):

    def __init__(self, vertices, tangents, grid=None):
        """

        list of vertices
        list of tangent vectors (two per segment, outgoing and incoming)
        list of times

        """
        if len(vertices) < 2:
            raise ValueError('At least 2 vertices are needed')
        if len(tangents) != 2 * (len(vertices) - 1):
            raise ValueError('Exactly 2 tangents per segment are needed')
        if grid is None:
            grid = range(len(vertices))
        if len(vertices) != len(grid):
            raise ValueError('As many grid times as vertices are needed')
        segments = [
            _cubic_hermite_matrix(t1 - t0) @ [x0, x1, v0, v1]
            for (x0, x1), (v0, v1), (t0, t1) in zip(
                zip(vertices, vertices[1:]),
                zip(tangents[::2], tangents[1::2]),
                zip(grid, grid[1:]))]
        PiecewiseCurve.__init__(self, segments, grid)


class CatmullRom(CubicHermite):

    # NB: Catmull-Rom could be implemented as special case of Kochanek-Bartels,
    #     but here we chose not to.

    # NB: We could use the characteristic matrix for Catmull-Rom splines, but
    #     this wouldn't work if only 3 vertices are given by the user.
    #     Since we have to handle this special case anyway, we use the same
    #     method for everything.  Apart from reducing the amount of code, this
    #     also allows us to define derived classes that overwrite
    #     _calculate_tangent(), e.g. FiniteDifference.

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
        """

        endconditions: 'closed', 'natural' or pair of tangent vectors
                        (a.k.a. "clamped")

        If 'closed', the first vertex is re-used as last vertex and an
        additional grid time has to be specified.

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
        if closed:
            # Move last (outgoing) tangent to the beginning:
            tangents = tangents[-1:] + tangents[:-1]
        else:
            tangents.insert(0, _end_tangent(
                start, vertices[:2], grid[:2], tangents[0]))
            tangents.append(_end_tangent(
                end, vertices[-2:], grid[-2:], tangents[-1]))
        CubicHermite.__init__(self, vertices, tangents, grid)


class FiniteDifference(CatmullRom):

    @staticmethod
    def _calculate_tangent(points, times):
        x_1, x0, x1 = _np.asarray(points)
        t_1, t0, t1 = times
        return ((x0 - x_1) / (t0 - t_1) + (x1 - x0) / (t1 - t0)) / 2


class KochanekBartels(CubicHermite):

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
        """

        TCB values can only be given for the interior vertices.

        """
        closed = endconditions == 'closed'
        if closed:
            interior = len(vertices)
        else:
            interior = len(vertices) - 2
        vertices = _check_vertices(vertices, closed=closed)
        grid = _check_grid(grid, alpha, vertices)
        tcb = _np.asarray(tcb)
        if tcb.ndim == 1 and len(tcb) == 3:
            tcb = _np.tile(tcb, (interior, 1))
        if len(tcb) != interior:
            raise ValueError('There must be two more vertices than TCB values '
                             '(except for closed curves)')
        start, end, zip_vertices, zip_grid = _check_endconditions(
            endconditions, vertices, grid)
        if closed:
            # Move first TCB value to the end:
            tcb = _np.roll(tcb, -1, axis=0)
        tangents = [
            tangent
            for points, times, tcb in zip(zip_vertices, zip_grid, tcb)
            for tangent in self._calculate_tangents(points, times, tcb)]
        if closed:
            # Move last (outgoing) tangent to the beginning:
            tangents = tangents[-1:] + tangents[:-1]
        else:
            tangents.insert(0, _end_tangent(
                start, vertices[:2], grid[:2], tangents[0]))
            tangents.append(_end_tangent(
                end, vertices[-2:], grid[-2:], tangents[-1]))
        CubicHermite.__init__(self, vertices, tangents, grid)


def _monotone_end_condition(inner_slope, chord_slope):
    """

    Return the "outer" (i.e. first or last) slope given the "inner"
    (i.e. second or penultimate) slope and the slope of the
    corresponding chord.

    """
    # NB: This is a very ad-hoc algorithm meant to minimize the change in slope
    # within the first/last curve segment.  Especially, this should avoid a
    # change from negative to positive acceleration (and vice versa).
    # There might be a better method available!?!
    if chord_slope >= 0:
        assert 0 <= inner_slope <= 3 * chord_slope
        if inner_slope <= chord_slope:
            return 3 * chord_slope - 2 * inner_slope
        else:
            return (3 * chord_slope - inner_slope) / 2
    else:
        assert 3 * chord_slope <= inner_slope <= 0
        if chord_slope <= inner_slope:
            return 3 * chord_slope - 2 * inner_slope
        else:
            return (3 * chord_slope - inner_slope) / 2


class ShapePreservingCubic1D(FiniteDifference):
    """

    similar to https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.PchipInterpolator.html

    """

    def __init__(self, values, grid=None, slopes=None, *, alpha=None,
                 closed=False):
        """

        *slopes* is a list of slopes or None if slope should be computed
        from neighboring values.

        For undefined slopes, _calculate_tangent() is called on the base
        class.

        If no slopes are given, the curve also preserves
        concavity/convexity, otherwise it only preserves monotonicity
        and local extrema.

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
            """Manipulate the slope to preserve shape.

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
        elif not final_slopes:
            chord = (values[1] - values[0]) / (grid[1] - grid[0])
            one, two = slopes

            def check_slope(slope):
                if slope != fix_slope(slope, chord, chord):
                    raise ValueError(f'Slope too steep or wrong sign: {slope}')

            if one is None:
                if two is None:
                    final_slopes = [chord] * 2
                else:
                    check_slope(two)
                    final_slopes = [_monotone_end_condition(two, chord), two]
            else:
                if two is None:
                    check_slope(one)
                    final_slopes = [one, _monotone_end_condition(one, chord)]
                else:
                    check_slope(one)
                    check_slope(two)
                    final_slopes = [one, two]
        else:

            def end_slope(outer, inner, chord):
                if outer is None:
                    outer = _monotone_end_condition(inner, chord)
                else:
                    if outer != fix_slope(outer, chord, chord):
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

    # TODO: rename to something with "solve"?
    def get_time(self, value):
        """Get the time instance for the given value.

        If the solution is not unique (i.e. non-monotonic or repeated
        values), None is returned.

        """
        values = self(self.grid)
        if not values[0] <= value <= values[-1]:
            raise ValueError(f'value outside allowed range: {value}')
        # First, check for exact matches to avoid accuracy issues
        matches, = _np.nonzero(values == value)
        if len(matches) > 1:
            return None  # non-monotonic or plateau
        if len(matches) == 1:
            return self.grid[matches[0]]

        idx = _bisect_right(values, value) - 1
        coeffs = self.segments[idx]
        # Solve for p - value = 0
        roots = (_np.poly1d(coeffs) - value).roots
        # Segment is only defined for t in [0, 1]
        roots = roots[_np.isreal(roots) & (roots >= 0) & (roots <= 1)]
        assert len(roots) == 1 and _np.isreal(roots)
        time, = roots.real
        t0, t1 = self.grid[idx:idx + 2]
        return time * (t1 - t0) + t0
