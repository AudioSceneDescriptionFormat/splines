"""Helper functions for plotting."""
import matplotlib.pyplot as _plt
from matplotlib.collections import LineCollection
import numpy as _np
import sympy as _sp


def plot_slopes_1d(slopes, values, grid, scale=1, ax=None, **kwargs):
    """Plot incoming and outging slopes for 1D spline."""
    if ax is None:
        ax = _plt.gca()
    slopes = _np.asarray(slopes)
    for slopes, values, grid, pivot in (
            [slopes[::2], values[:-1], grid[:-1], 'tail'],
            [slopes[1::2], values[1:], grid[1:], 'tip']):
        lengths = _np.sqrt(1 + slopes**2)
        ax.quiver(
            grid, values, 1 / lengths, slopes / lengths,
            scale=scale, scale_units='x', angles='xy', color='lightgrey',
            pivot=pivot, **kwargs)


def plot_spline_1d(spline, ax=None, samples=100, **kwargs):
    """Plot a one-dimensional spline."""
    if ax is None:
        ax = _plt.gca()
    times = _np.linspace(spline.grid[0], spline.grid[-1], samples)
    ax.plot(times, spline.evaluate(times), **kwargs)
    ax.scatter(spline.grid, spline.evaluate(spline.grid))


def plot_tangent_2d(tangent, vertex, color='lightgrey', outgoing=True,
                    scale=1, ax=None, **kwargs):
    """Plot outgoing or incoming 2D tangent."""
    if ax is None:
        ax = _plt.gca()
    ax.quiver(
        *vertex, *tangent,
        scale=scale, scale_units='xy', angles='xy', color=color,
        pivot='tail' if outgoing else 'tip', **kwargs)
    if outgoing:
        endpoint = _np.add(vertex, tangent)
    else:
        endpoint = _np.subtract(vertex, tangent)
    # Plot an invisible point at the end of the tangent vector
    # to make sure the vector is visible when the plot is autoscaled.
    # NB: Selecting a (unused) color to not disturb the color cycle.
    ax.scatter(*endpoint, marker='', color='green')


def plot_tangents_2d(tangents, vertices, color='lightgrey',
                     scale=1, ax=None, **kwargs):
    """Plot outgoing and incoming tangents for 2D spline."""
    if ax is None:
        ax = _plt.gca()
    tangents = _np.asarray(tangents)
    vertices = _np.asarray(vertices)
    for i in range(len(vertices) - 1):
        plot_tangent_2d(tangents[2 * i], vertices[i], color=color, **kwargs)
        plot_tangent_2d(
            tangents[2 * i + 1], vertices[i + 1], color=color, outgoing=False,
            **kwargs)


def plot_spline_2d(spline, dots_per_second=15, marker='.', linestyle='',
                   chords=True, ax=None, **kwargs):
    """Plot a two-dimensional spline."""
    if ax is None:
        ax = _plt.gca()
    total_duration = spline.grid[-1] - spline.grid[0]
    dots = int(total_duration * dots_per_second) + 1
    times = spline.grid[0] + _np.arange(dots) / dots_per_second
    ax.plot(
        *spline.evaluate(spline.grid).T,
        color='lightgrey',
        linestyle=(':' if chords else ''),
        marker='x',
        markeredgecolor='black',
    )
    ax.plot(
        *spline.evaluate(times).T,
        marker=marker,
        linestyle=linestyle,
        **kwargs)
    ax.axis('equal')


def grid_lines(x=None, y=None, ax=None):
    if ax is None:
        ax = _plt.gca()
    if x is not None:
        ax.set_xticks(x)
        ax.xaxis.grid(True)
    if y is not None:
        ax.set_yticks(y)
        ax.yaxis.grid(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')


def latexify(expr):
    """Convert SymPy expression to LaTeX string.

    Strings are passed through unchanged.

    """
    if isinstance(expr, str):
        return expr
    # \boldsymbol is not available, see:
    # https://github.com/matplotlib/matplotlib/issues/1366
    return _sp.latex(expr, mode='inline').replace(r'\boldsymbol', r'\mathbf')


def plot_sympy(*args, ax=None, **kwargs):
    """Plot a SymPy expression into a Matplotlib plot."""
    if ax is None:
        ax = _plt.gca()
    for line in _sp.plot(*args, show=False, **kwargs):
        # NB: line.get_points() works only for smooth plots
        segments = line.get_segments()
        # Dummy plot to use default color cycle
        dummy, = ax.plot([])
        ax.add_collection(LineCollection(segments, color=dummy.get_color()))
        dummy.remove()
    ax.autoscale()


def plot_basis(*args, ax=None, parameter=_sp.Symbol('t'), labels=None):
    """Plot a polynomial basis (given as SymPy expressions)."""
    if ax is None:
        ax = _plt.gca()
    # Alternatively, plot_sympy() could be used, but using LineCollection
    # would inhibit automatic placement of the legend
    for line in _sp.plot(*args, (parameter, 0, 1), show=False):
        x, y = line.get_points()
        # if the function is constant, SymPy only emits one value:
        x, y = _np.broadcast_arrays(x, y)
        ax.plot(x, y)
    ax.autoscale()
    grid_lines([0, 1], [0, 1], ax=ax)
    if labels:
        ax.legend([latexify(l) for l in labels])
    ax.set_xlabel(latexify(parameter))
    ax.set_ylabel('weight')
