"""Helper functions for plotting."""
import matplotlib.pyplot as _plt
import numpy as _np


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
