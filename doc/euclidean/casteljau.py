import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np


def affine_combination(pair, t):
    """Return intermediate point between two points."""
    x0, x1 = np.asarray(pair)
    return (1 - t) * x0 + t * x1


def casteljau_step(points, t):
    """One step of De Casteljau's algorithm.

    Takes a list of points.
    Returns a list of intermediate points.

    """
    return [affine_combination(pair, t)
            for pair in zip(points[:-1], points[1:])]


def plot_casteljau(points, t, ax=None):
    """Plot steps of De Casteljau's algorithm.

    t is supposed to be between 0 and 1.

    """
    if len(points) < 2:
        raise ValueError('We need at least two points')
    if ax is None:
        ax = plt.gca()
    ax.set_title(f'Bézier curve of degree {len(points) - 1}; t = {t:1.2f}')
    while len(points) >= 2:
        ax.plot(*np.asarray(points).T)
        points = casteljau_step(points, t)
        ax.scatter(*np.asarray(points).T)
    ax.axis('equal')
    result, = points
    return result


def create_animation(points, frames, ax=None, **kwargs):
    """Create matplotlib animation for De Casteljau's algorithm.

    ``**kwargs`` are passed to ``FuncAnimation()``.

    """
    if ax is None:
        ax = plt.gca()
    partial_curve = []

    def animation_func(t):
        ax.clear()
        if partial_curve:
            ax.plot(*np.asarray(partial_curve).T, '.', c='lightgrey')
        point = plot_casteljau(points, t, ax=ax)
        ax.scatter(*np.asarray(points).T, marker='x', c='black')
        partial_curve.append(point)

    times = np.linspace(0, 1, frames)
    return FuncAnimation(ax.figure, animation_func, frames=times, **kwargs)
