import matplotlib.pyplot as plt
import numpy as np

import splines

from helper import plot_spline_2d


def plot_tcb(*tcb, ax=None):
    """Plot four TCB examples."""
    if ax is None:
        ax = plt.gca()
    vertices = [
        (-3.5, 0),
        (-1, 1.5),
        (0, 0.1),
        (1, 1.5),
        (3.5, 0),
        (1, -1.5),
        (0, -0.1),
        (-1, -1.5),
    ]
    for idx, tcb in zip([1, 7, 3, 5], tcb):
        all_tcb = np.zeros((len(vertices), 3))
        all_tcb[idx] = tcb
        s = splines.KochanekBartels(
            vertices, tcb=all_tcb, endconditions='closed')
        label = ', '.join(
            f'{name} = {value}'
            for name, value in zip('TCB', tcb)
            if value)
        plot_spline_2d(s, chords=False, label=label, ax=ax)
    plot_spline_2d(
        splines.KochanekBartels(vertices, endconditions='closed'),
        color='lightgrey', chords=False, ax=ax)
    lines = [l for l in ax.get_lines() if not l.get_label().startswith('_')]
    # https://matplotlib.org/tutorials/intermediate/legend_guide.html#multiple-legends-on-the-same-axes
    ax.add_artist(ax.legend(
        handles=lines[:2], bbox_to_anchor=(0, 0., 0.5, 1),
        loc='center', numpoints=3))
    ax.legend(
        handles=lines[2:], bbox_to_anchor=(0.5, 0., 0.5, 1),
        loc='center', numpoints=3)
