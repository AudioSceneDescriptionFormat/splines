"""Visualization of the Barry/Goldman (1988) algorithm."""
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
import numpy as np

from helper import plot_x_3_to_6


def lerp(xs, ts, t):
    """Linear interpolation."""
    x_begin, x_end = map(np.asarray, xs)
    t_begin, t_end = ts
    if not np.isscalar(t):
        t = np.expand_dims(t, axis=-1)
    return (x_begin * (t_end - t) + x_end * (t - t_begin)) / (t_end - t_begin)


def animation(points, times, frames=30, interval=200):
    points = np.asarray(points)
    x_1, x0, x1, x2 = points
    t_1, t0, t1, t2 = times

    def p_10(t):
        return lerp((x_1, x0), (t_1, t0), t)

    def p01(t):
        return lerp((x0, x1), (t0, t1), t)

    def p12(t):
        return lerp((x1, x2), (t1, t2), t)

    def p_101(p_10, p01, t):
        return lerp((p_10, p01), (t_1, t1), t)

    def p012(p01, p12, t):
        return lerp((p01, p12), (t0, t2), t)

    def x01(p_101, p012, t):
        return lerp((p_101, p012), (t0, t1), t)

    c0, c1, c2 = plt.rcParams['axes.prop_cycle'].by_key()['color'][:3]

    fig, ax = plt.subplots()
    ax.plot(*p_10([t0, t1]).T, color=c0)
    ax.plot(*p01([t0, t1]).T, color=c0)
    ax.plot(*p12([t0, t1]).T, color=c0)
    plot_x_3_to_6(points, ax)
    three_dots = ax.scatter([], [], color=c0)
    line_101, = ax.plot([], color=c1)
    line012, = ax.plot([], color=c1)
    two_dots = ax.scatter([], [], color=c1)
    line, = ax.plot([], color=c2)
    one_dot = ax.scatter([], [], color=c2)
    ax.axis('equal')
    plt.close(fig)

    def ani_func(t):
        p_10_v = p_10(t)
        p01_v = p01(t)
        p12_v = p12(t)
        three_dots.set_offsets([p_10_v, p01_v, p12_v])
        line_101.set_data(p_101(p_10_v, p01_v, [t0, t1]).T)
        line012.set_data(p012(p01_v, p12_v, [t0, t1]).T)
        p_101_v = p_101(p_10_v, p01_v, t)
        p012_v = p012(p01_v, p12_v, t)
        two_dots_data = np.column_stack([p_101_v, p012_v])
        two_dots.set_offsets(two_dots_data.T)
        line.set_data(two_dots_data)
        one_dot_data = x01(p_101_v, p012_v, t)
        one_dot.set_offsets(one_dot_data)
        ax.plot(*one_dot_data, '.', color='lightgrey')
        ax.set_title(f'Barry–Goldman algorithm; t = {t:1.2f}')

    frames = np.linspace(t0, t1, frames)

    return FuncAnimation(fig, ani_func, frames=frames, interval=interval)
