"""Visualization of Catmull--Rom spline evaluation."""
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


def animation_2_1(points, times, frames=30, interval=200, blit=True):
    """Animation for quadratic Lagrange followed by linear blending."""
    points = np.asarray(points)
    x_1, x0, x1, x2 = points
    t_1, t0, t1, t2 = times

    def neville_the_first(t):
        return lerp(
            [lerp([x_1, x0], [t_1, t0], t), lerp([x0, x1], [t0, t1], t)],
            [t_1, t1], t)

    def neville_the_second(t):
        return lerp(
            [lerp([x0, x1], [t0, t1], t), lerp([x1, x2], [t1, t2], t)],
            [t0, t2], t)

    c0, c1 = plt.rcParams['axes.prop_cycle'].by_key()['color'][:2]

    fig, ax = plt.subplots()
    ax.plot(*neville_the_first(np.linspace(t_1, t1, 30)).T, color=c0)
    ax.plot(*neville_the_second(np.linspace(t0, t2, 30)).T, color=c0)
    plot_x_3_to_6(points, ax)
    two_dots = ax.scatter([], [], color=c0)
    line, = ax.plot([], color=c1)
    one_dot = ax.scatter([], [], color=c1)
    dots, = ax.plot([], '.', color='lightgrey')
    partial_curve = []
    ax.axis('equal')
    plt.close(fig)

    def ani_func(t):
        first_dot = neville_the_first(t)
        second_dot = neville_the_second(t)
        third_dot = lerp([first_dot, second_dot], [t0, t1], t)
        two_dots_data = np.column_stack([first_dot, second_dot])
        two_dots.set_offsets(two_dots_data.T)
        line.set_data(two_dots_data)
        one_dot.set_offsets(third_dot)
        partial_curve.append(third_dot)
        dots.set_data(np.array(partial_curve).T)
        return line, two_dots, one_dot, dots

    frames = np.linspace(t0, t1, frames)

    return FuncAnimation(
        fig, ani_func, frames=frames, interval=interval, blit=blit)


def animation_1_2(points, times, frames=30, interval=200, blit=True):
    """Animation for linear interpolations followed by quadratic B-spline."""
    points = np.asarray(points)
    x_1, x0, x1, x2 = points
    t_1, t0, t1, t2 = times

    def p_10(t):
        return lerp((x_1, x0), (t_1, t0), t)

    def p01(t):
        return lerp((x0, x1), (t0, t1), t)

    def p12(t):
        return lerp((x1, x2), (t1, t2), t)

    def de_boor(xs, t):
        return lerp(
            [lerp(xs[:2], [t_1, t1], t), lerp(xs[1:], [t0, t2], t)],
            [t0, t1], t)

    c0, c1 = plt.rcParams['axes.prop_cycle'].by_key()['color'][:2]

    fig, ax = plt.subplots()
    ax.plot(*p_10([t0, t1]).T, color=c0)
    ax.plot(*p01([t0, t1]).T, color=c0)
    ax.plot(*p12([t0, t1]).T, color=c0)
    plot_x_3_to_6(points, ax)
    three_dots = ax.scatter([], [], color=c0)
    chords, = ax.plot([], linestyle='dashed', linewidth=1, color='lightgrey')
    b_spline, = ax.plot([], color=c1)
    one_dot = ax.scatter([], [], color=c1)
    dots, = ax.plot([], '.', color='lightgrey')
    partial_curve = []
    ax.axis('equal')
    plt.close(fig)

    def ani_func(t):
        p_10_v = p_10(t)
        p01_v = p01(t)
        p12_v = p12(t)
        three_dots_data = np.array([p_10_v, p01_v, p12_v])
        chords.set_data(three_dots_data.T)
        b_spline.set_data(*de_boor(three_dots_data, np.linspace(t0, t1, 30)).T)
        three_dots.set_offsets(three_dots_data)
        one_dot_data = de_boor(three_dots_data, t)
        one_dot.set_offsets(one_dot_data)
        partial_curve.append(one_dot_data)
        dots.set_data(np.array(partial_curve).T)
        return three_dots, b_spline, one_dot, dots

    frames = np.linspace(t0, t1, frames)

    return FuncAnimation(
        fig, ani_func, frames=frames, interval=interval, blit=blit)
