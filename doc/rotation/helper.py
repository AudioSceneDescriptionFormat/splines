"""Helper functions for plotting rotations."""
from functools import partial
from math import radians, degrees, atan2, asin

import matplotlib
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LightSource
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np

from splines.quaternion import UnitQuaternion


shade_colors = partial(Axes3D._shade_colors, 'dummy')
generate_normals = partial(Axes3D._generate_normals, 'dummy')


def faces():
    """Quadratic faces for an F-shaped object."""
    top = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1], [1, -1, 1]])
    north = np.array([[1,  1, -1], [-1,  1, -1], [-1,  1,  1], [1,  1,  1]])
    east = np.array([[1,  1, -1], [1,  1,  1], [1, -1,  1], [1, -1, -1]])
    south = np.array([[1, -1,  1], [-1, -1,  1], [-1, -1, -1], [1, -1, -1]])
    west = np.array([[-1,  1,  1], [-1,  1, -1], [-1, -1, -1], [-1, -1,  1]])
    bottom = np.array([[-1,  1, -1], [1,  1, -1], [1, -1, -1], [-1, -1, -1]])
    origin = np.array([0, 0, 0])
    p = origin + [-2, 0, 0]
    yield p + top
    yield p + west
    yield p + bottom
    p += [2, 0, 0]
    yield p + top
    yield p + north
    yield p + bottom
    yield p + south
    p += [2, 0, 0]
    yield p + top
    yield p + north
    yield p + bottom
    yield p + south
    yield p + east
    p = origin + [-2, 2, 0]
    yield p + top
    yield p + west
    yield p + bottom
    yield p + east
    p += [0, 2, 0]
    yield p + top
    yield p + west
    yield p + bottom
    yield p + north
    p += [2, 0, 0]
    yield p + top
    yield p + north
    yield p + bottom
    yield p + south
    p += [2, 0, 0]
    yield p + top
    yield p + north
    yield p + bottom
    yield p + south
    yield p + east
    p = origin + [-2, -2, 0]
    yield p + top
    yield p + west
    yield p + bottom
    yield p + east
    p += [0, -2, 0]
    yield p + top
    yield p + west
    yield p + bottom
    yield p + east
    yield p + south


def create_polys(rot, *, ls=None):
    if not isinstance(rot, UnitQuaternion):
        rot = UnitQuaternion.from_unit_xyzw(rot)
    polys = np.array([list(map(rot.rotate_vector, face)) for face in faces()])
    if ls is None:
        ls = LightSource()
    color = 'white'
    facecolors = shade_colors(color, generate_normals(polys), ls)
    return polys, facecolors


def create_empty_collection(ax):
    alpha = 1
    linewidth = 0.5
    edgecolor = 'black'
    coll = Poly3DCollection(
        [],
        closed=True,
        alpha=alpha,
        linewidth=linewidth,
        edgecolor=edgecolor,
    )
    ax.add_collection3d(coll)
    return coll


def prepare_axes(ax):
    size = 12
    _, _, x1, y1 = ax.bbox.bounds
    aspect = x1 / y1
    if x1 > y1:
        # landscape
        height = size
        width = height * aspect
    else:
        width = size
        height = width / aspect
    ax.set_xlim(-width / 2, width / 2)
    ax.set_ylim(-height / 2, height / 2)
    return create_empty_collection(ax)


def plot_rotation(rot, figsize=None, ls=None):
    if not isinstance(rot, dict):
        rot = {'': rot}
    collections = prepare_figure(rot.keys(), figsize=figsize)
    collections = update_collections(collections, rot.values(), ls=ls)
    return collections


def prepare_figure(titles='', **kwargs):
    if isinstance(titles, str):
        titles = [titles]
    fig, (axs,) = plt.subplots(
        ncols=len(titles),
        squeeze=False,
        subplot_kw=dict(projection='dumb3d'),
        **kwargs)
    collections = []
    for ax, title in zip(axs, titles):
        collections.append(prepare_axes(ax))
        ax.set_title(title)
    return collections


def update_collections(collections, rotations, *, ls=None):
    if ls is None:
        ls = LightSource()
    for coll, rot in zip(collections, rotations):
        polys, facecolors = create_polys(rot, ls=ls)
        coll.set_verts(polys)
        coll.set_facecolors(facecolors)
    return collections


def plot_rotations(rotations, *, figsize=None, ls=None):
    if not isinstance(rotations, dict):
        rotations = {'': rotations}
    fig, axs = plt.subplots(
        nrows=len(rotations),
        squeeze=False,
        figsize=figsize,
        subplot_kw=dict(projection='dumb3d'),
    )
    if ls is None:
        ls = LightSource()
    for [ax], [title, rotations] in zip(axs, rotations.items()):
        ax.set_title(title)
        object_width = 12
        shift_x = object_width
        _, _, x1, y1 = ax.bbox.bounds
        aspect = x1 / y1
        total_width = object_width * len(rotations)
        total_height = total_width / aspect
        ax.set_xlim(0, total_width)
        ax.set_ylim(0, total_height)
        x = object_width / 2
        y = total_height / 2
        z = 0
        for i, rot in enumerate(rotations):
            coll = create_empty_collection(ax)
            offset = [x + i * shift_x, y, z]
            polys, facecolors = create_polys(rot, ls=ls)
            polys += offset
            coll.set_verts(polys)
            coll.set_facecolors(facecolors)
    plt.show(fig)


def animate_rotations(rotations, *, figsize=None, ls=None, interval=40,
                      **kwargs):
    if not isinstance(rotations, dict):
        rotations = {'': rotations}
    if ls is None:
        ls = LightSource()
    collections = prepare_figure(rotations.keys(), figsize=figsize)
    plt.close(collections[0].axes.figure)

    def ani_func(rot):
        return update_collections(collections, rot, ls=ls)

    return FuncAnimation(
        collections[0].axes.figure,
        ani_func,
        frames=list(zip(*rotations.values())),
        interval=interval,
        **kwargs)


def display_animation(ani, **kwargs):
    from IPython.display import display
    display({
        'text/html': ani.to_jshtml(**kwargs),
        'text/plain': 'Animations can only be shown in HTML output, sorry!',
    }, raw=True)


class DumbAxes3D(Axes3D):

    name = 'dumb3d'

    def __init__(self, figure, rect=None, sharex=None, sharey=None):
        if sharex is not None:
            raise TypeError('sharex not supported')
        if sharey is not None:
            raise TypeError('sharey not supported')
        try:
            # "auto_add_to_figure" was added in Matplotlib 3.4.0,
            # its default will change to False in Matplotlib 3.5.
            # Once Matplotlib 3.4.x is no longer supported,
            # this work-around can be removed.
            super().__init__(figure, auto_add_to_figure=False, rect=rect)
        except AttributeError:
            super().__init__(figure, rect=rect)
        self.set_axis_off()
        self.set_figure(figure)
        self.disable_mouse_rotation()

    @matplotlib.artist.allow_rasterization
    def draw(self, renderer):
        self.patch.draw(renderer)  # background (axes.facecolor)
        xmin, xmax = self.get_xlim3d()
        ymin, ymax = self.get_ylim3d()
        zmin, zmax = self.get_zlim3d()
        # NB: z limits are deliberately swapped to switch to right-handed
        #     coordinates, z pointing out of the figure (towards the viewer)
        self.M = proj3d.world_transformation(xmin, xmax,
                                             ymin, ymax,
                                             zmax, zmin)
        self.vvec = NotImplemented
        self.eye = NotImplemented

        renderer.M = self.M
        renderer.vvec = self.vvec
        renderer.eye = self.eye
        renderer.get_axis_position = NotImplemented

        for coll in self.collections:
            try:
                # The "renderer" argument has been deprecated
                coll.do_3d_projection()
            except TypeError:
                coll.do_3d_projection(renderer)
        for patch in self.patches:
            try:
                # The "renderer" argument has been deprecated
                patch.do_3d_projection()
            except TypeError:
                patch.do_3d_projection(renderer)

        super(Axes3D, self).draw(renderer)

    def apply_aspect(self, position=None):
        pass


matplotlib.projections.register_projection(DumbAxes3D)


def angles2quat(azimuth, elevation, roll):
    return (
        UnitQuaternion.from_axis_angle((0, 0, 1), radians(azimuth)) *
        UnitQuaternion.from_axis_angle((1, 0, 0), radians(elevation)) *
        UnitQuaternion.from_axis_angle((0, 1, 0), radians(roll))
    )


# See https://AudioSceneDescriptionFormat.readthedocs.io/quaternions.html
def quat2angles(quat):
    a, b, c, d = quat.wxyz
    sin_elevation = 2 * (a * b + c * d)
    if 0.999999 < sin_elevation:
        # elevation ~= 90°
        return (
            degrees(atan2(2 * (a * c + b * d), 2 * (a * b - c * d))),
            90,
            0)
    elif sin_elevation < -0.999999:
        # elevation ~= -90°
        return (
            degrees(atan2(-2 * (a * c + b * d), 2 * (c * d - a * b))),
            -90,
            0)
    return (
        degrees(atan2(2 * (a * d - b * c), 1 - 2 * (b**2 + d**2))),
        degrees(asin(sin_elevation)),
        degrees(atan2(2 * (a * c - b * d), 1 - 2 * (b**2 + c**2))))
