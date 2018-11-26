import ipyvolume as ipv
import numpy as np


def mesh_coordinates(l1, l2):
    """Coordinates of an L-shaped mesh.

    Length of shorter part: *l1*, length of longer part: *l2*.
    Thickness: 1 unit.

    """
    return np.array([
        # L:
        [l1 - 0.5, -0.5,     -0.5],
        [l1 - 0.5, -0.5,      0.5],
        [    -0.5, -0.5,     -0.5],
        [     0.5, -0.5,      0.5],
        [    -0.5, -0.5, l2 - 0.5],
        [     0.5, -0.5, l2 - 0.5],
        # top:
        [    -0.5, -0.5, l2 - 0.5],
        [    -0.5,  0.5, l2 - 0.5],
        # down and along bottom:
        [    -0.5, -0.5,     -0.5],
        [    -0.5,  0.5,     -0.5],
        [l1 - 0.5, -0.5,     -0.5],
        [l1 - 0.5,  0.5,     -0.5],
        # other end:
        [l1 - 0.5,  0.5,     -0.5],
        [l1 - 0.5,  0.5,      0.5],
        # other side of L:
        [    -0.5,  0.5,     -0.5],
        [     0.5,  0.5,      0.5],
        [    -0.5,  0.5, l2 - 0.5],
        [     0.5,  0.5, l2 - 0.5],
        # top again:
        [     0.5, -0.5, l2 - 0.5],
        [     0.5,  0.5, l2 - 0.5],
        # down the inner edge:
        [     0.5, -0.5,      0.5],
        [     0.5,  0.5,      0.5],
        [l1 - 0.5, -0.5,      0.5],
        [l1 - 0.5,  0.5,      0.5],
        # final triangle:
        [l1 - 0.5, -0.5,      0.5],
        [l1 - 0.5, -0.5,     -0.5],
    ])


def transform_with_frenet_frames(coords, spline):
    xs = []
    ys = []
    zs = []
    # TODO: argument for time step
    times = np.linspace(spline.grid[0], spline.grid[-1], 50)
    for t in times:
        transformed = np.array(coords, copy=True)
        T, N, B = spline.tangent(t), spline.normal(t), spline.binormal(t)
        if not np.isfinite(B).all():
            # Frenet frame is undefined, we reduce the object to a single point
            transformed *= 0
            B = np.zeros_like(B)
        transformed = (
            spline.evaluate(t)
            + transformed[:, 0:1] * T
            + transformed[:, 1:2] * N
            + transformed[:, 2:3] * B
        )
        x, y, z = [coord.reshape(-1, 2) for coord in transformed.T]
        xs.append(x)
        ys.append(y)
        zs.append(z)
    # TODO: more elegant way to stack arrays?
    x = np.array(xs)
    y = np.array(ys)
    z = np.array(zs)
    return x, y, z


def frenet_frame_animation(mesh, spline):
    # TODO: args for width and height
    fig = ipv.figure(width=640, height=480)
    surface = ipv.plot_surface(*transform_with_frenet_frames(mesh, spline))
    # TODO: argument for grid spacing?
    grid = np.linspace(spline.grid[0], spline.grid[-1], 50, endpoint=True)
    # TODO: dots instead of line?
    ipv.plot(*spline.evaluate(grid).T)
    ipv.scatter(*spline.evaluate(spline.grid).T, color='black', marker='sphere')
    ipv.style.box_off()
    ipv.squarelim()
    ipv.animation_control(surface)
    ipv.show()
