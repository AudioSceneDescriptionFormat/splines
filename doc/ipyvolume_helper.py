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
