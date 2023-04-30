End Conditions
==============

Most spline types that are defined
by a sequence of control points to be interpolated
need some additional information
to be able to draw their segments at the beginning and at the end.
For example, cubic `Catmull--Rom splines`__
need four consecutive control points to define the segment between the middle two.
For the very first and last segment, the fourth control point is missing.
Another example are `natural splines`__,
which would require to solve an
underdetermined system of equations when only the control points are given.

__ catmull-rom.ipynb
__ natural.ipynb

There are many ways to provide this missing information,
here we will mention only a few of them.

clamped
    This means providing a fixed tangent (i.e. first derivative)
    at the beginning and end of a cubic spline.
    For higher degree splines, additional derivatives have to be be specified.

natural
    For a cubic spline, this means setting the second derivative
    at the beginning and end of the spline to zero
    and calculating the first derivative from that constraint,
    see :doc:`end-conditions-natural`.

closed
    This problem can also be solved by simply not having a begin and an end.
    When reaching the last control point,
    the spline can just continue at the first control point.
    For non-uniform splines an additional parameter interval has to be specified
    for the segment that's inserted between the end and the beginning.

For most splines in the `splines module`__,
*clamped*, *natural* and *closed* end conditions are available
via the ``endconditions`` argument.
Except for *closed*,
the end conditions can differ between the beginning and end of the spline.

__ ../python-module.rst

Additional information is available for
`end conditions of natural splines`__ and
`monotone end conditions`__.

__ natural-uniform.ipynb#End-Conditions
__ piecewise-monotone.ipynb#End-Conditions


.. toctree::

    end-conditions-natural
