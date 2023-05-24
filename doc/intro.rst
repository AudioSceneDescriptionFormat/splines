Introduction
============

This is the documentation for the
|splines|__ module for Python.
However, instead of a Python module with a bit of documentation,
this project is mostly documentation,
with a bit of Python module at the side.
The goal is not so much to provide a turn-key software for using splines,
but rather to provide the background and mathematical derivations
for fully understanding the presented types of splines
and their inter-relations.
The Python module serves mostly for experimenting further with the presented
ideas and methods.
Therefore, the implementation is not focused on efficiency.

.. |splines| replace:: ``splines``
__ https://pypi.org/project/splines/

The documentation consists of two main parts.
The :doc:`first part <euclidean/index>` investigates some *polynomial splines*
in their natural habitat, the Euclidean space.
In the unlikely case you are reading this
and don't know what "spline" means,
the first part also contains
:doc:`a definition of the term <euclidean/splines>`
and a description of some of the common properties of splines.
The :doc:`second part <rotation/index>` leaves the comfort zone of flat space
and tries to apply some of the approaches from the first part to
the curved space of rotations.
The Python module is similarly split into two parts
whose API documentation is available at
:mod:`splines` and :mod:`splines.quaternion`, respectively.

This project was originally inspired by :cite:t:`millington2009matrices`,
who concisely lists the *basis matrices*
(a.k.a. *characteristic matrices*)
of a few common types of splines
and also provides matrices that can be used to convert *control points*
between those different types.
However, the derivation of those matrices is not shown.
Furthermore, the paper only considers *uniform* curves,
where all parameter intervals have a size of 1.
One goal of this documentation is to show the derivation of
all equations and matrices.
The derivations often utilize SymPy_ to make them more reproducible
and to ease further experimentation.
A special focus is put on *non-uniform* splines,
which seem to have been neglected in some of the literature
and especially in some online resources.

.. _SymPy: https://www.sympy.org/

Another focus is the speed along curves.
In many applications only the shape
(a.k.a. the image_)
of a curve matters.
However, sometimes it is important how fast a point travels along a spline
when changing the parameter (which can be interpreted as *time*).
The "timing" of a spline is not visible in a static line plot
(as is often used by default).
That's why most of the plots in the following sections will instead use
dots at regular parameter intervals, for example 15 dots per second.
If a spline already has the desired image but the wrong timing,
this can be fixed by :doc:`euclidean/re-parameterization`.

.. _image: https://en.wikipedia.org/wiki/Image_(mathematics)

A non-goal of this Python module
and its documentation is
to cover all possible types of splines.
Maybe some additional types will be added in the future,
but the list will always stay incomplete.
One of the most glaring omissions for now are B-splines_,
which are mentioned a few times but not properly derived nor implemented.
Another family of splines that is missing are *rational splines*,
and therefore also their most popular member NURBS_.
Spline surfaces are not covered either.

.. _B-splines: https://en.wikipedia.org/wiki/B-spline
.. _NURBS: https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline
