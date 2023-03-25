Splines
=======

The term *spline* for the mathematical description of a smooth piecewise curve
was introduced in :cite:`schoenberg1946smoothing`,
with reference to a drawing tool called spline_.

.. _spline: https://en.wiktionary.org/wiki/spline

..

    A spline is a simple mechanical device for drawing smooth curves.
    It is a slender flexible bar made of wood or some other elastic material.
    The spline is place[d] on the sheet of graph paper and held in place
    at various points by means of certain heavy objects
    (called "dogs" or "rats")
    such as to take the shape of the curve we wish to draw.

    -- :cite:`schoenberg1946smoothing`, page 67

The term is defined in the context of what is nowadays known as
`natural splines`_, especially cubic natural splines
(i.e. of degree 3; i.e. of order 4),
which have :math:`C^2` continuity.

..

    For :math:`k=4` they represent approximately the curves drawn by
    means of a spline and for this reason we propose to call them
    *spline curves of order* :math:`k`.

    -- :cite:`schoenberg1946smoothing`, page 48


Definition
----------

Different authors use different definitions for the term *spline*,
here is ours:
*splines* are composite parametric curves.
Splines are typically used for defining curves in
one-, two- or three-dimensional Euclidean space.
Such splines will be described in the following sections.
Later, we will also have a look at `rotation splines`__.

__ ../rotation/index.ipynb

Sometimes it is not obvious whether the term *spline*
refers to the composite curve or to one of its segments,
especially when talking about `Bézier splines`_.
In the rest of this text we are using the term *spline*
to refer to the entire composite curve.


Properties
----------

Different types of splines have different properties.
In the following, we list the most important properties,
focusing on the types of splines that will be described
in more detail in later sections.

piecewise
    Arguably the most important property of splines is that they are
    composed of somewhat independent pieces.
    This allows using simpler mathematical objects for the pieces,
    while still being able to construct a more or less arbitrarily
    complicated composite curve.
    For example,
    as shown in the previous section about `Lagrange interpolation`__,
    using a curve with a high polynomial degree can lead to
    unintended behavior like `Runge's phenomenon`__.
    This can be avoided by using multiple polynomial pieces of lower degrees.

    __ lagrange.ipynb
    __ lagrange.ipynb#Runge's-Phenomenon

parametric
    Here we are only talking about *univariate* curves,
    i.e. curves with one parameter,
    i.e. a single real number,
    but similar approaches can be used to describe
    surfaces with two parameters.
    We are normally using the symbol :math:`t` for the free parameter.
    This parameter can often intuitively be interpreted as *time*,
    but it doesn't have to.

    The individual segments of a spline are of course also parametric,
    and they may have their own parameter ranges
    (often, but not necessarily, the so-called *unit interval* from 0 to 1),
    which have to be calculated from the appropriate sub-range
    of the main spline parameter.

    A spline can also be re-parameterized, see :doc:`re-parameterization`.

    The sequence of parameter values at the start and end of segments
    is sometimes (e.g. in :cite:`gordon1974bspline`) called the *knot vector*.
    In the accompanying :doc:`/python-module`, however,
    it is called ``grid``.

non-uniform
    The parameter range of a spline can be uniquely separated
    into the parameter ranges of its segments.
    If those sub-ranges all have the same length,
    the spline is called *uniform*.

    When a uniform spline has curve segments of very different lengths,
    the speed along the curve
    (assuming that the parameter :math:`t` is interpreted as time)
    varies strongly.
    By using *non-uniform* parameter intervals, this can be avoided.

continuous
    Splines are not necessarily continuous.
    The segments of a spline might be defined by discontinuous functions,
    but for most practical applications
    it is more common to use continuous functions.
    Often, some derivatives of these functions are continuous as well.
    If the spline segments are polynomials, they are always continuous,
    and so are all their derivatives.

    However, even if its segments are continuous,
    that doesn't automatically mean that the whole spline is continuous.
    The transitions between segments can still be discontinuous.
    But again, in most practical applications the transitions are continuous.
    If that's the case, the spline is said to have
    :math:`C^0` *continuity*.
    The spline is called :math:`C^1` continuous
    if the first derivatives of the two neighboring segments
    at each transition are equal and
    :math:`C^2` continuous if also the second derivatives match.

control points
    Splines are fully defined by the mathematical functions of their segments
    and the corresponding parameter ranges.
    However, those functions (and their coefficients) have to be chosen somehow.
    And that's what differentiates different types of splines.

    For some applications it is desired to specify
    a sequence of *control points*
    (sometimes also called *vertex/vertices*)
    where the curve is supposed to pass through.
    Based on those points, the appropriate functions for the
    spline segments are constructed.
    The `Catmull--Rom splines`_ and `natural splines`_ are examples
    where segments are derived from such a sequence of control points.

    Some splines, most notably `Bézier splines`_,
    only pass through some of their control points and
    the remaining control points affect the shape of the curve
    between those points.

    The set of all control points, connected by straight lines,
    is sometimes called *control polygon*.

    Some splines have a set of control points where they pass through
    and additional values that are not points at all.
    We call them *control values*.
    For example, `Hermite splines`_ pass through a set of control points,
    but they need additional information about the tangent vectors
    (i.e. the first derivatives) at the transitions between segments.
    For higher-order splines they also need the second and higher derivatives.

interpolating
    Splines are called *interpolating* if they
    pass through all of their aforementioned control points.
    If a spline is not interpolating,
    it is called *approximating*.

    Here we are almost exclusively talking about interpolating splines.
    A notable special case are `Bézier splines`_,
    which pass through a sequence of control points,
    but between each pair of those interpolated control points
    there are :math:`d - 1` (where :math:`d` is the degree)
    additional control points that are only approximated by the curve
    (and they can be used to control the shape of the curve).

local control
    For some types of splines,
    when changing a single control value, the shape of the whole curve changes.
    These splines are said to have *global control*.
    For many applications, however,
    it is preferable, when a control value is changed,
    that the shape of the curve only changes in the immediate vicinity of that
    control value.
    This is called *local control*.

additional parameters
    Some types of splines have additional parameters,
    either separately for each vertex, or the same one(s) for all vertices.
    An example are `Kochanek--Bartels splines`_ with their
    *tension*, *continuity* and *bias* parameters.

polynomial
    The curve segments that make up a spline can have an arbitrary mathematical
    description.
    Very often, polynomial curve segments are used,
    and that's also what we will be mostly using here.
    The polynomials will be defined by their basis functions and
    corresponding coefficients, as described in
    the `notebook about polynomial parametric curves`__.

    __ polynomials.ipynb

    The following properties are only relevant for polynomial splines.

degree
    The degree of a polynomial spline is the highest degree among its segments.
    Splines of degree 3, a.k.a *cubic* splines,
    are very common for drawing smooth curves.
    Old-school references like :cite:`de_boor1978splines`
    might use the term *order*, which is one more than the degree,
    which means that cubic splines are of order 4.
    We will mostly consider cubic splines,
    but some of the presented algorithms allow arbitrary degree,
    for example `De Casteljau's algorithm`__.

    __ bezier-de-casteljau.ipynb

non-rational
    The splines discussed here are defined by one polynomial per segment.
    However, there are also splines whose segments are defined by
    ratios of polynomials instead.
    Those are called *rational* splines.
    Rational splines are invariant under perspective transformations
    (non-rational splines are only invariant under rotation/scale/translation),
    and they can precisely define conic sections (e.g. circles).
    They are also the foundation for NURBS__.

    __ https://en.wikipedia.org/wiki/Non-uniform_rational_B-spline


Types
-----

There are an infinite number of types of splines,
only very few of which will be presented in the following sections.
Some of them can create the same curve from different control values,
like `Hermite splines`_ and `Bézier splines`_.
Some create different curves from the same control values,
like `Catmull--Rom splines`_ and `natural splines`_.
Some have additional parameters to control the shape of the curve,
like `Kochanek--Bartels splines`_ with their TCB values.

Some spline types have certain constraints
on the transitions between segments,
for example, natural splines require :math:`C^2` continuity.
Other splines have no such constraints,
like for example Hermite splines,
which allow specifying arbitrary derivatives at their segment transitions.

Cubic splines cannot be interpolating *and* have :math:`C^2` continuity *and*
local control at the same time.

======================= ============= =========== =============
type                    local control continuity  interpolating
======================= ============= =========== =============
`Catmull--Rom splines`_ yes           :math:`C^1` yes
`natural splines`_      no            :math:`C^2` yes
B-splines_              yes           :math:`C^2` no
======================= ============= =========== =============

.. _Hermite splines: hermite.ipynb
.. _Bézier splines: bezier.ipynb
.. _Catmull--Rom splines: catmull-rom.ipynb
.. _natural splines: natural.ipynb
.. _Kochanek--Bartels splines: kochanek-bartels.ipynb
.. _B-splines: https://en.wikipedia.org/wiki/B-spline

Kochanek--Bartels splines with :math:`C = 0`
are in the same category as Catmull--Rom splines
(which are a subset of former).

From any polynomial segment of a certain degree
the control values according to any polynomial spline type
(of that same degree)
can be computed and vice versa.
This means that different types of polynomial splines
can be unambiguously (if using the same parameter intervals)
converted between each other
as long as the target spline has the same or weaker constraints.
For example, any natural spline can be converted into
its corresponding Bézier spline.
The reverse is not true.
Catmull--Rom splines and natural splines
can generally not be converted between each other
because they have mutually incompatible constraints.
