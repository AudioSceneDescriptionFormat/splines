from collections import Mapping
import sympy as sp
from sympy.core.containers import Dict


# https://github.com/sympy/sympy/issues/4986
# https://github.com/sympy/sympy/issues/5031


class NamedExpressionError(Exception):
    pass


class NamedExpression(sp.Equality):
    """Named expression.

    "name" is not very strictly a name ... can be a symbol or an
    arbitrary expression.

    """

    def __new__(cls, name, expr=None, **kwargs):
        self = object.__new__(cls)
        if isinstance(name, sp.Basic):
            assert not kwargs
            self._name = name
        else:
            self._name = sp.Symbol(name, **kwargs)
        self.expr = expr
        return self

    def __getattribute__(self, name):
        if name in ['__getitem__', '_mhash']:
            # We are hijacking these attribute accesses to spot usage errors
            raise NamedExpressionError('Use .name or .expr')
        return super().__getattribute__(name)

    @property
    def name(self):
        return self._name

    @property
    def expr(self):
        return self._expr

    @expr.setter
    def expr(self, value):
        self._expr = value

    lhs = name
    rhs = expr

    @property
    def _args(self):
        return self.name, self.expr

    def doit(self):
        return self.func(self.name.doit(), self.expr.doit())

    def subs(self, *args, **kwargs):
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, NamedExpression):
                args = [[(arg.name, arg.expr)]]
            else:
                new_arg = []

                if isinstance(arg, (Dict, Mapping)):
                    pass  # Nothing to do
                else:
                    try:
                        for s in arg:
                            if isinstance(s, NamedExpression):
                                new_arg.append((s.name, s.expr))
                            else:
                                new_arg.append(s)
                    except TypeError:
                        pass  # Nothing to do, SymPy will raise
                    else:
                        if isinstance(arg, set):
                            new_arg = set(new_arg)
                        args = [new_arg]
        elif len(args) == 2:
            pass  # Nothing to do
        return self.func(
            self.name.subs(*args, **kwargs),
            self.expr.subs(*args, **kwargs))

    def simplify(self, *args, **kwargs):
        return self.func(self.name, sp.simplify(self.expr, *args, **kwargs))

    def factor(self, *args, **kwargs):
        return self.func(self.name, sp.factor(self.expr, *args, **kwargs))

    def expand(self, *args, **kwargs):
        return self.func(self.name, sp.expand(self.expr, *args, **kwargs))

    def pull_out(self, expr):
        # NB: This ignores the subclass and just returns a NamedExpression:
        return NamedExpression(
            self.name,
            sp.Mul(expr, sp.simplify(self.expr / expr), evaluate=False))


class NamedMatrix(NamedExpression):

    def __new__(cls, *args, **kwargs):
        """

        NamedMatrix(name, n, m)
        NamedMatrix(name, n, m, expr)
        NamedMatrix(name)  # if name is already a matrix
        NamedMatrix(name, expr)  # if name or expr is already a matrix

        """
        self = object.__new__(cls)
        if not (1 <= len(args) <= 4):
            TypeError('1 to 4 positional arguments are required')
        args = list(args)
        name = args.pop(0)
        if isinstance(name, (sp.MatrixBase, sp.MatrixExpr)):
            if len(args) >= 2:
                n = args.pop(0)
                m = args.pop(0)
                if (n, m) != name.shape:
                    raise ValueError('Shape mismatch')
            if kwargs:
                raise TypeError('No kwargs allowed if name is already matrix')
            self._name = name
        else:
            if len(args) == 1 and isinstance(
                    args[0], (sp.MatrixBase, sp.MatrixExpr)):
                n, m = args[0].shape
            elif len(args) < 2:
                raise TypeError('Number of rows and columns are required')
            else:
                n = args.pop(0)
                m = args.pop(0)
            self._name = sp.MatrixSymbol(name, n, m, **kwargs)
        if not args:
            self._expr = None
        else:
            self.expr = args.pop(0)
        assert not args
        return self

    @NamedExpression.expr.setter
    def expr(self, value):
        if value is None:
            pass
        elif not isinstance(value, (sp.MatrixBase, sp.MatrixExpr)):
            raise TypeError('Expression must be a matrix')
        elif value.shape != self.shape:
            raise ValueError('Shape mismatch: {!r} vs {!r}'.format(
                value.shape, self.shape))
        self._expr = value

    @property
    def shape(self):
        return self._name.shape

    @property
    def T(self):
        expr = None if self.expr is None else self.expr.T
        return self.func(self.name.T, expr)

    @T.setter
    def T(self, value):
        self.expr = value.T

    @property
    def I(self):
        return self.func(*map(_inverse, self.args))

    @I.setter
    def I(self, value):
        self.expr = _inverse(value)

    def as_explicit(self):
        try:
            name = self.name.as_explicit()
        except AttributeError:
            name = self.name
        try:
            expr = self.expr.as_explicit()
        except AttributeError:
            expr = self.expr
        return self.func(name, expr)

    def det(self):
        return NamedExpression(sp.det(self.name), sp.det(self.expr))

    def simplify(self, *args, **kwargs):
        return self.func(self.name,
                         sp.MatrixBase.simplify(self.expr, *args, **kwargs))


def _inverse(expr):
    if expr is None:
        return None
    elif isinstance(expr, sp.MatrixBase):
        # .simplify() behaves differently on mutable and immutable matrices,
        # see https://github.com/sympy/sympy/issues/2647
        return sp.MatrixBase.simplify(expr.inv())
    elif isinstance(expr, sp.MatrixExpr):
        return expr.inverse()
    raise TypeError('Unable to invert')
