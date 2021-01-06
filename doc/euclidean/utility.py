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

    @classmethod
    def solve(cls, expr, name):
        result = cls(name)
        result.expr = sp.solve(expr, result.name)[0]
        return result

    def with_name(self, name):
        """Return a copy of the expression, but with a new name."""
        return self.func(name, self.expr)

    def subs_symbols(self, *args, **kwargs):
        substitutions = []
        for arg in args:
            if not isinstance(arg, NamedExpression):
                raise TypeError('args must be of type NamedExpression')
            substitutions.append((arg.lhs, arg.rhs))
        return self.subs(substitutions, **kwargs)

    def evaluated_at(self, old, new, symbols=()):
        new_expr = self.expr.subs(old, new)
        for symbol in symbols:
            new_expr = new_expr.subs(
                symbol.name, symbol.evaluated_at(old, new).name)
        return self.func(
            r'\left.{' + sp.latex(self.name) + r'}\right\rvert_{' +
            sp.latex(old) + '=' + sp.latex(new) + '}',
            new_expr
        )

    def pull_out(self, expr):
        # NB: This ignores the subclass and just returns a NamedExpression:
        return NamedExpression(
            self.name,
            sp.UnevaluatedExpr(expr) *
            sp.UnevaluatedExpr(sp.simplify(self.expr / expr)))

    def diff(self, *args, **kwargs):
        return self.func(
            sp.Derivative(self.name, *args, **kwargs),
            self.expr.diff(*args, **kwargs))

    def _call_function(self, func, *args, **kwargs):
        return self.func(
            func(self.name, *args, **kwargs),
            func(self.expr, *args, **kwargs))

    def simplify(self, *args, **kwargs):
        return self._call_function(sp.simplify, *args, **kwargs)

    def factor(self, *args, **kwargs):
        return self._call_function(sp.factor, *args, **kwargs)

    def expand(self, *args, **kwargs):
        return self._call_function(sp.expand, *args, **kwargs)

    def _call_method(self, name, *args, **kwargs):
        return self.func(
            getattr(self.name, name)(*args, **kwargs),
            getattr(self.expr, name)(*args, **kwargs))

    def doit(self, *args, **kwargs):
        return self._call_method('doit', *args, **kwargs)

    def subs(self, *args, **kwargs):
        return self._call_method('subs', *args, **kwargs)


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

    def factor(self, *args, **kwargs):
        return self.func(
            sp.factor(self.name, *args, **kwargs),
            # https://github.com/sympy/sympy/issues/19331
            self.expr.applyfunc(lambda elem: sp.factor(elem, *args, **kwargs)))

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
