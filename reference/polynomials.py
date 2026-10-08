"""Small immutable polynomials over an explicitly represented finite field."""

from dataclasses import dataclass

from .constructions import polynomial_eval
from .fields import FiniteField


@dataclass(frozen=True)
class Polynomial:
    """Coefficients in ascending degree, using canonical field encodings.

    Trailing zeros are removed; the zero polynomial has coefficients () and
    degree -1. Equality is formal coefficient equality, not equality of the
    functions induced on the finite field.
    """

    field: FiniteField
    coefficients: tuple[int, ...] = ()

    def __post_init__(self):
        coefficients = list(self.coefficients)
        for coefficient in coefficients:
            self.field.check(coefficient)
        while coefficients and coefficients[-1] == 0:
            coefficients.pop()
        object.__setattr__(self, "coefficients", tuple(coefficients))

    @classmethod
    def constant(cls, field, value):
        return cls(field, (value,))

    @classmethod
    def monomial(cls, field, degree, coefficient=1):
        if type(degree) is not int or degree < 0:
            raise ValueError("monomial degree must be a nonnegative integer")
        field.check(coefficient)
        if coefficient == 0:
            return cls(field)
        return cls(field, (0,) * degree + (coefficient,))

    @property
    def degree(self):
        return len(self.coefficients) - 1

    def coefficient(self, degree):
        if type(degree) is not int or degree < 0:
            raise ValueError("coefficient degree must be a nonnegative integer")
        return self.coefficients[degree] if degree < len(self.coefficients) else 0

    def evaluate(self, node):
        return polynomial_eval(self.field, self.coefficients, node)

    def _check_other(self, other):
        if not isinstance(other, Polynomial):
            raise TypeError("polynomial arithmetic requires another Polynomial")
        if self.field != other.field:
            raise ValueError("polynomials must use the same field representation")

    def __add__(self, other):
        self._check_other(other)
        return Polynomial(self.field, tuple(self.field.add(self.coefficient(i), other.coefficient(i))
                                           for i in range(max(self.degree, other.degree) + 1)))

    def __neg__(self):
        return Polynomial(self.field, tuple(self.field.neg(c) for c in self.coefficients))

    def __sub__(self, other):
        self._check_other(other)
        return self + (-other)

    def __mul__(self, other):
        self._check_other(other)
        if not self.coefficients or not other.coefficients:
            return Polynomial(self.field)
        result = [0] * (self.degree + other.degree + 1)
        for i, a in enumerate(self.coefficients):
            if a == 0:
                continue
            for j, b in enumerate(other.coefficients):
                if b != 0:
                    result[i + j] = self.field.add(result[i + j], self.field.mul(a, b))
        return Polynomial(self.field, tuple(result))

    def __pow__(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("polynomial exponent must be a nonnegative integer")
        result = Polynomial.constant(self.field, 1)
        base = self
        while exponent:
            if exponent & 1:
                result = result * base
            base = base * base
            exponent //= 2
        return result
