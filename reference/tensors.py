"""Explicit finite three-leg tensors over represented finite fields."""

from dataclasses import dataclass
from itertools import product

from .fields import FiniteField


def _shape(shape):
    shape = tuple(shape)
    if len(shape) != 3 or any(type(d) is not int or d < 0 for d in shape):
        raise ValueError("tensor shape must contain three nonnegative integers")
    return shape


def _vector(field, vector, size):
    vector = tuple(vector)
    if len(vector) != size:
        raise ValueError(f"expected a vector of length {size}")
    return tuple(field.check(x) for x in vector)


@dataclass(frozen=True)
class FiniteTensor:
    """T[x,y,z], stored in lexicographic order with z varying fastest.

    The axes are independent coordinate spaces, NOT matrix dimensions. Any
    finite sets can be represented by enumerating them as 0,...,dimension-1.
    Zero-length axes and the zero tensor are supported. Small dense storage is
    intentional; this is not an efficient representation for large powers.
    """

    field: FiniteField
    shape: tuple[int, int, int]
    coefficients: tuple[int, ...]

    def __post_init__(self):
        shape = _shape(self.shape)
        coefficients = tuple(self.coefficients)
        if len(coefficients) != shape[0]*shape[1]*shape[2]:
            raise ValueError("coefficient count does not match tensor shape")
        for coefficient in coefficients:
            self.field.check(coefficient)
        object.__setattr__(self, "shape", shape)
        object.__setattr__(self, "coefficients", coefficients)

    @classmethod
    def from_function(cls, field, shape, coefficient):
        shape = _shape(shape)
        return cls(field, shape, tuple(coefficient(i, j, k)
                                       for i, j, k in product(*(range(d) for d in shape))))

    @classmethod
    def zero(cls, field, shape):
        shape = _shape(shape)
        return cls(field, shape, (0,) * (shape[0]*shape[1]*shape[2]))

    @classmethod
    def unit(cls, field):
        return cls(field, (1, 1, 1), (1,))

    def coefficient(self, i, j, k):
        for index, size in zip((i, j, k), self.shape):
            if type(index) is not int or not 0 <= index < size:
                raise ValueError("tensor coordinate is out of range")
        return self.coefficients[(i*self.shape[1] + j)*self.shape[2] + k]

    def contract(self, left, right):
        """Directly compute output[z] = sum_{x,y} T[x,y,z]*left[x]*right[y]."""
        f = self.field
        x, y, z = self.shape
        left, right = _vector(f, left, x), _vector(f, right, y)
        return tuple(f.sum(f.mul(f.mul(self.coefficient(i, j, k), left[i]), right[j])
                           for i in range(x) for j in range(y)) for k in range(z))

    def tensor_product(self, other):
        """Pair each axis separately: (i1,i2) is encoded i1*dim2+i2."""
        if self.field != other.field:
            raise ValueError("tensor factors must use the same field representation")
        shape = tuple(a*b for a, b in zip(self.shape, other.shape))

        def coefficient(i, j, k):
            first, second = zip(*(divmod(index, size) for index, size in zip((i, j, k), other.shape)))
            return self.field.mul(self.coefficient(*first), other.coefficient(*second))

        return FiniteTensor.from_function(self.field, shape, coefficient)

    def tensor_power(self, exponent):
        if type(exponent) is not int or exponent < 0:
            raise ValueError("tensor exponent must be a nonnegative integer")
        result = FiniteTensor.unit(self.field)
        for _ in range(exponent):
            result = result.tensor_product(self)
        return result

    def project_to_prime_field(self):
        """Apply the constant-coordinate linear functional pi to every coefficient.

        This is a linear projection, not a field homomorphism. For a tensor
        already defined over F_p, it recovers its base-field coefficients.
        """
        return FiniteTensor(FiniteField(self.field.p), self.shape,
                            tuple(self.field.coordinates(c)[0] for c in self.coefficients))


def matrix_multiplication_tensor(field, n, m, k):
    """The matrix tensor, with row-major input and output vector coordinates."""
    if any(type(d) is not int or d < 1 for d in (n, m, k)):
        raise ValueError("matrix dimensions must be positive integers")

    def coefficient(x, y, z):
        i, j = divmod(x, m)
        j2, h = divmod(y, k)
        i2, h2 = divmod(z, k)
        return int(i == i2 and j == j2 and h == h2)

    return FiniteTensor.from_function(field, (n*m, m*k, n*k), coefficient)
