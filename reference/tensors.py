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

    def restrict(self, *maps):
        """Apply independent local maps: sum_ijk A[x,i] B[y,j] C[z,k] T[i,j,k]."""
        from .maps import validate_local_maps

        maps = validate_local_maps(self.field, self.shape, maps)
        shape = tuple(m.output_size for m in maps)
        coefficients = [0]*(shape[0]*shape[1]*shape[2])
        columns = tuple(m.column_supports for m in maps)
        for index, value in enumerate(self.coefficients):
            if not value:
                continue
            ij, k = divmod(index, self.shape[2])
            i, j = divmod(ij, self.shape[1])
            for (x, a), (y, b), (z, c) in product(columns[0][i], columns[1][j], columns[2][k]):
                out = (x*shape[1]+y)*shape[2]+z
                coefficients[out] = self.field.add(coefficients[out], self.field.mul(value, self.field.mul(a, self.field.mul(b, c))))
        return FiniteTensor(self.field, shape, coefficients)

    def permute_axes(self, order):
        """New axis i is old axis order[i]; this changes bilinear leg roles."""
        order = tuple(order)
        if len(order) != 3 or any(type(i) is not int for i in order) or set(order) != {0, 1, 2}:
            raise ValueError("axis order must be a permutation of 0,1,2")
        shape = tuple(self.shape[i] for i in order)

        def coefficient(*indices):
            old = [0, 0, 0]
            for axis, index in zip(order, indices):
                old[axis] = index
            return self.coefficient(*old)

        return FiniteTensor.from_function(self.field, shape, coefficient)

    @classmethod
    def direct_sum(cls, tensors):
        """Full independent blocks on ALL three axes, with concatenated coordinates."""
        tensors = _tensor_family(tensors)
        shape = tuple(sum(t.shape[i] for t in tensors) for i in range(3))
        coefficients = [0]*(shape[0]*shape[1]*shape[2])
        offsets = [0, 0, 0]
        for tensor in tensors:
            for i, j, k in product(*(range(d) for d in tensor.shape)):
                out = ((offsets[0]+i)*shape[1]+offsets[1]+j)*shape[2]+offsets[2]+k
                coefficients[out] = tensor.coefficient(i, j, k)
            offsets = [offset+d for offset, d in zip(offsets, tensor.shape)]
        return cls(tensors[0].field, shape, coefficients)

    def project_to_prime_field(self):
        """Apply the constant-coordinate linear functional pi to every coefficient.

        This is a linear projection, not a field homomorphism. For a tensor
        already defined over F_p, it recovers its base-field coefficients.
        """
        return FiniteTensor(FiniteField(self.field.p), self.shape,
                            tuple(self.field.coordinates(c)[0] for c in self.coefficients))


def _tensor_family(tensors, uniform=False):
    tensors = tuple(tensors)
    if not tensors or any(not isinstance(t, FiniteTensor) for t in tensors):
        raise ValueError("tensor family must be nonempty and contain finite tensors")
    if any(t.field != tensors[0].field for t in tensors):
        raise ValueError("tensor family must use the same field representation")
    if uniform and any(t.shape != tensors[0].shape for t in tensors):
        raise ValueError("branch tensors must have the same ambient shape")
    return tensors


def shared_first_tensor(branches):
    """Tag only Y and Z: X remains shared, matching Separation/Basic.lean."""
    branches = _tensor_family(branches, uniform=True)
    x, y, z = branches[0].shape

    def coefficient(i, j, k):
        h, j = divmod(j, y)
        h2, k = divmod(k, z)
        return branches[h].coefficient(i, j, k) if h == h2 else 0

    return FiniteTensor.from_function(branches[0].field, (x, len(branches)*y, len(branches)*z), coefficient)


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
