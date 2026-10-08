"""Independent local linear substitutions in [output][input] order."""

from dataclasses import dataclass
from functools import cached_property

from .fields import FiniteField
from .polynomials import Polynomial


@dataclass(frozen=True)
class LinearMap:
    field: FiniteField
    input_size: int
    rows: tuple[tuple[int, ...], ...]

    def __post_init__(self):
        if not isinstance(self.field, FiniteField):
            raise ValueError("map must use an explicitly represented finite field")
        if type(self.input_size) is not int or self.input_size < 0:
            raise ValueError("map input size must be a nonnegative integer")
        rows = tuple(tuple(row) for row in self.rows)
        if any(len(row) != self.input_size for row in rows):
            raise ValueError("map rows have the wrong input size")
        for row in rows:
            for value in row:
                self.field.check(value)
        object.__setattr__(self, "rows", rows)

    @property
    def output_size(self):
        return len(self.rows)

    @cached_property
    def row_supports(self):
        return tuple(tuple((i, value) for i, value in enumerate(row) if value) for row in self.rows)

    @cached_property
    def column_supports(self):
        columns = [[] for _ in range(self.input_size)]
        for out, row in enumerate(self.row_supports):
            for i, value in row:
                columns[i].append((out, value))
        return tuple(tuple(column) for column in columns)

    @classmethod
    def selection(cls, field, input_size, indices):
        """Pull back coordinates, allowing selection, repetition, and reordering."""
        if type(input_size) is not int or input_size < 0:
            raise ValueError("map input size must be a nonnegative integer")
        indices = tuple(indices)
        if any(type(i) is not int or not 0 <= i < input_size for i in indices):
            raise ValueError("selected input coordinate is out of range")
        return cls(field, input_size, tuple(tuple(int(j == i) for j in range(input_size)) for i in indices))

    @classmethod
    def identity(cls, field, size):
        if type(size) is not int or size < 0:
            raise ValueError("map size must be a nonnegative integer")
        return cls.selection(field, size, range(size))

    def apply(self, vector):
        vector = tuple(vector)
        if len(vector) != self.input_size:
            raise ValueError("map vector has the wrong input size")
        for value in vector:
            self.field.check(value)
        return tuple(self.field.sum(self.field.mul(value, vector[i]) for i, value in row)
                     for row in self.row_supports)

    def apply_polynomials(self, vector):
        vector = tuple(vector)
        if len(vector) != self.input_size or any(not isinstance(p, Polynomial) or p.field != self.field for p in vector):
            raise ValueError("map polynomial vector has the wrong size or field")
        result = []
        for row in self.row_supports:
            value = Polynomial(self.field)
            for i, scalar in row:
                value = value + Polynomial.constant(self.field, scalar)*vector[i]
            result.append(value)
        return tuple(result)

    def compose(self, previous):
        """Return self after previous, so apply(v) equals self(previous(v))."""
        if not isinstance(previous, LinearMap) or self.field != previous.field or self.input_size != previous.output_size:
            raise ValueError("map composition has incompatible sizes or fields")
        return LinearMap(self.field, previous.input_size,
                         tuple(tuple(self.field.sum(self.field.mul(value, previous.rows[i][j]) for i, value in row)
                                     for j in range(previous.input_size)) for row in self.row_supports))


def validate_local_maps(field, shape, maps):
    maps = tuple(maps)
    if len(maps) != 3 or any(not isinstance(m, LinearMap) or m.field != field or m.input_size != size
                             for m, size in zip(maps, shape)):
        raise ValueError("three local maps must match the tensor's field and input axes")
    return maps
