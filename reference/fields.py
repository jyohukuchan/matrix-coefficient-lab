"""Exact small finite fields F_p[X]/(f), with no external dependencies.

Elements are canonical integers encoding power-basis coordinates in base p.
Use field.embed(k) to interpret an ordinary integer as a field scalar;
in an extension, the encoding p denotes X, NOT the scalar p (which is zero).
"""

from dataclasses import dataclass
from math import isqrt


def _trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def _sub(a, b, p):
    return _trim([((a[i] if i < len(a) else 0)
                   - (b[i] if i < len(b) else 0)) % p
                  for i in range(max(len(a), len(b)))])


def _mul(a, b, p):
    if not a or not b:
        return []
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] = (c[i + j] + x * y) % p
    return _trim(c)


def _rem(a, b, p):
    a = _trim(a)
    b = _trim(b)
    if not b:
        raise ZeroDivisionError("zero polynomial")
    inverse = pow(b[-1], -1, p)
    while len(a) >= len(b):
        shift = len(a) - len(b)
        scale = a[-1] * inverse % p
        for j, x in enumerate(b):
            a[shift + j] = (a[shift + j] - scale * x) % p
        a = _trim(a)
    return a


def _powmod(a, n, modulus, p):
    result = [1]
    a = _rem(a, modulus, p)
    while n:
        if n & 1:
            result = _rem(_mul(result, a, p), modulus, p)
        a = _rem(_mul(a, a, p), modulus, p)
        n //= 2
    return result


def _gcd(a, b, p):
    while b:
        a, b = b, _rem(a, b, p)
    return a


def _irreducible(f, p):
    # Every reducible degree-d polynomial has a factor of degree <= d/2.
    # gcd(f, X^(p^i)-X) detects factors whose degree divides i.
    degree = len(f) - 1
    x = _rem([0, 1], f, p)
    frobenius = x
    for i in range(1, degree + 1):
        frobenius = _powmod(frobenius, p, f, p)
        difference = _sub(frobenius, x, p)
        if i <= degree // 2 and len(_gcd(f, difference, p)) > 1:
            return False
    return frobenius == x


@dataclass(frozen=True)
class FiniteField:
    """A validated finite field, intended for small reference computations.

    modulus is monic with coefficients ordered from constant to leading.
    The default (0, 1) gives F_p itself. No algebraic closure is represented.
    """

    p: int
    modulus: tuple[int, ...] = (0, 1)

    def __post_init__(self):
        if type(self.p) is not int or self.p < 2:
            raise ValueError("characteristic must be prime")
        if any(self.p % k == 0 for k in range(2, isqrt(self.p) + 1)):
            raise ValueError("characteristic must be prime")
        if any(type(c) is not int for c in self.modulus):
            raise ValueError("modulus coefficients must be integers")
        f = tuple(c % self.p for c in self.modulus)
        if len(f) < 2 or f[-1] != 1:
            raise ValueError("modulus must be monic of positive degree")
        if not _irreducible(list(f), self.p):
            raise ValueError("modulus is reducible over the prime field")
        object.__setattr__(self, "modulus", f)

    @property
    def degree(self):
        return len(self.modulus) - 1

    @property
    def order(self):
        return self.p ** self.degree

    def check(self, a):
        if type(a) is not int or not 0 <= a < self.order:
            raise ValueError(f"element must be a canonical encoding in [0, {self.order})")
        return a

    def embed(self, integer):
        if type(integer) is not int:
            raise ValueError("scalar must be an integer")
        return integer % self.p

    def coordinates(self, a):
        self.check(a)
        result = []
        for _ in range(self.degree):
            result.append(a % self.p)
            a //= self.p
        return tuple(result)

    def from_coordinates(self, coordinates):
        coordinates = tuple(coordinates)
        if len(coordinates) != self.degree:
            raise ValueError("wrong number of power-basis coordinates")
        if any(type(c) is not int for c in coordinates):
            raise ValueError("coordinates must be integers")
        return sum((c % self.p) * self.p ** i for i, c in enumerate(coordinates))

    def add(self, a, b):
        return self.from_coordinates((x + y) % self.p for x, y in
                                     zip(self.coordinates(a), self.coordinates(b)))

    def neg(self, a):
        return self.from_coordinates(-x for x in self.coordinates(a))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        product = _mul(self.coordinates(a), self.coordinates(b), self.p)
        reduced = _rem(product, list(self.modulus), self.p)
        return self.from_coordinates(reduced + [0] * (self.degree - len(reduced)))

    def pow(self, a, n):
        self.check(a)
        if type(n) is not int:
            raise ValueError("exponent must be an integer")
        if n < 0:
            return self.pow(self.inv(a), -n)
        result = 1
        while n:
            if n & 1:
                result = self.mul(result, a)
            a = self.mul(a, a)
            n //= 2
        return result

    def inv(self, a):
        self.check(a)
        if a == 0:
            raise ZeroDivisionError("zero has no inverse")
        return self.pow(a, self.order - 2)

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def sum(self, elements):
        result = 0
        for a in elements:
            result = self.add(result, a)
        return result

    def dot(self, left, right):
        left, right = tuple(left), tuple(right)
        if len(left) != len(right):
            raise ValueError("dot-product lengths differ")
        return self.sum(self.mul(a, b) for a, b in zip(left, right))
