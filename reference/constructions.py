"""Finite Fourier filtering and nonzero-node coefficient recovery.

Counterparts of AuxiliarySeparation/Separation/Fourier.lean and
AuxiliarySeparation/Polynomial/Interpolation.lean, for explicit finite fields.
"""


def separation_period(field, blocks):
    if type(blocks) is not int or blocks < 0:
        raise ValueError("block count must be a nonnegative integer")
    return 5 * blocks + int(field.embed(5 * blocks) == 0)


def _prime_divisors(n):
    result = []
    divisor = 2
    while divisor * divisor <= n:
        if n % divisor == 0:
            result.append(divisor)
            while n % divisor == 0:
                n //= divisor
        divisor += 1
    if n > 1:
        result.append(n)
    return result


def _check_period(period):
    if type(period) is not int or period < 1:
        raise ValueError("period must be a positive integer")


def is_primitive_root(field, root, period):
    _check_period(period)
    field.check(root)
    return (root != 0 and field.pow(root, period) == 1
            and all(field.pow(root, period // divisor) != 1
                    for divisor in _prime_divisors(period)))


def primitive_root(field, period):
    """Brute-force a root in the SUPPLIED field; no closure/extension search."""
    _check_period(period)
    if (field.order - 1) % period:
        raise ValueError("period does not divide field.order - 1; supply a suitable extension field")
    for root in range(1, field.order):
        if is_primitive_root(field, root, period):
            return root
    raise ValueError("no primitive root found")


def fourier_filter(field, period, root, exponent):
    """L^-1 sum_{r=0}^{L-1} root^(r*e), including negative e."""
    _check_period(period)
    if type(exponent) is not int:
        raise ValueError("Fourier exponent must be an integer")
    if field.embed(period) == 0:
        raise ValueError("period is zero in the field")
    if not is_primitive_root(field, root, period):
        raise ValueError("root does not have the required exact order")
    return field.div(field.sum(field.pow(root, r * exponent) for r in range(period)),
                     field.embed(period))


def nonzero_nodes(field, count):
    if type(count) is not int or not 1 <= count < field.order:
        raise ValueError("need between 1 and field.order - 1 distinct nonzero nodes")
    return tuple(range(1, count + 1))


def constant_weights(field, nodes):
    """Lagrange basis evaluated at zero: prod_{j!=i} -t_j/(t_i-t_j)."""
    nodes = tuple(nodes)
    for node in nodes:
        field.check(node)
    if not nodes or 0 in nodes or len(set(nodes)) != len(nodes):
        raise ValueError("nodes must be nonempty, distinct, and nonzero")
    weights = []
    for i, node in enumerate(nodes):
        weight = 1
        for j, other in enumerate(nodes):
            if i != j:
                weight = field.mul(weight, field.div(field.neg(other), field.sub(node, other)))
        weights.append(weight)
    return tuple(weights)


def polynomial_eval(field, coefficients, node):
    field.check(node)
    result = 0
    for coefficient in reversed(tuple(coefficients)):
        result = field.add(field.mul(result, node), coefficient)
    return result


def recover_leading_coefficient(field, coefficients, leading, nodes):
    """Recover [X^leading]p when p=X^leading*q and deg(q)<len(nodes).

    Checks the degree and vanishing hypotheses; evaluates p ONLY at nonzero
    nodes and divides out their leading power. The finite field must contain
    enough distinct nonzero nodes; the proof uses an infinite algebraic closure.
    """
    if type(leading) is not int or leading < 0:
        raise ValueError("leading degree must be a nonnegative integer")
    coefficients = tuple(coefficients)
    for coefficient in coefficients:
        field.check(coefficient)
    if any(coefficients[:leading]):
        raise ValueError("coefficients below the leading degree must vanish")
    normalized = list(coefficients[leading:])
    while normalized and normalized[-1] == 0:
        normalized.pop()
    nodes = tuple(nodes)
    weights = constant_weights(field, nodes)
    if len(normalized) > len(nodes):
        raise ValueError("not enough nodes for the normalized polynomial degree")
    return field.sum(field.mul(weight, field.div(polynomial_eval(field, coefficients, node),
                                                field.pow(node, leading)))
                     for node, weight in zip(nodes, weights))
