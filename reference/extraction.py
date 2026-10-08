"""Extract proof-separated branches and realize a small matrix tensor.

This finite boundary construction uses C(a,1), not the spectral 9/4 argument.
Its singleton-leg matrix factorization is explicit Gaussian elimination of the
generated coefficients; it is not a low-rank search for general tensors.
"""

from dataclasses import dataclass

from .convolution import convolution_scheme
from .fields import FiniteField
from .maps import LinearMap
from .schemes import BilinearScheme, descend
from .sectors import BRANCHES, ThreeSectorConstruction, convolution_tensor
from .separation import FiniteSeparation
from .tensor_schemes import TensorScheme


def simplify_scheme(scheme):
    """Drop zero terms and merge proportional A,B pairs, exactly over the field.

    Normalize each nonzero vector by its first nonzero coordinate and absorb
    the two scalars into C. This preserves the supplied decomposition's tensor;
    neither its rank nor its arithmetic cost is asserted to be minimal.
    """
    if not isinstance(scheme, TensorScheme):
        raise ValueError("simplification requires a tensor scheme")
    scheme.require_exact()
    f, merged = scheme.field, {}
    for a, b, c in zip(scheme.a, scheme.b, scheme.c):
        if not any(a) or not any(b) or not any(c):
            continue
        sa, sb = next(x for x in a if x), next(x for x in b if x)
        aa = tuple(f.div(x, sa) for x in a)
        bb = tuple(f.div(x, sb) for x in b)
        scale = f.mul(sa, sb)
        previous = merged.get((aa, bb), (0,)*scheme.shape[2])
        merged[aa, bb] = tuple(f.add(x, f.mul(scale, y)) for x, y in zip(previous, c))
    rows = tuple((a, b, c) for (a, b), c in merged.items() if any(c))
    result = TensorScheme(scheme.target, *(tuple(row[axis] for row in rows) for axis in range(3)))
    result.require_exact()
    return result


def factor_singleton_leg(scheme):
    """Factor the coefficient matrix of a verified scheme whose Y axis is 1.

    Reconstruct M[i,k] from the supplied rank-one families, then row-eliminate
    M to produce M=A*C. Thus this compacts actual generated coefficients and
    applies only where three-leg rank reduces to ordinary matrix rank.
    """
    if not isinstance(scheme, TensorScheme) or scheme.shape[1] != 1:
        raise ValueError("matrix factorization requires a singleton second leg")
    scheme.require_exact()
    f, (x, _, z) = scheme.field, scheme.shape
    matrix = [[0]*z for _ in range(x)]
    for a, b, c in zip(scheme.a, scheme.b, scheme.c):
        for i, value in enumerate(a):
            if value and b[0]:
                scale = f.mul(value, b[0])
                matrix[i] = [f.add(old, f.mul(scale, entry)) for old, entry in zip(matrix[i], c)]
    basis, pivots, coordinates = [], [], []
    for row in matrix:
        residual, coefficients = list(row), []
        for pivot, vector in zip(pivots, basis):
            scalar = residual[pivot]
            coefficients.append(scalar)
            residual = [f.sub(value, f.mul(scalar, entry)) for value, entry in zip(residual, vector)]
        if any(residual):
            pivot = next(i for i, value in enumerate(residual) if value)
            scalar = residual[pivot]
            pivots.append(pivot)
            basis.append(tuple(f.div(value, scalar) for value in residual))
            coefficients.append(scalar)
        coordinates.append(coefficients)
    aa = tuple(tuple(row[q] if q < len(row) else 0 for row in coordinates) for q in range(len(basis)))
    result = TensorScheme(scheme.target, aa, ((1,),)*len(basis), tuple(basis))
    result.require_exact()
    return result


def branch_maps(construction, separation, branch, dot_coordinate=0):
    """Select one branch and one diagonal entry of its auxiliary dot factor.

    Maps are in the original separated coordinate order. For the middle
    branch, these first return (X,Z_of_convolution,Y_of_convolution); the
    subsequent leg swap is essential to recover C(a,h).
    """
    if not isinstance(construction, ThreeSectorConstruction) or not isinstance(separation, FiniteSeparation):
        raise ValueError("branch extraction requires a sector construction and finite separation")
    construction.require_valid()
    if separation.branches != construction.branches:
        raise ValueError("separation must use these exact proof-generated branches")
    if branch not in BRANCHES:
        raise ValueError("unknown proof branch")
    if type(dot_coordinate) is not int or not 0 <= dot_coordinate < separation.count:
        raise ValueError("auxiliary dot coordinate is out of range")
    label, m = BRANCHES.index(branch), separation.count
    a, h, (_, y, z) = construction.a, construction.h, construction.shape
    width = max(0, a+h-1)
    if branch == "left":
        axes = (range(a), range(h), range(width))
    elif branch == "middle":
        axes = (range(a-1, -1, -1), range(h, h+width),
                range(construction.middle_z_start, construction.middle_z_start+h))
    else:
        start = construction.right_start
        axes = (range(a), range(start, start+h), range(start, start+width))
    indices = (tuple(i*m+label for i in axes[0]),
               tuple((label*y+j)*m+dot_coordinate for j in axes[1]),
               tuple((label*z+k)*m+dot_coordinate for k in axes[2]))
    return tuple(LinearMap.selection(construction.field, size, selected)
                 for size, selected in zip(separation.shape, indices))


def extract_branch(construction, separation, scheme, branch, dot_coordinate=0):
    """Recover compact C(a,h) from a certified separated scheme."""
    maps = branch_maps(construction, separation, branch, dot_coordinate)
    if not isinstance(scheme, TensorScheme) or scheme.target != separation.target:
        raise ValueError("scheme must target this exact separated tensor")
    result = scheme.restrict(*maps)
    if branch == "middle":
        result = result.permute_axes((0, 2, 1))
    if result.target != convolution_tensor(construction.field, construction.a, construction.h):
        raise ValueError("extracted branch is not the claimed convolution tensor")
    result.require_exact()
    return result


def square_from_boundary(branch_schemes):
    """Tensor three cyclic variants of C(a,1) and reindex to matrix products.

    The product axes have labels X=(i,j), Y=(k,j), Z=(i,k). Transposing
    the Y coordinate pairs realizes ordinary row-major A[i,j] B[j,k].
    Singleton factorization BEFORE tensoring prevents a cubic blowup of the
    interpolated families. All three input schemes must already be exact.
    """
    branches = tuple(branch_schemes)
    if len(branches) != 3 or any(not isinstance(s, TensorScheme) for s in branches):
        raise ValueError("three boundary convolution schemes are required")
    f, a = branches[0].field, branches[0].shape[0]
    if a < 1 or any(s.target != convolution_tensor(f, a, 1) for s in branches):
        raise ValueError("all three schemes must target the same positive C(a,1)")
    compact = tuple(factor_singleton_leg(simplify_scheme(s)) for s in branches)
    orders = ((0, 1, 2), (1, 2, 0), (2, 0, 1))
    factors = tuple(s.permute_axes(order) for s, order in zip(compact, orders))
    tensor = factors[0].tensor_product(factors[1]).tensor_product(factors[2])
    result = tensor.restrict(LinearMap.identity(f, a*a),
                             LinearMap.selection(f, a*a, (k*a+j for j in range(a) for k in range(a))),
                             LinearMap.identity(f, a*a))
    return result.to_matrix_scheme(a, a, a)


@dataclass(frozen=True)
class ProofMatrixConstruction:
    """Concrete scheme plus counts recording every finite proof stage."""

    matrix_scheme: BilinearScheme
    branch_schemes: tuple[TensorScheme, ...]
    merged_branches: tuple[TensorScheme, ...]
    factored_branches: tuple[TensorScheme, ...]
    stage_terms: tuple[tuple[str, int], ...]
    period: int
    root: int

    def descend_power(self, exponent):
        """Take the matrix power over the fixed field, then descend ONCE."""
        return descend(self.matrix_scheme.tensor_power(exponent))


def proof_matrix_pipeline(a=2, field=None):
    """Execute evaluation -> regions -> tags -> Fourier -> branch -> matrix.

    The default F_16 provides order 15 and the required nonzero nodes. This
    boundary example produces a^3 matrix terms; it establishes a working
    construction path, not the 9/4 rank bound or a speed improvement.
    """
    if type(a) is not int or a < 1:
        raise ValueError("matrix size must be a positive integer")
    field = FiniteField(2, (1, 1, 0, 0, 1)) if field is None else field
    sectors = ThreeSectorConstruction(field, a, 1)
    source = convolution_scheme(field, a, sectors.source_width)
    retained = sectors.recover_power(1, source_scheme=source)
    tagged = sectors.tag_scheme(retained)
    separation = FiniteSeparation(sectors.branches)
    separated = separation.recover_scheme(source_scheme=tagged)
    branches = tuple(extract_branch(sectors, separation, separated, branch) for branch in BRANCHES)
    merged = tuple(simplify_scheme(s) for s in branches)
    factored = tuple(factor_singleton_leg(s) for s in merged)
    matrix = square_from_boundary(factored)
    matrix.require_exact()
    return ProofMatrixConstruction(matrix, branches, merged, factored,
                                   (("convolution", source.terms), ("retained", retained.terms),
                                    ("tagged", tagged.terms), ("projected", separation.period*tagged.terms),
                                    ("separated", separated.terms), ("matrix", matrix.terms)),
                                   separation.period, separation.root)
