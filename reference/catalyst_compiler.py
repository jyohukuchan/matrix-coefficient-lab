"""Compile actual powered catalytic restrictions and coefficient absorption.

Every restriction contains three explicit local maps and is independently
checked against its full target tensor. No arrays are inferred from a cost
planner. Dense limits report incomplete construction, not nonexistence.
"""

from dataclasses import dataclass
from functools import cached_property

from .catalysts import CatalyticCertificate
from .fields import FiniteField
from .maps import LinearMap, validate_local_maps
from .schemes import BilinearScheme, descend, naive_scheme
from .scalar_summands import eliminate_scalar_summands
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, matrix_multiplication_tensor


class CompilerResourceLimit(ValueError):
    """The requested dense construction is incomplete under supplied limits."""


@dataclass(frozen=True)
class CompilerLimits:
    max_tensor_coefficients: int = 2_000_000
    max_map_entries: int = 2_000_000
    max_scheme_scalars: int = 2_000_000
    max_terms: int = 100_000

    def __post_init__(self):
        for value in (self.max_tensor_coefficients, self.max_map_entries,
                      self.max_scheme_scalars, self.max_terms):
            if type(value) is not int or value < 0:
                raise ValueError("compiler limits must be nonnegative integers")

    def tensor(self, shape):
        if shape[0]*shape[1]*shape[2] > self.max_tensor_coefficients:
            raise CompilerResourceLimit("dense tensor limit reached; construction incomplete")

    def maps(self, source, target):
        self.tensor(source)
        self.tensor(target)
        if sum(x*y for x, y in zip(source, target)) > self.max_map_entries:
            raise CompilerResourceLimit("local-map allocation limit reached; construction incomplete")

    def scheme(self, shape, terms):
        self.tensor(shape)
        if terms > self.max_terms or terms*sum(shape) > self.max_scheme_scalars:
            raise CompilerResourceLimit("coefficient-family limit reached; construction incomplete")


DEFAULT_LIMITS = CompilerLimits()


def _sum_shape(shapes):
    return tuple(sum(shape[axis] for shape in shapes) for axis in range(3))


def _direct(tensors, limits):
    tensors = tuple(tensors)
    limits.tensor(_sum_shape(tuple(t.shape for t in tensors)))
    return FiniteTensor.direct_sum(tensors)


def _identity(field, shape, limits=DEFAULT_LIMITS):
    limits.maps(shape, shape)
    return tuple(LinearMap.identity(field, size) for size in shape)


def _block_maps(field, source_shapes, target_shapes, components, limits):
    """Place independent restrictions in explicitly chosen direct-sum slots."""
    source_shapes, target_shapes = tuple(source_shapes), tuple(target_shapes)
    source, target = _sum_shape(source_shapes), _sum_shape(target_shapes)
    limits.maps(source, target)
    inputs = tuple(i for _, positions, _ in components for i in positions)
    outputs = tuple(i for _, _, positions in components for i in positions)
    if sorted(inputs) != list(range(len(source_shapes))) or sorted(outputs) != list(range(len(target_shapes))):
        raise ValueError("block restrictions must partition all input and output slots")
    result = []
    for axis in range(3):
        source_offsets, target_offsets = [], []
        for shapes, offsets in ((source_shapes, source_offsets), (target_shapes, target_offsets)):
            total = 0
            for shape in shapes:
                offsets.append(total)
                total += shape[axis]
        rows = [[0]*source[axis] for _ in range(target[axis])]
        for maps, positions_in, positions_out in components:
            columns = tuple(j for pos in positions_in
                            for j in range(source_offsets[pos], source_offsets[pos]+source_shapes[pos][axis]))
            selected = tuple(j for pos in positions_out
                             for j in range(target_offsets[pos], target_offsets[pos]+target_shapes[pos][axis]))
            mapping = maps[axis]
            if mapping.field != field or mapping.input_size != len(columns) or mapping.output_size != len(selected):
                raise ValueError("block map sizes disagree with selected slots")
            for output, row in zip(selected, mapping.rows):
                for column, value in zip(columns, row):
                    rows[output][column] = value
        result.append(LinearMap(field, source[axis], rows))
    return tuple(result)


@dataclass(frozen=True)
class CompiledRestriction:
    source: FiniteTensor
    target: FiniteTensor
    maps: tuple[LinearMap, LinearMap, LinearMap]

    def __post_init__(self):
        if self.source.field != self.target.field:
            raise ValueError("restriction tensors must use one fixed field")
        maps = validate_local_maps(self.source.field, self.source.shape, self.maps)
        if tuple(m.output_size for m in maps) != self.target.shape:
            raise ValueError("restriction output axes have the wrong sizes")
        object.__setattr__(self, "maps", maps)

    @cached_property
    def mapped(self):
        return self.source.restrict(*self.maps)

    def verify(self):
        return self.mapped == self.target

    def require_valid(self):
        if not self.verify():
            raise ValueError("compiled restriction fails full coefficient equality")


def unit_from_nonzero(tensor):
    """Select one nonzero coefficient and normalize it to an exact scalar unit."""
    index = next((i for i, value in enumerate(tensor.coefficients) if value), None)
    if index is None:
        raise ValueError("cannot extract a unit from the zero tensor")
    xy, z = divmod(index, tensor.shape[2])
    x, y = divmod(xy, tensor.shape[1])
    scalar = tensor.field.inv(tensor.coefficients[index])
    maps = tuple(LinearMap(tensor.field, size, (tuple(scale if i == chosen else 0 for i in range(size)),))
                 for size, chosen, scale in zip(tensor.shape, (x, y, z), (1, 1, scalar)))
    if tensor.restrict(*maps) != FiniteTensor.unit(tensor.field):
        raise ValueError("nonzero coefficient selector does not yield unit")
    return maps


def coordinate_catalyst(d=2, *, D=None, S=None, m=1, limits=DEFAULT_LIMITS):
    """A genuine conservative catalyst with k=d^3+m, not a useful 9/4 witness."""
    if type(d) is not int or d < 2 or type(m) is not int or m < 1:
        raise ValueError("require integer d>=2 and scalar gain m>=1")
    if S is None:
        S = FiniteTensor.unit(D.field if D is not None else FiniteField(2))
    if D is None:
        D = FiniteTensor.zero(S.field, (0, 0, 0))
    if D.field != S.field:
        raise ValueError("D and S must use the same field")
    unit = unit_from_nonzero(S)
    k, f = d**3+m, S.field
    source_shape = tuple(c+k*s for c, s in zip(D.shape, S.shape))
    target_shape = tuple(c+m+d*d*s for c, s in zip(D.shape, S.shape))
    limits.maps(source_shape, target_shape)
    limits.scheme((d*d,)*3, d**3)
    matrix = naive_scheme(f, d, d, d)
    families = matrix.a, matrix.b, matrix.c
    maps = []
    for axis, (catalyst, auxiliary) in enumerate(zip(D.shape, S.shape)):
        rows = [[0]*source_shape[axis] for _ in range(target_shape[axis])]
        for i in range(catalyst):
            rows[i][i] = 1
        for gain in range(m):
            for j, value in enumerate(unit[axis].rows[0]):
                rows[catalyst+gain][catalyst+gain*auxiliary+j] = value
        for output in range(d*d):
            for s in range(auxiliary):
                for q in range(d**3):
                    rows[catalyst+m+output*auxiliary+s][catalyst+(m+q)*auxiliary+s] = families[axis][q][output]
        maps.append(LinearMap(f, source_shape[axis], rows))
    certificate = CatalyticCertificate(d, k, m, D, S, tuple(maps))
    certificate.require_valid()
    return certificate


@dataclass(frozen=True)
class PoweredCatalyst:
    certificate: CatalyticCertificate
    exponent: int
    catalyst_scheme: TensorScheme
    matrix_size: int
    copies: int
    restriction: CompiledRestriction

    @property
    def catalyst(self):
        return self.catalyst_scheme.target


def _lift_core(certificate, previous_size, limits):
    """Tensor the core maps with M_previous and implement matrix pair reindexing."""
    f, d, u = certificate.field, certificate.d, previous_size
    p = u*u
    new_matrix = (u*d)**2
    source = tuple(c*p+certificate.k*p*s for c, s in zip(certificate.D.shape, certificate.S.shape))
    target = tuple(c*p+new_matrix*s for c, s in zip(certificate.D.shape, certificate.S.shape))
    limits.maps(source, target)
    core = certificate.strip_gain()
    maps = []
    for axis, (c, s) in enumerate(zip(certificate.D.shape, certificate.S.shape)):
        # Canonical source D*M_u + k*(M_u*S), versus old (D+k*S)*M_u.
        inputs = [(i*p+j) for i in range(c) for j in range(p)]
        inputs += [(c+h*s+z)*p+j for h in range(certificate.k) for j in range(p) for z in range(s)]
        outputs = [(i*p+j) for i in range(c) for j in range(p)]
        for q in range(new_matrix):
            i, j = divmod(q, u*d)
            old_p = (i//d)*u+j//d
            old_t = (i%d)*d+j%d
            outputs.extend((c+old_t*s+z)*p+old_p for z in range(s))
        rows = []
        for output in outputs:
            old_out, spectator = divmod(output, p)
            rows.append(tuple(core[axis].rows[old_out][old_in//p] if old_in % p == spectator else 0
                              for old_in in inputs))
        maps.append(LinearMap(f, source[axis], rows))
    return tuple(maps)


def compile_powered_catalyst(certificate, catalyst_scheme, exponent=1, *, limits=DEFAULT_LIMITS):
    """Generate actual maps D_b+M_(d^b)*S <= D_b+k^b*S, keeping S fixed."""
    if type(exponent) is not int or exponent < 1:
        raise ValueError("powered catalyst exponent must be a positive integer")
    if not isinstance(certificate, CatalyticCertificate):
        raise ValueError("an actual catalytic certificate is required")
    if not isinstance(catalyst_scheme, TensorScheme) or catalyst_scheme.target != certificate.D:
        raise ValueError("supply an exact decomposition of this catalyst D")
    # S is nonzero, so every S axis is positive. The target contains
    # M_(d^exponent) tensor S and therefore at least d^(6*exponent)
    # coefficients. Since d>=2, this lower bound rejects huge exponents
    # before computing k^exponent, d^exponent, or iterating exponent times.
    if 6*exponent >= limits.max_tensor_coefficients.bit_length():
        raise CompilerResourceLimit("dense tensor limit reached; construction incomplete")
    # Preflight the largest planned source before constructing any powers.
    copies = certificate.k**exponent
    dims = [0, 0, 0]
    terms = 0
    for i in range(exponent):
        count = certificate.k**(exponent-1-i)
        for axis in range(3):
            dims[axis] += count*certificate.D.shape[axis]*(certificate.d**(2*i))
        terms += count*catalyst_scheme.terms*certificate.d**(3*i)
    largest_source = tuple(c+copies*s for c, s in zip(dims, certificate.S.shape))
    largest_target = tuple(c+certificate.d**(2*exponent)*s for c, s in zip(dims, certificate.S.shape))
    limits.maps(largest_source, largest_target)
    limits.scheme(tuple(dims), terms)
    certificate.require_valid()
    catalyst_scheme.require_exact()
    core = CompiledRestriction(certificate.source,
                               _direct((certificate.D, matrix_multiplication_tensor(certificate.field, certificate.d,
                                        certificate.d, certificate.d).tensor_product(certificate.S)), limits),
                               certificate.strip_gain())
    core.require_valid()
    result = PoweredCatalyst(certificate, 1, catalyst_scheme, certificate.d, certificate.k, core)
    for b in range(2, exponent+1):
        f, k, S = certificate.field, certificate.k, certificate.S
        previous = result
        limits.scheme((previous.matrix_size**2,)*3, previous.matrix_size**3)
        P = naive_scheme(f, previous.matrix_size, previous.matrix_size, previous.matrix_size).as_tensor_scheme()
        DP_scheme = catalyst_scheme.tensor_product(P)
        DP, C = DP_scheme.target, previous.catalyst
        PS = P.target.tensor_product(S)
        catalyst = TensorScheme.direct_sum((DP_scheme,)+(previous.catalyst_scheme,)*k)
        source = _direct((catalyst.target,)+(S,)*(k*previous.copies), limits)
        source_shapes = (DP.shape,)+(C.shape,)*k+(S.shape,)*(k*previous.copies)
        intermediate_shapes = (DP.shape,)+(C.shape,)*k+(PS.shape,)*k
        components = [(_identity(f, DP.shape, limits), (0,), (0,))]
        for h in range(k):
            positions = (1+h,)+tuple(1+k+h*previous.copies+j for j in range(previous.copies))
            components.append((previous.restriction.maps, positions, (1+h, 1+k+h)))
        first = _block_maps(f, source_shapes, intermediate_shapes, components, limits)
        new_size = previous.matrix_size*certificate.d
        TS = matrix_multiplication_tensor(f, new_size, new_size, new_size).tensor_product(S)
        target_shapes = (DP.shape,)+(C.shape,)*k+(TS.shape,)
        components = [(_lift_core(certificate, previous.matrix_size, limits),
                       (0,)+tuple(range(1+k, 1+2*k)), (0, 1+k))]
        components += [(_identity(f, C.shape, limits), (1+h,), (1+h,)) for h in range(k)]
        second = _block_maps(f, intermediate_shapes, target_shapes, components, limits)
        maps = tuple(after.compose(before) for after, before in zip(second, first))
        target = _direct((catalyst.target, TS), limits)
        restriction = CompiledRestriction(source, target, maps)
        restriction.require_valid()
        catalyst.require_exact()
        result = PoweredCatalyst(certificate, b, catalyst, new_size, k*previous.copies, restriction)
    return result


def _scheme_from_auxiliary(scheme, S, unit, limits):
    """Promote an exact C-term decomposition of D to an actual D <= C*S map."""
    source = tuple(scheme.terms*s for s in S.shape)
    limits.maps(source, scheme.shape)
    return tuple(LinearMap(S.field, size,
                           tuple(tuple(S.field.mul(row[out], value) for row in family for value in selector.rows[0])
                                 for out in range(output)))
                 for size, output, family, selector in zip(source, scheme.shape,
                                                          (scheme.a, scheme.b, scheme.c), unit))


@dataclass(frozen=True)
class CompiledAbsorption:
    matrix_scheme: BilinearScheme
    restriction: CompiledRestriction
    input_terms: int
    auxiliary_copies: int
    catalyst_terms: int

    def descend(self):
        """Descend only this final matrix scheme over its fixed supplied field."""
        return descend(self.matrix_scheme)


def compile_absorption(powered, seed, auxiliary_scheme, *, limits=DEFAULT_LIMITS):
    """Absorb ONE catalyst repeatedly and construct the next matrix coefficient arrays."""
    if not isinstance(powered, PoweredCatalyst) or not isinstance(seed, BilinearScheme):
        raise ValueError("supply a compiled powered catalyst and a square matrix seed")
    if not seed.n == seed.m == seed.k or seed.field != powered.certificate.field:
        raise ValueError("seed must be square over the fixed coefficient field")
    S, C, f = powered.certificate.S, powered.catalyst, seed.field
    if not isinstance(auxiliary_scheme, TensorScheme) or auxiliary_scheme.target != S:
        raise ValueError("supply an exact scheme of the same auxiliary S")
    r, K, cb = seed.terms, powered.copies, powered.catalyst_scheme.terms
    N = cb+r*K
    v, u = powered.matrix_size, seed.n
    source_shape = tuple(N*s for s in S.shape)
    final_shape = ((u*v)**2,)*3
    limits.maps(source_shape, final_shape)
    limits.scheme(source_shape, N*auxiliary_scheme.terms)
    limits.scheme(final_shape, N*auxiliary_scheme.terms)
    powered.restriction.require_valid()
    powered.catalyst_scheme.require_exact()
    seed.require_exact()
    auxiliary_scheme.require_exact()
    unit = unit_from_nonzero(S)
    promoted = _scheme_from_auxiliary(powered.catalyst_scheme, S, unit, limits)
    source = _direct((S,)*N, limits)
    source_shapes = (S.shape,)*N
    KS_shape = tuple(K*s for s in S.shape)
    current_shapes = (C.shape,)+(KS_shape,)*r
    components = [(promoted, tuple(range(cb)), (0,))]
    for q in range(r):
        identity = _identity(f, KS_shape, limits)
        components.append((identity, tuple(cb+q*K+j for j in range(K)), (q+1,)))
    maps = _block_maps(f, source_shapes, current_shapes, components, limits)
    TS_shape = tuple(v*v*s for s in S.shape)
    # Replace each K*S slot in turn using the same catalyst, retaining its
    # cross-slot local maps. There is no cancellation or fresh catalyst copy.
    for q in range(r):
        following = (C.shape,)+(TS_shape,)*(q+1)+(KS_shape,)*(r-q-1)
        components = [(powered.restriction.maps, (0, q+1), (0, q+1))]
        components += [(_identity(f, shape, limits), (j,), (j,))
                       for j, shape in enumerate(current_shapes) if j not in (0, q+1)]
        step = _block_maps(f, current_shapes, following, components, limits)
        maps = tuple(after.compose(before) for after, before in zip(step, maps))
        current_shapes = following
    final_maps = []
    for axis, (family, selector, catalyst_size, s) in enumerate(zip((seed.a, seed.b, seed.c), unit, C.shape, S.shape)):
        input_size = sum(shape[axis] for shape in current_shapes)
        limits.maps(_sum_shape(current_shapes), final_shape)
        rows = []
        for output in range((u*v)**2):
            i, j = divmod(output, u*v)
            old_seed = (i//v)*u+j//v
            old_matrix = (i%v)*v+j%v
            row = [0]*input_size
            for q in range(r):
                for z, value in enumerate(selector.rows[0]):
                    row[catalyst_size+q*v*v*s+old_matrix*s+z] = f.mul(family[q][old_seed], value)
            rows.append(tuple(row))
        final_maps.append(LinearMap(f, input_size, rows))
    maps = tuple(after.compose(before) for after, before in zip(final_maps, maps))
    target = matrix_multiplication_tensor(f, u*v, u*v, u*v)
    restriction = CompiledRestriction(source, target, maps)
    restriction.require_valid()
    supplied = TensorScheme.direct_sum((auxiliary_scheme,)*N)
    result = supplied.restrict(*maps).to_matrix_scheme(u*v, u*v, u*v)
    result.require_exact()
    return CompiledAbsorption(result, restriction, r, N, cb)


@dataclass(frozen=True)
class CompiledGainAbsorption:
    """Verified original-certificate absorption retaining every scalar gain."""

    raw_gain_scheme: TensorScheme
    matrix_scheme: BilinearScheme
    restriction: CompiledRestriction
    input_terms: int
    auxiliary_copies: int
    catalyst_terms: int
    gain_removed: int

    def descend(self):
        """Descend the final matrix arrays once over the fixed coefficient field."""
        return descend(self.matrix_scheme)


def compile_gain_absorption(certificate, seed, auxiliary_scheme, catalyst_scheme,
                            *, limits=DEFAULT_LIMITS):
    """Retain scalar gains and construct C_D+(k*R_S-m)*r matrix terms.

    Use ONE actual D and r applications of the ORIGINAL positive-gain
    certificate. Decompose D directly; do not promote it to copies of S.
    The intermediate target is M_(u*d) plus r*m independent scalar units.
    Only those scalar summands are removed by exact linear elimination.
    Accept an original certificate or a verified PoweredGainCatalyst record.
    """
    if isinstance(certificate, PoweredGainCatalyst):
        certificate = certificate.certificate
    if not isinstance(certificate, CatalyticCertificate) or not isinstance(seed, BilinearScheme):
        raise ValueError("supply an actual positive-gain certificate and square matrix seed")
    if not seed.n == seed.m == seed.k or seed.field != certificate.field:
        raise ValueError("seed must be square over the fixed coefficient field")
    for scheme, target in ((auxiliary_scheme, certificate.S), (catalyst_scheme, certificate.D)):
        if not isinstance(scheme, TensorScheme) or scheme.target != target:
            raise ValueError("supply exact schemes of this certificate's D and S")
    f, D, S = certificate.field, certificate.D, certificate.S
    r, k, d, m = seed.terms, certificate.k, certificate.d, certificate.m
    copies, gains = r*k, r*m
    raw_terms = catalyst_scheme.terms+copies*auxiliary_scheme.terms
    size = seed.n*d
    source_shape = tuple(c+copies*s for c, s in zip(D.shape, S.shape))
    matrix_shape = (size*size,)*3
    target_shape = tuple(coordinate+gains for coordinate in matrix_shape)
    final_intermediate = tuple(c+gains+r*d*d*s for c, s in zip(D.shape, S.shape))
    limits.maps(source_shape, target_shape)
    limits.tensor(final_intermediate)
    limits.scheme(source_shape, raw_terms)
    limits.scheme(target_shape, raw_terms)
    limits.scheme(matrix_shape, max(0, raw_terms-gains))
    certificate.require_valid()
    seed.require_exact()
    auxiliary_scheme.require_exact()
    catalyst_scheme.require_exact()
    unit = unit_from_nonzero(S)
    scalar = FiniteTensor.unit(f)
    KS_shape = tuple(k*s for s in S.shape)
    TS_shape = tuple(d*d*s for s in S.shape)
    current_shapes = (D.shape,)+(KS_shape,)*r
    maps = None
    # Gain blocks and completed matrix blocks are spectators in subsequent
    # uses of the same D. Preserve all cross-slot entries of certificate.maps.
    for q in range(r):
        following = (D.shape,)+(scalar.shape,)*((q+1)*m)+(TS_shape,)*(q+1)+(KS_shape,)*(r-q-1)
        active_input = 1+q*m+q
        new_gains = tuple(range(1+q*m, 1+(q+1)*m))
        active_output = 1+(q+1)*m+q
        components = [(certificate.maps, (0, active_input), (0,)+new_gains+(active_output,))]
        components += [(_identity(f, scalar.shape, limits), (1+j,), (1+j,)) for j in range(q*m)]
        components += [(_identity(f, TS_shape, limits), (1+q*m+j,), (1+(q+1)*m+j,))
                       for j in range(q)]
        components += [(_identity(f, KS_shape, limits), (active_input+j,), (active_output+j,))
                       for j in range(1, r-q)]
        step = _block_maps(f, current_shapes, following, components, limits)
        maps = step if maps is None else tuple(after.compose(before) for after, before in zip(step, maps))
        current_shapes = following
    # Core-first output order is exactly the scalar-elimination API contract.
    limits.maps(_sum_shape(current_shapes), target_shape)
    final_maps = []
    for axis, (family, selector, catalyst_size, s) in enumerate(zip((seed.a, seed.b, seed.c), unit, D.shape, S.shape)):
        input_size = sum(shape[axis] for shape in current_shapes)
        rows = []
        for output in range(size*size):
            i, j = divmod(output, size)
            old_seed = (i//d)*seed.n+j//d
            old_matrix = (i%d)*d+j%d
            row = [0]*input_size
            for q in range(r):
                for z, value in enumerate(selector.rows[0]):
                    row[catalyst_size+gains+q*d*d*s+old_matrix*s+z] = f.mul(family[q][old_seed], value)
            rows.append(tuple(row))
        rows += [tuple(int(column == catalyst_size+gain) for column in range(input_size)) for gain in range(gains)]
        final_maps.append(LinearMap(f, input_size, rows))
    maps = tuple(after.compose(before) for after, before in zip(final_maps, maps))
    source = _direct((D,)+(S,)*copies, limits)
    core = matrix_multiplication_tensor(f, size, size, size)
    target = _direct((core,)+(scalar,)*gains, limits)
    restriction = CompiledRestriction(source, target, maps)
    restriction.require_valid()
    supplied = TensorScheme.direct_sum((catalyst_scheme,)+(auxiliary_scheme,)*copies)
    raw = supplied.restrict(*maps)
    raw.require_exact()
    reduced = eliminate_scalar_summands(raw, core, gains)
    matrix = reduced.to_matrix_scheme(size, size, size)
    matrix.require_exact()
    if matrix.terms != raw_terms-gains:
        raise ValueError("scalar elimination did not realize the supplied gain count")
    return CompiledGainAbsorption(raw, matrix, restriction, r, copies, catalyst_scheme.terms, gains)


@dataclass(frozen=True)
class PoweredGainCatalyst:
    """An actual powered positive-gain certificate and supplied D_b arrays."""

    original_certificate: CatalyticCertificate
    exponent: int
    certificate: CatalyticCertificate
    catalyst_scheme: TensorScheme
    restriction: CompiledRestriction

    @property
    def matrix_size(self):
        return self.certificate.d

    @property
    def copies(self):
        return self.certificate.k

    @property
    def gain(self):
        return self.certificate.m

    @property
    def catalyst(self):
        return self.certificate.D

    @property
    def source(self):
        return self.certificate.source

    @property
    def target(self):
        return self.certificate.target

    @property
    def gain_maps(self):
        return self.certificate.maps


def _lift_positive(certificate, previous_size, limits):
    """Tensor ORIGINAL positive maps with M_u, retaining m actual matrix gains."""
    f, d, u = certificate.field, certificate.d, previous_size
    p, new_matrix = u*u, (u*d)**2
    source = tuple(c*p+certificate.k*p*s for c, s in zip(certificate.D.shape, certificate.S.shape))
    target = tuple(c*p+certificate.m*p+new_matrix*s for c, s in zip(certificate.D.shape, certificate.S.shape))
    limits.maps(source, target)
    maps = []
    for axis, (c, s) in enumerate(zip(certificate.D.shape, certificate.S.shape)):
        inputs = [i*p+j for i in range(c) for j in range(p)]
        inputs += [(c+h*s+z)*p+j for h in range(certificate.k) for j in range(p) for z in range(s)]
        outputs = [i*p+j for i in range(c) for j in range(p)]
        outputs += [(c+h)*p+j for h in range(certificate.m) for j in range(p)]
        for q in range(new_matrix):
            i, j = divmod(q, u*d)
            old_p = (i//d)*u+j//d
            old_t = (i%d)*d+j%d
            outputs.extend((c+certificate.m+old_t*s+z)*p+old_p for z in range(s))
        rows = []
        for output in outputs:
            old_out, spectator = divmod(output, p)
            rows.append(tuple(certificate.maps[axis].rows[old_out][old_in//p]
                              if old_in % p == spectator else 0 for old_in in inputs))
        maps.append(LinearMap(f, source[axis], rows))
    return tuple(maps)


def compile_gain_powered_catalyst(certificate, catalyst_scheme, exponent=1,
                                 *, limits=DEFAULT_LIMITS):
    """Compile D_b+g_b*unit+M_(d^b)*S <= D_b+k^b*S with unchanged S.

    g_b=m*sum_i k^(b-1-i)d^i counts independent scalar summands. Induction
    retains k*g_previous old scalar gains and m copies of M_previous, then
    restricts each gain matrix to its u independent diagonal entries. This
    is an explicit coordinate restriction, without rank additivity or
    cancellation of D_b. The supplied D_b length stays conservative.
    """
    if type(exponent) is not int or exponent < 1:
        raise ValueError("powered gain exponent must be a positive integer")
    if not isinstance(certificate, CatalyticCertificate):
        raise ValueError("an actual positive-gain catalytic certificate is required")
    if not isinstance(catalyst_scheme, TensorScheme) or catalyst_scheme.target != certificate.D:
        raise ValueError("supply an exact decomposition of this original D")
    if 6*exponent >= limits.max_tensor_coefficients.bit_length():
        raise CompilerResourceLimit("powered gain tensor limit reached; construction incomplete")
    k, d, m, f, S = certificate.k, certificate.d, certificate.m, certificate.field, certificate.S
    copies = k**exponent
    gain = m*sum(k**(exponent-1-i)*d**i for i in range(exponent))
    dims = tuple(sum(k**(exponent-1-i)*certificate.D.shape[axis]*d**(2*i)
                     for i in range(exponent)) for axis in range(3))
    terms = catalyst_scheme.terms*sum(k**(exponent-1-i)*d**(3*i) for i in range(exponent))
    source_shape = tuple(c+copies*s for c, s in zip(dims, S.shape))
    target_shape = tuple(c+gain+d**(2*exponent)*s for c, s in zip(dims, S.shape))
    limits.maps(source_shape, target_shape)
    limits.scheme(dims, terms)
    certificate.require_valid()
    catalyst_scheme.require_exact()
    restriction = CompiledRestriction(certificate.source, certificate.target, certificate.maps)
    restriction.require_valid()
    result = PoweredGainCatalyst(certificate, 1, certificate, catalyst_scheme, restriction)
    scalar = FiniteTensor.unit(f)
    for b in range(2, exponent+1):
        previous = result
        u, old_copies = previous.matrix_size, previous.copies
        limits.scheme((u*u,)*3, u**3)
        P = naive_scheme(f, u, u, u).as_tensor_scheme()
        DP_scheme = catalyst_scheme.tensor_product(P)
        DP, C = DP_scheme.target, previous.catalyst
        PS = P.target.tensor_product(S)
        catalyst = TensorScheme.direct_sum((DP_scheme,)+(previous.catalyst_scheme,)*k)
        new_size, old_gains = u*d, k*previous.gain
        source_shapes = (DP.shape,)+(C.shape,)*k+(S.shape,)*(k*old_copies)
        first_shapes = (DP.shape,)+(C.shape,)*k+(scalar.shape,)*old_gains+(PS.shape,)*k
        components = [(_identity(f, DP.shape, limits), (0,), (0,))]
        for h in range(k):
            positions_in = (1+h,)+tuple(1+k+h*old_copies+j for j in range(old_copies))
            gains_out = tuple(1+k+h*previous.gain+j for j in range(previous.gain))
            positions_out = (1+h,)+gains_out+(1+k+old_gains+h,)
            components.append((previous.certificate.maps, positions_in, positions_out))
        first = _block_maps(f, source_shapes, first_shapes, components, limits)
        TS_shape = tuple(new_size*new_size*s for s in S.shape)
        second_shapes = (DP.shape,)+(C.shape,)*k+(scalar.shape,)*old_gains+(P.shape,)*m+(TS_shape,)
        components = [(_lift_positive(certificate, u, limits),
                       (0,)+tuple(range(1+k+old_gains, 1+2*k+old_gains)),
                       (0,)+tuple(range(1+k+old_gains, 2+k+old_gains+m)))]
        components += [(_identity(f, shape, limits), (pos,), (pos,))
                       for pos, shape in enumerate(first_shapes[:1+k+old_gains]) if pos]
        second = _block_maps(f, first_shapes, second_shapes, components, limits)
        new_gain = old_gains+m*u
        third_shapes = (DP.shape,)+(C.shape,)*k+(scalar.shape,)*new_gain+(TS_shape,)
        components = [(_identity(f, shape, limits), (pos,), (pos,))
                      for pos, shape in enumerate(second_shapes[:1+k+old_gains])]
        limits.maps(P.shape, (u, u, u))
        diagonal = tuple(LinearMap.selection(f, u*u, (i*u+i for i in range(u))) for _ in range(3))
        for h in range(m):
            out = tuple(1+k+old_gains+h*u+j for j in range(u))
            components.append((diagonal, (1+k+old_gains+h,), out))
        components.append((_identity(f, TS_shape, limits), (1+k+old_gains+m,), (1+k+new_gain,)))
        third = _block_maps(f, second_shapes, third_shapes, components, limits)
        maps = tuple(last.compose(middle).compose(start) for start, middle, last in zip(first, second, third))
        # Constructing CatalyticCertificate supplies exact target/source
        # representatives and verifies all mixed direct-sum coefficients.
        powered = CatalyticCertificate(new_size, k*old_copies, new_gain, catalyst.target, S, maps)
        powered.require_valid()
        catalyst.require_exact()
        restriction = CompiledRestriction(powered.source, powered.target, maps)
        restriction.require_valid()
        result = PoweredGainCatalyst(certificate, b, powered, catalyst, restriction)
    return result
