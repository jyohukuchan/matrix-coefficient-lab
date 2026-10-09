"""Checked construction DAGs with bounded expansion and coefficient queries.

Leaves and supplied restriction identities are checked coefficient by
coefficient. Product, direct sum, permutation, restriction and fixed-field
projection are exact construction rules. No numeric rank oracle is a node.
This Python checker is not a Lean proof or a practical coefficient discovery
algorithm. Scalar-gain elimination is still handled by the dense compiler.
"""

from dataclasses import dataclass
from itertools import product
from math import prod

from .fields import FiniteField
from .maps import LinearMap
from .tensor_schemes import TensorScheme
from .tensors import FiniteTensor, _shape


class GraphResourceLimit(ValueError):
    """The declared graph/query/expansion budget was exceeded."""


@dataclass(frozen=True)
class GraphLimits:
    max_nodes: int = 10_000
    max_query_work: int = 1_000_000
    max_tensor_entries: int = 2_000_000
    max_scheme_scalars: int = 2_000_000
    max_terms: int = 100_000
    max_map_entries: int = 2_000_000

    def __post_init__(self):
        if any(type(v) is not int or v < 0 for v in vars(self).values()):
            raise ValueError("graph limits must be nonnegative integers")


@dataclass(frozen=True)
class ConstructionNode:
    op: str
    args: tuple[int, ...]
    data: tuple
    kind: str
    field: FiniteField
    shape: tuple[int, ...]
    terms: int | None = None
    target: int | None = None
    base_valued: bool = False


def _freeze(value):
    if isinstance(value, (tuple, list)):
        return tuple(_freeze(x) for x in value)
    if not isinstance(value, (int, str)):
        raise ValueError("node data must contain only integer/string sequences")
    return value


def _natural(value):
    if type(value) is not int or value < 0:
        raise ValueError("expected a nonnegative integer")
    return value


class _Work:
    def __init__(self, cap):
        self.remaining = cap

    def use(self, count=1):
        self.remaining -= count
        if self.remaining < 0:
            raise GraphResourceLimit("coefficient-query work cap reached; query incomplete")


class CertificateGraph:
    """Append-only, hash-consed typed nodes; saved metadata is never trusted."""

    def __init__(self, field, limits=GraphLimits()):
        if not isinstance(field, FiniteField) or not isinstance(limits, GraphLimits):
            raise ValueError("supply a represented field and graph limits")
        self.field, self.limits = field, limits
        self._nodes, self._interned = [], {}

    @property
    def nodes(self):
        return tuple(self._nodes)

    def node(self, index):
        if type(index) is not int or not 0 <= index < len(self._nodes):
            raise ValueError("node must reference an earlier graph entry")
        return self._nodes[index]

    def _cap(self, count, name):
        if count > getattr(self.limits, name):
            raise GraphResourceLimit(f"{name} cap reached; expansion incomplete")

    def _add(self, op, args=(), data=()):
        if not isinstance(op, str):
            raise ValueError("construction rule must be a string")
        args, data = tuple(args), _freeze(data)
        if not isinstance(data, tuple):
            raise ValueError("node data must be a sequence")
        children = tuple(self.node(i) for i in args)
        key = op, args, data
        if key in self._interned:
            return self._interned[key]
        self._cap(len(self._nodes)+1, "max_nodes")
        field, terms, target, base = self.field, None, None, False

        def expect(kind, count=None):
            if count is not None and len(children) != count:
                raise ValueError("wrong node arity")
            if not children or any(c.kind != kind for c in children):
                raise ValueError(f"operation requires {kind} children")
            if any(c.field != children[0].field for c in children):
                raise ValueError("operation mixes field representations")
            return children[0].field

        def no_children():
            if children:
                raise ValueError("leaf cannot have children")

        if op == "tensor_literal":
            no_children()
            shape, coefficients = data
            shape = _shape(shape)
            self._cap(prod(shape), "max_tensor_entries")
            literal = FiniteTensor(field, shape, coefficients)
            kind, base = "tensor", all(x < field.p for x in literal.coefficients)
        elif op == "tensor_matrix":
            no_children()
            if len(data) != 3 or any(type(x) is not int or x < 1 for x in data):
                raise ValueError("matrix dimensions must be positive integers")
            n, m, k = data
            kind, shape, base = "tensor", (n*m, m*k, n*k), True
        elif op in ("tensor_product", "tensor_sum", "tensor_repeat", "tensor_permute"):
            field = expect("tensor", 2 if op == "tensor_product" else
                           1 if op in ("tensor_repeat", "tensor_permute") else None)
            kind, base = "tensor", all(c.base_valued for c in children)
            if op == "tensor_product":
                if data:
                    raise ValueError("product takes no data")
                shape = tuple(a*b for a, b in zip(children[0].shape, children[1].shape))
            elif op == "tensor_sum":
                if data:
                    raise ValueError("sum takes no data")
                shape = tuple(sum(c.shape[a] for c in children) for a in range(3))
            elif op == "tensor_repeat":
                count, = data
                shape = tuple(_natural(count)*x for x in children[0].shape)
            else:
                if len(data) != 3 or any(type(x) is not int for x in data) or set(data) != {0, 1, 2}:
                    raise ValueError("axis permutation must contain 0,1,2")
                shape = tuple(children[0].shape[a] for a in data)
        elif op in ("tensor_restrict", "tensor_witness"):
            offset = 2 if op == "tensor_witness" else 1
            if len(children) != offset+3 or any(c.kind != "tensor" for c in children[:offset]) or any(c.kind != "map" for c in children[offset:]):
                raise ValueError("restriction requires tensors and three maps")
            source, maps = children[0], children[offset:]
            field = source.field
            if any(c.field != field for c in children) or any(m.shape[0] != x for m, x in zip(maps, source.shape)):
                raise ValueError("restriction fields or input axes disagree")
            kind, shape = "tensor", tuple(m.shape[1] for m in maps)
            base = source.base_valued and all(m.base_valued for m in maps)
            if op == "tensor_witness":
                label, = data
                if not isinstance(label, str) or not label:
                    raise ValueError("witness label must be nonempty")
                if children[1].shape != shape:
                    raise ValueError("witness target shape disagrees")
                self._cap(prod(shape), "max_tensor_entries")
                actual = self.materialize_tensor(args[0]).restrict(*(self.materialize_map(i) for i in args[2:]))
                if actual != self.materialize_tensor(args[1]):
                    raise ValueError("restriction witness fails exact coefficient equality")
                base = children[1].base_valued
            elif data:
                raise ValueError("restriction takes no data")
        elif op == "tensor_prime":
            field = expect("tensor", 1)
            if data or not children[0].base_valued:
                raise ValueError("unchanged descent requires a certified base-valued target")
            field, kind, shape, base = FiniteField(field.p), "tensor", children[0].shape, True
        elif op in ("map_literal", "map_identity", "map_selection"):
            no_children()
            kind = "map"
            if op == "map_literal":
                width, rows = data
                self._cap(_natural(width)*len(rows), "max_map_entries")
                mapping = LinearMap(field, width, rows)
                shape, base = (width, mapping.output_size), all(v < field.p for row in rows for v in row)
            elif op == "map_identity":
                width, = data
                shape, base = (_natural(width), width), True
            else:
                width, indices = data
                _natural(width)
                if any(type(i) is not int or not 0 <= i < width for i in indices):
                    raise ValueError("selected map coordinate out of range")
                shape, base = (width, len(indices)), True
        elif op in ("map_product", "map_compose"):
            field = expect("map", 2)
            if data:
                raise ValueError("map operation takes no data")
            kind, base = "map", all(c.base_valued for c in children)
            a, b = children
            if op == "map_product":
                shape = tuple(x*y for x, y in zip(a.shape, b.shape))
            else:
                if a.shape[0] != b.shape[1]:
                    raise ValueError("map composition dimensions disagree")
                shape = b.shape[0], a.shape[1]
        elif op == "scheme_literal":
            field = expect("tensor", 1)
            if field != self.field:
                raise ValueError("literal scheme field disagrees with graph field")
            if len(data) != 3:
                raise ValueError("scheme leaf needs three coefficient families")
            target, shape = args[0], children[0].shape
            self._cap(len(data[0]), "max_terms")
            self._cap(len(data[0])*sum(shape), "max_scheme_scalars")
            TensorScheme(self.materialize_tensor(target), *data).require_exact()
            kind, terms = "scheme", len(data[0])
        elif op in ("scheme_product", "scheme_sum", "scheme_repeat", "scheme_permute"):
            field = expect("scheme", 2 if op == "scheme_product" else
                           1 if op in ("scheme_repeat", "scheme_permute") else None)
            kind = "scheme"
            tensor_op = op.replace("scheme", "tensor")
            target = self._add(tensor_op, tuple(c.target for c in children), data)
            shape = self.node(target).shape
            if op == "scheme_product":
                terms = children[0].terms*children[1].terms
            elif op == "scheme_sum":
                terms = sum(c.terms for c in children)
            elif op == "scheme_repeat":
                terms = _natural(data[0])*children[0].terms
            else:
                terms = children[0].terms
        elif op == "scheme_restrict":
            if len(children) != 4 or children[0].kind != "scheme":
                raise ValueError("scheme restriction requires scheme and three maps")
            if data:
                raise ValueError("scheme restriction takes no data")
            field, kind, terms = children[0].field, "scheme", children[0].terms
            target = self._add("tensor_restrict", (children[0].target,)+args[1:])
            shape = self.node(target).shape
        elif op == "scheme_descend":
            extension = expect("scheme", 1)
            if data:
                raise ValueError("descent takes no data")
            target = self._add("tensor_prime", (children[0].target,))
            field, kind, shape = self.node(target).field, "scheme", children[0].shape
            terms = extension.degree**2*children[0].terms
        else:
            raise ValueError(f"unknown construction rule: {op}")
        self._cap(len(self._nodes)+1, "max_nodes")
        node = ConstructionNode(op, args, data, kind, field, shape, terms, target, base)
        index = len(self._nodes)
        self._nodes.append(node)
        self._interned[key] = index
        return index

    def tensor(self, tensor):
        if not isinstance(tensor, FiniteTensor) or tensor.field != self.field:
            raise ValueError("literal tensor field disagrees with graph field")
        return self._add("tensor_literal", data=(tensor.shape, tensor.coefficients))

    def matrix(self, n, m=None, k=None):
        return self._add("tensor_matrix", data=(n, n if m is None else m, n if k is None else k))

    def mapping(self, mapping):
        if not isinstance(mapping, LinearMap) or mapping.field != self.field:
            raise ValueError("literal map field disagrees with graph field")
        return self._add("map_literal", data=(mapping.input_size, mapping.rows))

    def identity(self, size):
        return self._add("map_identity", data=(size,))

    def selection(self, size, indices):
        return self._add("map_selection", data=(size, tuple(indices)))

    def map_product(self, first, second):
        return self._add("map_product", (first, second))

    def compose(self, after, before):
        return self._add("map_compose", (after, before))

    def scheme(self, scheme):
        if not isinstance(scheme, TensorScheme) or scheme.field != self.field:
            raise ValueError("literal scheme field disagrees with graph field")
        return self._add("scheme_literal", (self.tensor(scheme.target),), (scheme.a, scheme.b, scheme.c))

    def tensor_product(self, first, second):
        return self._add("tensor_product", (first, second))

    def scheme_product(self, first, second):
        return self._add("scheme_product", (first, second))

    def direct_sum(self, indices):
        indices = tuple(indices)
        kind = self.node(indices[0]).kind if indices else None
        if kind not in ("tensor", "scheme"):
            raise ValueError("direct sum requires tensors or schemes")
        return self._add(kind+"_sum", indices)

    def repeat(self, index, count):
        kind = self.node(index).kind
        if kind not in ("tensor", "scheme"):
            raise ValueError("repeat requires tensor or scheme")
        return self._add(kind+"_repeat", (index,), (count,))

    def power(self, index, exponent):
        _natural(exponent)
        node = self.node(index)
        if node.kind not in ("tensor", "scheme"):
            raise ValueError("power requires tensor or scheme")
        if node.field != self.field:
            raise ValueError("take powers before final descent")
        unit = FiniteTensor.unit(self.field)
        result = self.tensor(unit) if node.kind == "tensor" else self.scheme(TensorScheme.from_tensor(unit))
        factor = index
        operation = self.tensor_product if node.kind == "tensor" else self.scheme_product
        while exponent:
            if exponent & 1:
                result = operation(result, factor)
            exponent //= 2
            if exponent:
                factor = operation(factor, factor)
        return result

    def permute(self, index, order):
        kind = self.node(index).kind
        if kind not in ("tensor", "scheme"):
            raise ValueError("permutation requires tensor or scheme")
        return self._add(kind+"_permute", (index,), tuple(order))

    def restrict(self, index, maps):
        kind = self.node(index).kind
        if kind not in ("tensor", "scheme"):
            raise ValueError("restriction requires tensor or scheme")
        return self._add(kind+"_restrict", (index,)+tuple(maps))

    def witness(self, source, target, maps, label):
        return self._add("tensor_witness", (source, target)+tuple(maps), (label,))

    def descend(self, scheme):
        return self._add("scheme_descend", (scheme,))

    def _tensor_value(self, index, coordinates, work):
        work.use()
        node = self.node(index)
        f, op = node.field, node.op
        if op == "tensor_literal":
            i, j, k = coordinates
            return node.data[1][(i*node.shape[1]+j)*node.shape[2]+k]
        if op == "tensor_matrix":
            n, m, k = node.data
            i, j = divmod(coordinates[0], m)
            j2, h = divmod(coordinates[1], k)
            i2, h2 = divmod(coordinates[2], k)
            return int(i == i2 and j == j2 and h == h2)
        if op == "tensor_product":
            a, b = node.args
            pairs = tuple(divmod(x, d) for x, d in zip(coordinates, self.node(b).shape))
            return f.mul(self._tensor_value(a, tuple(x[0] for x in pairs), work),
                         self._tensor_value(b, tuple(x[1] for x in pairs), work))
        if op in ("tensor_sum", "tensor_repeat"):
            if op == "tensor_repeat":
                source = node.args[0]
                pairs = tuple(divmod(x, d) for x, d in zip(coordinates, self.node(source).shape))
                return self._tensor_value(source, tuple(x[1] for x in pairs), work) if len({x[0] for x in pairs}) == 1 else 0
            offsets = [0, 0, 0]
            for child in node.args:
                work.use()
                shape = self.node(child).shape
                local = tuple(x-o for x, o in zip(coordinates, offsets))
                if all(0 <= x < d for x, d in zip(local, shape)):
                    return self._tensor_value(child, local, work)
                offsets = [o+d for o, d in zip(offsets, shape)]
            return 0
        if op == "tensor_permute":
            old = [0, 0, 0]
            for axis, value in zip(node.data, coordinates):
                old[axis] = value
            return self._tensor_value(node.args[0], tuple(old), work)
        if op == "tensor_witness":
            return self._tensor_value(node.args[1], coordinates, work)
        if op == "tensor_prime":
            child = node.args[0]
            return self.node(child).field.coordinates(self._tensor_value(child, coordinates, work))[0]
        if op == "tensor_restrict":
            source, *maps = node.args
            shape = self.node(source).shape
            work.use(prod(shape))
            total = 0
            for old in product(*(range(d) for d in shape)):
                scale = 1
                for mapping, output, input_ in zip(maps, coordinates, old):
                    scale = f.mul(scale, self._map_value(mapping, output, input_, work))
                if scale:
                    total = f.add(total, f.mul(scale, self._tensor_value(source, old, work)))
            return total
        raise ValueError("node is not a tensor")

    def _map_value(self, index, output, input_, work):
        work.use()
        node = self.node(index)
        f = node.field
        if node.op == "map_literal":
            return node.data[1][output][input_]
        if node.op == "map_identity":
            return int(output == input_)
        if node.op == "map_selection":
            return int(node.data[1][output] == input_)
        a, b = node.args
        if node.op == "map_product":
            o1, o2 = divmod(output, self.node(b).shape[1])
            i1, i2 = divmod(input_, self.node(b).shape[0])
            return f.mul(self._map_value(a, o1, i1, work), self._map_value(b, o2, i2, work))
        if node.op == "map_compose":
            width = self.node(a).shape[0]
            work.use(width)
            return f.sum(f.mul(self._map_value(a, output, h, work), self._map_value(b, h, input_, work)) for h in range(width))
        raise ValueError("node is not a map")

    def _scheme_value(self, index, axis, term, coordinate, work):
        work.use()
        node = self.node(index)
        f, op = node.field, node.op
        if op == "scheme_literal":
            return node.data[axis][term][coordinate]
        if op == "scheme_product":
            a, b = node.args
            q, r = divmod(term, self.node(b).terms)
            i, j = divmod(coordinate, self.node(b).shape[axis])
            return f.mul(self._scheme_value(a, axis, q, i, work), self._scheme_value(b, axis, r, j, work))
        if op == "scheme_repeat":
            child = node.args[0]
            copy, term = divmod(term, self.node(child).terms)
            slot, local = divmod(coordinate, self.node(child).shape[axis])
            return self._scheme_value(child, axis, term, local, work) if copy == slot else 0
        if op == "scheme_sum":
            offset = 0
            for child in node.args:
                work.use()
                c = self.node(child)
                if term < c.terms:
                    return self._scheme_value(child, axis, term, coordinate-offset, work) if offset <= coordinate < offset+c.shape[axis] else 0
                term -= c.terms
                offset += c.shape[axis]
            raise ValueError("term outside sum")
        if op == "scheme_permute":
            return self._scheme_value(node.args[0], node.data[axis], term, coordinate, work)
        if op == "scheme_restrict":
            child, *maps = node.args
            width = self.node(child).shape[axis]
            work.use(width)
            return f.sum(f.mul(self._map_value(maps[axis], coordinate, h, work), self._scheme_value(child, axis, term, h, work)) for h in range(width))
        if op == "scheme_descend":
            child = node.args[0]
            extension = self.node(child).field
            e = extension.degree
            original, pair = divmod(term, e*e)
            i, j = divmod(pair, e)
            value = self._scheme_value(child, axis, original, coordinate, work)
            if axis < 2:
                return extension.coordinates(value)[i if axis == 0 else j]
            basis = extension.mul(extension.p**i, extension.p**j)
            return extension.coordinates(extension.mul(basis, value))[0]
        raise ValueError("node is not a scheme")

    def tensor_coefficient(self, index, *coordinates):
        node = self.node(index)
        if node.kind != "tensor" or len(coordinates) != 3 or any(type(x) is not int or not 0 <= x < d for x, d in zip(coordinates, node.shape)):
            raise ValueError("tensor coordinate out of range")
        return self._tensor_value(index, coordinates, _Work(self.limits.max_query_work))

    def scheme_coefficient(self, index, axis, term, coordinate):
        node = self.node(index)
        if node.kind != "scheme" or type(axis) is not int or not 0 <= axis < 3 or type(term) is not int or not 0 <= term < node.terms or type(coordinate) is not int or not 0 <= coordinate < node.shape[axis]:
            raise ValueError("scheme coordinate out of range")
        return self._scheme_value(index, axis, term, coordinate, _Work(self.limits.max_query_work))

    def materialize_tensor(self, index):
        node = self.node(index)
        if node.kind != "tensor":
            raise ValueError("expected tensor node")
        self._cap(prod(node.shape), "max_tensor_entries")
        return FiniteTensor.from_function(node.field, node.shape, lambda *x: self.tensor_coefficient(index, *x))

    def materialize_map(self, index):
        node = self.node(index)
        if node.kind != "map":
            raise ValueError("expected map node")
        width, height = node.shape
        self._cap(width*height, "max_map_entries")
        return LinearMap(node.field, width, tuple(tuple(self._map_value(index, o, i, _Work(self.limits.max_query_work)) for i in range(width)) for o in range(height)))

    def materialize_scheme(self, index):
        node = self.node(index)
        if node.kind != "scheme":
            raise ValueError("expected scheme node")
        self._cap(node.terms, "max_terms")
        self._cap(node.terms*sum(node.shape), "max_scheme_scalars")
        target = self.materialize_tensor(node.target)
        families = tuple(tuple(tuple(self.scheme_coefficient(index, axis, q, i) for i in range(size)) for q in range(node.terms)) for axis, size in enumerate(node.shape))
        result = TensorScheme(target, *families)
        result.require_exact()
        return result

    def scheme_support_summary(self, index):
        """Exact support/addition counts where construction rules determine them.

        Linear substitution and extension projection can cause cancellation;
        those cases return unknown, rather than guessing from rank bounds.
        Counts describe the stored bilinear evaluation, not CPU/GPU runtime.
        """
        if self.node(index).kind != "scheme":
            raise ValueError("expected scheme node")
        memo = {}

        def support(i):
            if i in memo:
                return memo[i]
            node = self.node(i)
            if node.op == "scheme_literal":
                nnz = tuple(sum(bool(v) for row in family for v in row) for family in node.data)
                rows = tuple(sum(any(row) for row in family) for family in node.data)
                columns = tuple(sum(any(row[j] for row in family) for j in range(size))
                                for family, size in zip(node.data, node.shape))
                result = nnz, rows, columns
            elif node.op in ("scheme_product", "scheme_sum", "scheme_repeat", "scheme_permute"):
                children = tuple(support(c) for c in node.args)
                if any(c is None for c in children):
                    result = None
                elif node.op == "scheme_product":
                    result = tuple(tuple(x*y for x, y in zip(a, b))
                                   for a, b in zip(children[0], children[1]))
                elif node.op == "scheme_sum":
                    result = tuple(tuple(sum(c[group][axis] for c in children) for axis in range(3))
                                   for group in range(3))
                elif node.op == "scheme_repeat":
                    result = tuple(tuple(node.data[0]*x for x in group) for group in children[0])
                else:
                    result = tuple(tuple(group[axis] for axis in node.data) for group in children[0])
            else:
                result = None
            memo[i] = result
            return result

        counts = support(index)
        if counts is None:
            return {"status": "unknown_after_linear_substitution", "terms": self.node(index).terms}
        nnz, rows, columns = counts
        additions = nnz[0]-rows[0]+nnz[1]-rows[1]+nnz[2]-columns[2]
        return {"status": "exact", "terms": self.node(index).terms,
                "coefficient_nonzeros": nnz, "active_rows": rows,
                "active_coordinates": columns, "linear_additions": additions}

    def to_dict(self, root):
        self.node(root)
        return {"format": "matrix-coefficient-dag-v1", "field": {"p": self.field.p, "modulus": self.field.modulus},
                "root": root, "nodes": [{"op": n.op, "args": n.args, "data": n.data} for n in self._nodes]}

    @classmethod
    def from_dict(cls, payload, limits=GraphLimits()):
        if not isinstance(payload, dict) or set(payload) != {"format", "field", "root", "nodes"} or payload["format"] != "matrix-coefficient-dag-v1":
            raise ValueError("unsupported certificate format")
        field = payload["field"]
        if not isinstance(field, dict) or set(field) != {"p", "modulus"}:
            raise ValueError("invalid field record")
        graph = cls(FiniteField(field["p"], tuple(field["modulus"])), limits)
        if not isinstance(payload["nodes"], (tuple, list)):
            raise ValueError("nodes must be a sequence")
        graph._cap(len(payload["nodes"]), "max_nodes")
        translated = []
        for record in payload["nodes"]:
            if not isinstance(record, dict) or set(record) != {"op", "args", "data"} or not isinstance(record["args"], (tuple, list)):
                raise ValueError("invalid construction record")
            if any(type(i) is not int or not 0 <= i < len(translated) for i in record["args"]):
                raise ValueError("forward/cyclic node reference")
            translated.append(graph._add(record["op"], tuple(translated[i] for i in record["args"]), record["data"]))
        root = payload["root"]
        if type(root) is not int or not 0 <= root < len(translated):
            raise ValueError("invalid root node")
        return graph, translated[root]
