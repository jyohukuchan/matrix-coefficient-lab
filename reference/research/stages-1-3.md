# Cost screening, structured certificates, and finite proof-row integration

This checkpoint implements the first three stages of the revised coefficient
extraction plan. It supplies a bounded, reproducible foundation for discovery;
**no useful d=2,k=5 certificate or new exponent-2.4 scheme was found**.
The new Python checker has not been formalized in Lean. Existing Lean sources
and their protected specification are unchanged.

## Reproduce

From the repository root, with Python 3.10 or newer and no third-party packages:

```sh
python3 -m reference.proof_pipeline --output /tmp/matrix-coefficient-control
python3 -m reference.proof_pipeline --k 5 --tau 12/5 --output /tmp/matrix-coefficient-k5
python3 -m reference.proof_pipeline --proof-only
python3 -m reference.structured_demo
python3 -m unittest reference.test_resource_planning reference.test_structured_certificates reference.test_proof_pipeline -v
```

The default pipeline produces `report.json` and `construction-dag.json`.
The latter includes all generated ordinary witnesses and a final exact
49-term known 4-by-4 matrix control. Its root is a supplied, checked scheme.
Replay it independently of the search:

```python
import json
from pathlib import Path
from reference.structured_certificates import CertificateGraph

graph, root = CertificateGraph.from_dict(json.loads(
    Path('/tmp/matrix-coefficient-control/construction-dag.json').read_text()))
node = graph.node(root)
assert node.kind == 'scheme' and node.terms == 49
assert graph.materialize_scheme(root).verify()
```

Saving and loading checks every node, including nodes outside the root's
ancestor set. Field representation, node types, coordinate order, decomposition
lengths and shapes are reconstructed, rather than read as asserted metadata.
Forward/cyclic references, unknown rules, invalid leaf decompositions and false
restriction identities are rejected. A cap produces an incomplete result.

Checkpoint verification on 2026-10-08: all **219 reference tests passed** in
158.442 seconds, including 23 new focused tests. The protected-source check
also passed for all eleven files and the frozen Challenge model. No Lean
sources were changed or rebuilt for this Python checkpoint.

## Stage 1: numerical screening before coefficient allocation

`reference.resource_planning` evaluates each capped block size using

```text
g_(b+1) = k*g_b + m*d^b
C_(b+1) = k*C_b + R_D*(d^b)^3
Q_b     = k^b*R_S - g_b
r_next  = C_b + Q_b*r
```

It records the first strict-gap iteration within the cap, or the last attempted
iteration when no gap is reached. Integer comparisons use rational exponents;
no floating-point approximation decides a rank inequality. One final descent
charges `extension_degree**2` once.

Each block records final matrix size, extension/prime-field decomposition
lengths, final dense tensor entries (`n**6`), coefficient-family scalars
(`3*r*n**2`), ideal bit-packed family bytes, and powered catalyst-map entries.
Packed bytes are not Python peak memory. Budget checks cover these counts;
intermediate absorption maps may require still more space. These are allocation
screening measures, not speed forecasts or a complete resource feasibility proof.
All capped blocks are retained for comparison rather than selecting the first
sufficient block alone.

`forecast_candidate(CandidateCosts(...), Fraction(...))` accepts **unverified
numerical assumptions**. `forecast_verified_inputs(certificate, S_scheme,
D_scheme, Fraction(...))` first checks actual maps and exact decompositions.
Both return `output_status='numerical_only'`: neither supplies final coefficients.

For the hypothetical values d=2,k=5,m=1,R_S=4,R_D=0, the first sufficient
block for tau=12/5 is b=25, with n=33,554,432 and
Q=1,092,751,820,893,346,269. This is only a size warning about the current
recurrence; it does not establish that such a certificate exists, or exclude
better seed decompositions or another conversion strategy.

## Stage 2: exact construction DAG

`reference.structured_certificates.CertificateGraph` supports:

- Analytic matrix tensors and checked finite tensor/decomposition leaves.
- Tensor products, direct sums, integer repetition, powers and axis permutation.
- Dense map leaves, symbolic identity/selection, Kronecker maps and composition.
- Local tensor/scheme substitutions and checked ordinary restriction identities.
- Unchanged base-field descent after powering in one fixed extension.
- Individual tensor/map-family coefficient queries and bounded dense expansion.
- Exact nonzero-support and linear-addition counts for operations where support
  factorizes; cancellation after general substitution/projection is marked unknown.

The correctness argument is inductive: every scheme leaf satisfies all target
coefficient identities, and each allowed construction acts on the target and
its coefficient families by the same exact formula. An ordinary witness is
accepted only after its full map equality has been checked on bounded finite
inputs. In particular, a detector row is never stored as an ordinary restriction.
Decomposition lengths are rank **upper bounds**, not minimum rank or measured
operation counts. Addition counts omit coefficient-scaling multiplications and
data movement.

The huge-power demo stores a supplied seven-term scheme's 25th tensor power in
18 nodes and about 1.6 KB of JSON. Its axes each have size 4**25 and it has 7**25
terms. It replays the construction and queries specified coefficients without
allocating the full target or families. This uses product-coordinate ordering;
identifying it with row-major coordinates of a larger matrix needs the usual
coordinate reindexing. It reconstructs a known scheme, with no new rank gain.

Current limits are deliberate:

- Interpolation, Fourier and determinant rows enter as **small, fully checked
  witness seeds** with their actual finite maps and overhead. Their formulas
  are not yet separate scalable symbolic operation types.
- General witness equalities still require bounded dense verification. The DAG
  can take huge powers of a checked seed; it cannot certify an arbitrary giant
  supplied map from dimensions or rank estimates alone.
- Scalar-gain elimination still uses the existing dense compiler. This graph
  does not automatically compress the whole giant catalytic compilation route.
- Query work, graph nodes and expansion counts have separate caps. A compact
  description can still be expensive to evaluate. Python verification is not
  a Lean/kernel proof of this checker.

## Stage 3: generated proof rows, exact search, attribution

`reference.proof_pipeline` automatically generates four witnessed finite row
families over F_11: sector interpolation, exact-type word pulls, Fourier/square
separation and determinant filtration. It preserves node/copy overhead and
interns literally equal tensors across the row registries. The exact-type target
and Fourier source thereby share a state coordinate and really connect.

The exact-type/Fourier test uses supplied scalar branches with coefficients 1
and 3. The sector and determinant rows are separate finite proof-family
controls. This is integration of finite constructors, not an implementation of
the spectral existence argument. `PipelineConfig` exposes bounded sector/type/
determinant parameters; F_11 is fixed for the control, so Fourier root/node
availability and existing allocation caps still constrain valid configurations.

The optional known control is a checked seven-term Strassen decomposition.
An M_2 detector on the unit is hypothetical as a state inequality, with its
registered product verified exactly. The pipeline then searches rational
independent supports, clears denominators over the integers, assembles actual
catalytic maps, forecasts cost and runs the bounded coefficient compiler.

Default results:

| Item | Outcome |
| --- | --- |
| Generated ordinary proof rows | Four, with exact local maps |
| Full finite search, d=2,k=8 | Found after 21 supports; gain 1 |
| Used rows | Known-rank upper row and unit detector only |
| Remove any one proof family | Known-control dual still found |
| Remove the known control | Selected finite family exhausted |
| Controls alone | Same known-control route found |
| Dense known reconstruction | 2-by-2: 15 terms simplify to 7; 4-by-4: 63 simplify to 49 |
| Gain-preserving compiler, tau=4 | Exact 2-by-2 decomposition, 14 supplied terms |
| Same selected family, d=2,k=5 | Finite family exhausted |

The report contains used row indices and integer weights, plus a search with
each family removed and a controls-only search. A zero weight only means the
chosen dual did not use that row. A family-removal result compares existence
within this capped finite problem; it is not a universal claim about the
family's usefulness. Search caps are reported separately from finite exhaustion.

## Next frontier

The first three stages now have an executable bounded implementation and
independent dense/negative controls. Discovery remains open. The next research
step is to supply richer **connected** proof-derived row families and compare
bounded exploration with direct low-rank search. Both need explicit budgets
and switching criteria. Scalable symbolic interpolation/Fourier rules, symbolic
scalar elimination, and Lean formalization of the new finite catalytic argument
are further verification/construction work; they should not be confused with
finding a useful witness.
