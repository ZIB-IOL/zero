# Certified periodic almost five-coloring

Removed area: **3.485287603537%**.
Exact removed fraction: **136749079/3923609600**.
The five real colors avoid distance exactly 1.

## Verify

Requires Python, NumPy, and a C++17 compiler for the fast integer checker.
No solver, neural model, GPU, or network is required.

```
python3 audit_sparse.py artifacts/sparse.json --native --output my_audit.json
```

Omitting `--native` uses the independent exhaustive Python reference. It is
substantially slower. The geometric predicates use exact rational/integer
arithmetic in both modes, with no floating-point distance tolerance.

## Complete construction

`artifacts/parent_chain.json` defines the root coloring, rational periodic
lattice, and dense parent coloring. Colors 0..4 are real; label 5 is removed.
Every real parent cell is retained unchanged. Subdivide each parent cell by
the refinement factor specified by the two manifests. `artifacts/paint.npz`
(array `paint`) assigns real colors to the listed fine children of removed
parent cells. Every other such child remains removed. Rows are integer
`(row, column, color)` triples, with row/column measured in the final grid
of shape (124480, 126080).

Use the same half-open parallelogram convention throughout. There are no
overlapping assignments or unspecified cells. The checker rejects duplicate
triples, out-of-range coordinates, and any overwrite of a real parent cell.
Removed area is the parent's removed count times the squared subdivision
factor, less the number of distinct additions, divided by the total number
of fine cells. There is no sampled area estimate.

The dense parent is checked directly against the root. Every newly colored fine cell is then
checked against all possible forbidden-distance neighbors, including every
relevant periodic translation, inherited real cells, and other additions.
Unchanged cells inherit validity from the certified parent. Closed cells are
checked conservatively, which suffices for the half-open construction.

The fine grid is represented sparsely and need not be allocated to verify.
This proves feasibility and removed area; it does not establish a globally
optimal coloring.

Paint SHA256: `06bbc295b3f062303b8d2406a651b3dd131c4af6db948c1991ed2c9152c1f800`.
