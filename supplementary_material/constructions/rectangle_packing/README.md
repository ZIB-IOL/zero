# Rectangle-packing constructions

This directory contains the machine-readable constructions supporting the paper's rectangle-packing study.

- `constructions/2d`: the `n=15` construction and all seven claimed `n=18,...,24` constructions.
- `constructions/3d`: the saved axis-aligned cuboid constructions.
- `constructions/equilateral`: the saved equilateral-triangle constructions.
- `constructions/tetrahedra`: the saved tetrahedron constructions.
- `verify_all.py`: an independent standard-library checker for dimensions, labels, containment, pairwise interior disjointness, and elementary area/volume lower bounds.
- `evidence/validation_results.json`: the complete bundled validation report.
- `checkers/`: the standalone single-layout 2D checker.

Run all bundled checks from this directory:

```bash
python3 verify_all.py
```

The layouts are frozen copies of the constructions evaluated in the paper. The manifest at the bundle root records file hashes without exposing private source locations or repository history.

Passing a feasibility check does not by itself prove optimality. `optimal_by_bound=true` means the displayed side matches the elementary area or volume lower bound implemented by the checker.
