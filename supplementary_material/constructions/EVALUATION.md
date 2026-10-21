# Verification

Run the checks from the indicated directories.

## Rectangle packing

```bash
cd rectangle_packing
python3 verify_all.py
```

## Almost-five-coloring

Requires Python, NumPy, and a C++17 compiler.

```bash
cd colorings/almost5
python3 audit_sparse.py artifacts/sparse.json --native --output /tmp/almost5-audit.json
```

## Almost-14-coloring

Requires Python's standard library and a C++17 compiler.

```bash
cd colorings/almost14
python3 verify.py --out /tmp/almost14-verification
```

## Off-diagonal six-colorings

The first three constructions use rational coordinates and can be checked with the standard library. The final two use SymPy for exact algebraic-number arithmetic. See `colorings/offdiagonal6/README.md`.
