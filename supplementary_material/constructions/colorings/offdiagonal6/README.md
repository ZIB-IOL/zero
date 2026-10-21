# Fixed off-diagonal six-colorings

`constructions/` contains the five fixed periodic colorings used in the paper. The first three have rational coordinates, the fourth reconstructs the published algebraic baseline, and the fifth certifies the new algebraic upper endpoint.

Run the rational checkers from this directory:

```bash
python3 checkers/check.py constructions/lower.json
python3 checkers/check.py constructions/lower-wide.json
python3 checkers/check.py constructions/bridge.json
```

The remaining checks require SymPy:

```bash
python3 checkers/check_reference.py --output /tmp/published-fixed-check.json
python3 checkers/check_algebraic.py constructions/upper.json --output /tmp/upper-check.json
```
