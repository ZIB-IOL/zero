# Certified periodic almost-14-coloring

**3.4622037578% → 1.9851143373% removed volume**, a **42.6633% reduction**.

The complete construction is given by `artifacts/candidate2.json` and `artifacts/candidate2.labels.bin`.

Verify with Python's standard library and a C++17 compiler:

```bash
python3 verify.py --out /tmp/almost14-verification
```

The verifier derives the exact geometric conflict mask, checks every label, and confirms the reported removed fraction. The construction is feasible; no global optimality claim is made.
