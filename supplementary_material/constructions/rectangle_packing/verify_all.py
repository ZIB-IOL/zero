"""Independently validate every bundled packing construction using integer arithmetic."""

from __future__ import annotations

import itertools
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
CONSTRUCTIONS = ROOT / "constructions"


def separate_polygons(a, b):
    for polygon in (a, b):
        for p, q in zip(polygon, polygon[1:] + polygon[:1]):
            axis = (q[1] - p[1], p[0] - q[0])
            pa = [sum(x * y for x, y in zip(v, axis)) for v in a]
            pb = [sum(x * y for x, y in zip(v, axis)) for v in b]
            if max(pa) <= min(pb) or max(pb) <= min(pa):
                return True
    return False


def validate(data, family):
    n, side = data["n"], data["side"]
    assert type(n) is int and n > 0 and type(side) is int and side > 0
    lower_bound = None

    if family in ("2d", "3d"):
        key = {"2d": "rectangles", "3d": "boxes"}[family]
        pieces = data[key]
        dims = 3 if family == "3d" else 2
        labels = "abc"[:dims]
        positions = "xyz"[:dims]
        lengths = "whd"[:dims]
        assert len(pieces) == n
        assert sorted(r[k] for r in pieces for k in labels) == list(range(1, dims * n + 1))
        for piece in pieces:
            assert all(type(piece[k]) is int for k in labels + positions + lengths)
            assert sorted(piece[k] for k in lengths) == sorted(piece[k] for k in labels)
            assert all(piece[k] > 0 for k in lengths)
            assert all(
                0 <= piece[p] and piece[p] + piece[length] <= side
                for p, length in zip(positions, lengths)
            )
        assert all(
            any(
                a[p] + a[length] <= b[p] or b[p] + b[length] <= a[p]
                for p, length in zip(positions, lengths)
            )
            for a, b in itertools.combinations(pieces, 2)
        )
        if family == "2d":
            area = sum(i * (2 * n + 1 - i) for i in range(1, n + 1))
            lower_bound = math.isqrt(area - 1) + 1

    elif family == "equilateral":
        pieces = data["triangles"]
        assert sorted(r["side"] for r in pieces) == list(range(1, n + 1))
        polygons = []
        for piece in pieces:
            a, b = piece["anchor"]
            k = piece["side"]
            assert piece["orientation"] in ("U", "D")
            expected = (
                [[a, b], [a + k, b], [a, b + k]]
                if piece["orientation"] == "U"
                else [[a + k, b], [a, b + k], [a + k, b + k]]
            )
            assert piece["vertices"] == expected
            assert all(
                type(x) is int and type(y) is int and x >= 0 and y >= 0 and x + y <= side
                for x, y in expected
            )
            polygons.append(expected)
        assert all(separate_polygons(a, b) for a, b in itertools.combinations(polygons, 2))
        lower_bound = math.isqrt(sum(k * k for k in range(1, n + 1)) - 1) + 1

    elif family == "tetrahedra":
        pieces = data["tetrahedra"]
        assert sorted(r["side"] for r in pieces) == list(range(1, n + 1))
        bounds = []
        for piece in pieces:
            k, anchor = piece["side"], piece["anchor"]
            assert piece["orientation"] in ("U", "D") and len(anchor) == 4
            up = piece["orientation"] == "U"
            expected = [
                [anchor[j] + (k if (j == i if up else j != i) else 0) for j in range(4)]
                for i in range(4)
            ]
            assert piece["vertices"] == expected
            assert all(
                all(type(v) is int and v >= 0 for v in point) and sum(point) == side
                for point in expected
            )
            bounds.append((anchor if up else [0] * 4, [side] * 4 if up else [v + k for v in anchor]))
        for (lo_a, hi_a), (lo_b, hi_b) in itertools.combinations(bounds, 2):
            lo = [max(a, b) for a, b in zip(lo_a, lo_b)]
            hi = [min(a, b) for a, b in zip(hi_a, hi_b)]
            overlap = all(a < b for a, b in zip(lo, hi)) and sum(lo) < side < sum(hi)
            assert not overlap
        volume = sum(k**3 for k in range(1, n + 1))
        lower_bound = 1
        while lower_bound**3 < volume:
            lower_bound += 1

    else:
        raise ValueError(f"Unknown family: {family}")

    return {
        "n": n,
        "side": side,
        "valid": True,
        "area_or_volume_lower_bound": lower_bound,
        "optimal_by_this_bound": lower_bound == side if lower_bound is not None else None,
    }


def main():
    results = []
    for path in sorted(CONSTRUCTIONS.glob("*/*.json")):
        family = path.parent.name
        result = validate(json.loads(path.read_text()), family)
        result["path"] = str(path.relative_to(ROOT))
        results.append(result)
        print(
            f"PASS {result['path']}: side={result['side']}, "
            f"bound={result['area_or_volume_lower_bound']}, "
            f"optimal_by_bound={result['optimal_by_this_bound']}"
        )
    destination = ROOT / "evidence" / "validation_results.json"
    destination.write_text(json.dumps({"results": results}, indent=2) + "\n")
    print(f"Validated {len(results)} bundled constructions; wrote {destination.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
