"""Check a complete integer rectangle layout and optionally time verification."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import statistics
import timeit


def verify(layout: dict) -> dict:
    n, side = layout["n"], layout["side"]
    if type(n) is not int or type(side) is not int or min(n, side) <= 0:
        raise ValueError("n and side must be positive integers")
    rectangles = layout["rectangles"]
    if len(rectangles) != n or layout.get("missing"):
        raise ValueError("Layout is incomplete")
    labels, area = [], 0
    for i, r in enumerate(rectangles):
        if any(type(r[k]) is not int for k in ("x", "y", "w", "h", "a", "b")):
            raise ValueError("Coordinates and dimensions must be integers")
        x, y, w, h, a, b = (r[k] for k in ("x", "y", "w", "h", "a", "b"))
        if min(w, h) <= 0 or sorted((w, h)) != sorted((a, b)):
            raise ValueError(f"Invalid dimensions for rectangle {i}")
        if x < 0 or y < 0 or x + w > side or y + h > side:
            raise ValueError(f"Rectangle {i} is outside the square")
        labels.extend((a, b))
        area += w * h
        for j in range(i):
            s = rectangles[j]
            if (x < s["x"] + s["w"] and s["x"] < x + w
                    and y < s["y"] + s["h"] and s["y"] < y + h):
                raise ValueError(f"Rectangles {j} and {i} overlap")
    if sorted(labels) != list(range(1, 2 * n + 1)):
        raise ValueError("Labels 1..2n must appear exactly once")
    minimum_area = sum(i * (2 * n + 1 - i) for i in range(1, n + 1))
    return dict(valid=True, n=n, side=side, area=area, empty_area=side * side - area,
                pairwise_checks=n * (n - 1) // 2,
                minimum_pairing_area=minimum_area,
                integer_square_side_optimal=(side - 1) ** 2 < minimum_area)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("layout", type=Path)
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    if args.iterations < 1:
        parser.error("--iterations must be positive")
    raw = args.layout.read_bytes()
    layout = json.loads(raw)
    result = verify(layout)
    compute = timeit.repeat(lambda: verify(layout), number=args.iterations, repeat=5)
    with_io = timeit.repeat(lambda: verify(json.loads(args.layout.read_bytes())),
                            number=args.iterations, repeat=5)
    result.update(source=str(args.layout.resolve()), sha256=hashlib.sha256(raw).hexdigest(),
                  benchmark_iterations_per_repeat=args.iterations, benchmark_repeats=5,
                  median_verification_seconds=statistics.median(compute) / args.iterations,
                  median_read_parse_verify_seconds=statistics.median(with_io) / args.iterations,
                  timing_excludes_python_startup=True,
                  read_timing_note="Repeated reads; filesystem cache is warm")
    text = json.dumps(result, indent=2) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
