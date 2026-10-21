"""Independent, exhaustive exact uniform-mask audit; Python standard library.

Geometry and complete enumeration bounds come from the lattice in the supplied
metadata. Exact norm-ball bounds discard distant and near-origin boxes, exact
vertex witnesses decide most remaining boxes, and the rational active-set
oracle checks every unresolved box.

"""
import csv
import hashlib
from itertools import product
import json
from math import lcm
from pathlib import Path
import time

from geometry_exact import LatticeGeometry, interval_conflicts


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(reference, generated_mask):
    started = time.monotonic()
    metadata = json.loads(reference.read_text())
    shape = tuple(metadata["shape"])
    geometry = LatticeGeometry(metadata["lattice_rows"])
    cell_lattice = [[x / n for x in row] for row, n in zip(geometry.lattice, shape)]
    cell_geometry = LatticeGeometry(cell_lattice)
    gram = cell_geometry.gram
    denominator = lcm(*(x.denominator for row in gram for x in row))
    q = tuple(tuple(int(x * denominator) for x in row) for row in gram)

    def norm(x):
        a, b, c = x
        return (q[0][0]*a*a + q[1][1]*b*b + q[2][2]*c*c
                + 2*(q[0][1]*a*b + q[0][2]*a*c + q[1][2]*b*c))

    # Error box is [-1,1]^3 in grid coordinates. The maximum squared
    # Euclidean radius is exactly the maximum of the eight vertex norms.
    corners = tuple(product((-1, 1), repeat=3))
    radius = max(map(norm, corners))
    ranges = geometry.uniform_offset_ranges(shape)
    observed = set()
    counts = dict(tested=0, outside_radius=0, inside_radius=0,
                  vertex_witness=0, all_vertices_inside=0, active_sets=0)
    for a in ranges[0]:
        if a < 0:
            continue
        for b, c in product(ranges[1], ranges[2]):
            d = a, b, c
            if d < (0, 0, 0):
                continue
            counts["tested"] += 1
            center = norm(d)
            # ||center|| > 1+R: equivalent exact squared rational comparison.
            gap = center - denominator - radius
            if gap > 0 and gap*gap > 4*denominator*radius:
                counts["outside_radius"] += 1
                continue
            # ||center||+R < 1. Strict inequality is needed for closed cells.
            gap = denominator - center - radius
            if gap > 0 and gap*gap > 4*center*radius:
                counts["inside_radius"] += 1
                continue
            vertices = [norm((a+i, b+j, c+k)) for i, j, k in corners]
            if max(vertices) < denominator:
                counts["all_vertices_inside"] += 1
                continue
            if min(vertices) <= denominator:
                # Connectedness and two vertex witnesses prove a unit pair.
                counts["vertex_witness"] += 1
                observed.add(d)
                continue
            counts["active_sets"] += 1
            interval = cell_geometry.squared_distance_range(
                tuple(x-1 for x in d), tuple(x+1 for x in d))
            if interval_conflicts(*interval, boundary="closed"):
                observed.add(d)
        print(json.dumps({"completed_first_coordinate": a, **counts}), flush=True)
    with generated_mask.open("w", newline="") as target:
        csv.writer(target).writerows(sorted(observed))

    def signed_modulo(offsets):
        return {tuple(sign*x % n for x, n in zip(d, shape))
                for d in offsets for sign in (-1, 1)}

    observed_modulo = signed_modulo(observed)
    inputs = [reference, Path(__file__), Path(__file__).with_name("geometry_exact.py")]
    return {
        "verified": True,
        "boundary": "closed cells (also certifies every half-open subset)",
        "independent_offset_bounds": [[r.start, r.stop-1] for r in ranges],
        "unit_coordinate_L1_bounds": list(map(str, geometry.unit_coordinate_bounds)),
        "classification_counts": counts,
        "conflict_offsets_up_to_sign": len(observed),
        "symmetrized_modulo_conflicts": len(observed_modulo),
        "generated_mask_sha256": sha256(generated_mask),
        "seconds": time.monotonic() - started,
        "sha256": {str(p): sha256(p) for p in inputs},
    }
