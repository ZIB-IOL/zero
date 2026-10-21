"""Standalone exact Euclidean whole-box geometry in a rational 3-D lattice.

Lattice vectors are ROWS: x in lattice coordinates means sum(x[i]*B[i]).
Decimal strings are parsed exactly. Floats are deliberately rejected: use the
decimal text that defines the construction. Cells are (lower, upper) boxes in
lattice coordinates, of any positive widths. No third-party packages required.

The closure's squared distance range is obtained from eight vertices and all
27 active sets (free coordinates / lower / upper). This includes face interiors.
Closed boxes conflict when min <= 1 <= max. Open boxes conflict when min < 1 <
max; removing every cell boundary has zero volume and gives this latter exact
predicate. The open convention is NOT a certificate for half-open boxes.
"""

from fractions import Fraction
from itertools import product
import argparse
import json


def rational(value):
    if isinstance(value, float):
        raise TypeError("Specify exact rational numbers or decimal strings, not floats")
    return Fraction(value)


def vector(values):
    result = tuple(map(rational, values))
    if len(result) != 3:
        raise ValueError("Expected three coordinates")
    return result


def inverse(matrix):
    """Exact Gauss-Jordan inverse; also supports principal 1x1 and 2x2 blocks."""
    n = len(matrix)
    a = [[rational(x) for x in row] + [Fraction(i == j) for j in range(n)]
         for i, row in enumerate(matrix)]
    for j in range(n):
        pivot = next((i for i in range(j, n) if a[i][j]), None)
        if pivot is None:
            raise ValueError("Singular lattice or matrix")
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [x / scale for x in a[j]]
        for i in range(n):
            if i != j:
                scale = a[i][j]
                a[i] = [x - scale * y for x, y in zip(a[i], a[j])]
    return tuple(tuple(row[n:]) for row in a)


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def floor_fraction(value):
    return value.numerator // value.denominator


def cell_box(cell):
    lower, upper = map(vector, cell)
    if any(l >= u for l, u in zip(lower, upper)):
        raise ValueError("A cell must have positive width on every axis")
    return lower, upper


def box_volume(cell):
    """Volume as a fraction of a unit lattice fundamental domain."""
    lower, upper = cell_box(cell)
    result = Fraction(1)
    for l, u in zip(lower, upper):
        result *= u - l
    return result


class LatticeGeometry:
    def __init__(self, lattice):
        self.lattice = tuple(map(vector, lattice))
        if len(self.lattice) != 3:
            raise ValueError("Expected three lattice vectors")
        # Coordinates = inverse(B transpose) * physical column vector.
        physical_matrix = tuple(zip(*self.lattice))
        self.inverse_coordinates = inverse(physical_matrix)
        self.unit_coordinate_bounds = tuple(
            sum(map(abs, row)) for row in self.inverse_coordinates)
        self.gram = tuple(tuple(sum(a * b for a, b in zip(x, y))
                                for y in self.lattice) for x in self.lattice)
        self.active_sets = []
        for status in product((-1, 0, 1), repeat=3):
            free = tuple(i for i in range(3) if status[i] == 0)
            fixed = tuple(i for i in range(3) if status[i] != 0)
            inv = inverse([[self.gram[i][j] for j in free] for i in free]) if free else ()
            self.active_sets.append((status, free, fixed, inv))

    def squared_norm(self, x):
        return sum(x[i] * self.gram[i][j] * x[j] for i in range(3) for j in range(3))

    def squared_distance_range(self, lower, upper):
        """Exact min/max on the CLOSED lattice-coordinate box [lower, upper]."""
        lower, upper = vector(lower), vector(upper)
        if any(l > u for l, u in zip(lower, upper)):
            raise ValueError("Reversed box bounds")
        maximum = max(self.squared_norm(x) for x in product(*zip(lower, upper)))
        minimum = None
        for status, free, fixed, inv in self.active_sets:
            x = [lower[i] if status[i] == -1 else upper[i] for i in range(3)]
            rhs = [-sum(self.gram[i][j] * x[j] for j in fixed) for i in free]
            for k, i in enumerate(free):
                x[i] = sum(inv[k][j] * rhs[j] for j in range(len(free)))
            if any(x[i] < lower[i] or x[i] > upper[i] for i in free):
                continue
            norm = self.squared_norm(x)
            if minimum is None or norm < minimum:
                minimum = norm
        assert minimum is not None  # All vertices are feasible active sets.
        return minimum, maximum

    def translation_ranges(self, cell_a, cell_b):
        """Complete integer bounds for copies b+t that can reach distance <= 1.

        If |physical|_2 <= 1, every physical coordinate has absolute value <= 1,
        so |lattice_coordinate_i| <= L1(row_i(inverse(B transpose))). With
        difference coordinate in [b_low-a_high+t, b_high-a_low+t], interval
        intersection gives the bounds below. No short-vector assumption occurs.
        """
        a_low, a_high = cell_box(cell_a)
        b_low, b_high = cell_box(cell_b)
        return tuple(range(ceil_fraction(-r - bh + al),
                           floor_fraction(r - bl + ah) + 1)
                     for r, al, ah, bl, bh in zip(self.unit_coordinate_bounds,
                                                 a_low, a_high, b_low, b_high))

    def cell_pair_distance_range(self, cell_a, cell_b, translation=(0, 0, 0)):
        a_low, a_high = cell_box(cell_a)
        b_low, b_high = cell_box(cell_b)
        t = vector(translation)
        return self.squared_distance_range(
            [bl - ah + ti for bl, ah, ti in zip(b_low, a_high, t)],
            [bh - al + ti for bh, al, ti in zip(b_high, a_low, t)])

    def periodic_conflicts(self, cell_a, cell_b, *, boundary="closed"):
        """Yield (integer translation, exact min squared, exact max squared).

        Same-cell and zero-translation interactions are intentionally included.
        A large individual cell can contain a forbidden pair on its own.
        """
        check_boundary(boundary)
        for t in product(*self.translation_ranges(cell_a, cell_b)):
            minimum, maximum = self.cell_pair_distance_range(cell_a, cell_b, t)
            if interval_conflicts(minimum, maximum, boundary=boundary):
                yield t, minimum, maximum

    def uniform_offset_distance_range(self, shape, offset):
        """Distance range for grid index difference offset, without modulo."""
        shape = validate_shape(shape)
        offset = vector(offset)
        return self.squared_distance_range(
            [(d - 1) / n for d, n in zip(offset, shape)],
            [(d + 1) / n for d, n in zip(offset, shape)])

    def uniform_offset_ranges(self, shape):
        """Complete NON-period-reduced grid-offset ranges for unit conflicts."""
        shape = validate_shape(shape)
        return tuple(range(ceil_fraction(-n * r - 1), floor_fraction(n * r + 1) + 1)
                     for n, r in zip(shape, self.unit_coordinate_bounds))


def validate_shape(shape):
    shape = tuple(shape)
    if len(shape) != 3 or any(type(n) is not int or n <= 0 for n in shape):
        raise ValueError("Grid shape must contain three positive integers")
    return shape


def check_boundary(boundary):
    if boundary not in ("closed", "open"):
        raise ValueError("Boundary convention must be 'closed' or 'open'")


def interval_conflicts(minimum, maximum, *, boundary="closed"):
    check_boundary(boundary)
    return minimum <= 1 <= maximum if boundary == "closed" else minimum < 1 < maximum


def squared_distance_range(lattice, lower, upper):
    return LatticeGeometry(lattice).squared_distance_range(lower, upper)


def periodic_conflicts(lattice, cell_a, cell_b, *, boundary="closed"):
    return LatticeGeometry(lattice).periodic_conflicts(cell_a, cell_b, boundary=boundary)


def uniform_offset_conflicts(lattice, shape, offset, *, boundary="closed"):
    interval = LatticeGeometry(lattice).uniform_offset_distance_range(shape, offset)
    return interval_conflicts(*interval, boundary=boundary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lattice", required=True, help='JSON 3x3 row lattice, decimals as strings')
    parser.add_argument("--shape", nargs=3, type=int, required=True)
    parser.add_argument("--offset", nargs=3, type=int, required=True)
    parser.add_argument("--boundary", choices=("closed", "open"), default="closed")
    args = parser.parse_args()
    lattice = json.loads(args.lattice, parse_float=str)
    geometry = LatticeGeometry(lattice)
    lo, hi = geometry.uniform_offset_distance_range(args.shape, args.offset)
    print(json.dumps({"minimum_squared": str(lo), "maximum_squared": str(hi),
                      "conflict": interval_conflicts(lo, hi, boundary=args.boundary),
                      "boundary": args.boundary}))


if __name__ == "__main__":
    main()
