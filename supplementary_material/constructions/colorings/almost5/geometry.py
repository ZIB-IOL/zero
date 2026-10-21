"""Exact unit-distance mask for periodic parallelogram cells.

The CSV decimal lattice coordinates DEFINE the geometry as rational numbers.
All predicates and enumeration bounds use Python integers, without tolerances.
Vectors are rows; cell (i,j) is i*a+j*b+[0,1)*a+[0,1)*b.
Conflicts are conservatively checked on cell closures.
"""

import csv
from fractions import Fraction
from math import gcd, isqrt, lcm
from pathlib import Path

import numpy as np


def read_lattice(path, refinement=1):
    rows = {r['vector_id']: r for r in csv.DictReader(Path(path).open())}
    vectors = [[Fraction(rows[k][axis]) for axis in ('x', 'y')]
               for k in ('v1', 'v2')]
    shape = tuple(int(rows[k]['subdivisions']) * refinement for k in ('v1', 'v2'))
    return vectors, shape


def integer_basis(vectors, shape):
    cells = [[Fraction(x) / n for x in v] for v, n in zip(vectors, shape)]
    denominator = lcm(*(x.denominator for v in cells for x in v))
    a, b = [tuple(int(x * denominator) for x in v) for v in cells]
    if a[0] * b[1] == a[1] * b[0]:
        raise ValueError('Singular lattice')
    return a, b, denominator


def ceil_div(a, b):
    return -((-a) // b)


def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def intersects_unit_circle(q, r, a, b, denominator):
    """Whether distance range of the two CLOSED cells contains exactly 1."""
    center = (q*a[0]+r*b[0], q*a[1]+r*b[1])
    vertices = [(center[0]+s*a[0]+t*b[0], center[1]+s*a[1]+t*b[1])
                for s, t in ((-1,-1), (1,-1), (1,1), (-1,1))]
    radius2 = denominator**2
    if max(dot(v,v) for v in vertices) < radius2:
        return False
    if abs(q) <= 1 and abs(r) <= 1:
        return True  # origin is in the difference parallelogram
    for v, w in zip(vertices, vertices[1:]+vertices[:1]):
        edge = (w[0]-v[0], w[1]-v[1])
        ee = dot(edge, edge)
        projection = -dot(v, edge)
        if projection <= 0:
            close = dot(v, v) <= radius2
        elif projection >= ee:
            close = dot(w, w) <= radius2
        else:
            cross = v[0]*edge[1] - v[1]*edge[0]
            close = cross**2 <= radius2 * ee
        if close:
            return True
    return False


def unit_distance_offsets(vectors, shape):
    """Complete exact mask, with all lattice copies and self-conflicts.

    A difference-cell displacement from its center has Euclidean norm at most
    R = ||a||_1+||b||_1. Only centers in [D-R,D+R] can intersect radius D.
    Enumerate this annulus using integer square roots of row quadratics.
    Rounded outer bounds expand the disk; inner bounds only remove points
    strictly inside the inner disk. Then apply exact segment predicates.
    """
    a, b, denominator = integer_basis(vectors, shape)
    reach = sum(map(abs, a+b))
    outer = denominator + reach
    inner = max(0, denominator - reach)
    det = abs(a[0]*b[1]-a[1]*b[0])
    qmax = ceil_div(outer * (abs(b[0])+abs(b[1])), det)
    aa, ab, bb = dot(a,a), dot(a,b), dot(b,b)
    offsets = set()
    for q in range(-qmax, qmax+1):
        qb = q*ab
        disc = qb*qb - bb*(q*q*aa-outer*outer)
        if disc < 0:
            continue
        root = isqrt(disc)
        root += root*root != disc
        lo, hi = ceil_div(-qb-root, bb), (-qb+root)//bb
        inner_disc = qb*qb-bb*(q*q*aa-inner*inner)
        intervals = [(lo,hi)]
        if inner > 0 and inner_disc >= 0:
            inner_root = isqrt(inner_disc)
            ilo = (-qb-inner_root)//bb + 1
            ihi = ceil_div(-qb+inner_root, bb)-1
            if ilo <= ihi:
                intervals = [(lo,min(hi,ilo-1)), (max(lo,ihi+1),hi)]
        for first,last in intervals:
            for r in range(first,last+1):
                if intersects_unit_circle(q,r,a,b,denominator):
                    offsets.add((q % shape[0], r % shape[1]))
    return np.array(sorted(offsets), dtype=np.int64).reshape(-1,2)


def neighbor_counts(grid, offsets, colors=range(5)):
    """Fast SEARCH queries only. Certification uses integer comparisons."""
    kernel = np.zeros(grid.shape, dtype=np.float64)
    kernel[(-offsets[:,0]) % grid.shape[0], (-offsets[:,1]) % grid.shape[1]] = 1
    transform = np.fft.rfft2(kernel)
    for color in colors:
        raw = np.fft.irfft2(np.fft.rfft2(grid == color)*transform, s=grid.shape)
        rounded = np.rint(raw)
        error = float(np.max(np.abs(raw-rounded)))
        if error > 1e-5 or rounded.min() < 0:
            raise ArithmeticError(f'FFT residual {error}')
        yield color, rounded.astype(np.int32)


def direct_conflicts(grid, offsets):
    """Full integer label audit, independent of FFT rounding and optimizers."""
    if grid.ndim != 2 or grid.dtype.kind not in 'iu' or grid.min() < 0 or grid.max() > 5:
        raise ValueError('Expected a 2D integer label grid, colors 0..5')
    h,w = grid.shape
    real = grid < 5
    total = 0
    for index,(i,j) in enumerate(offsets):
        i,j = int(i),int(j)
        opposite = ((-i)%h,(-j)%w)
        if (i,j) > opposite:
            continue
        count = int(np.count_nonzero((grid == np.roll(grid,(i,j),(0,1))) & real))
        total += count * (1 if (i,j) == opposite else 2)
        if index % 1000 == 0:
            print(f'audit offset {index}/{len(offsets)}, conflicts={total}', flush=True)
    return total


def extension_conflicts(parent, child, offsets, native=False):
    """Check every changed real cell against ALL neighbors in the final grid.

    Precondition: parent is valid for the same physical lattice. Each unchanged
    child cell is a subset of its parent, so any new conflict must have at
    least one endpoint in the changed set. This also covers simultaneous edits
    and real-to-real recolorings. The count is incident witnesses, not the
    full-grid directed edge count (which is unnecessary for the zero test).
    """
    if child.dtype.kind not in 'iu' or child.min()<0 or child.max()>5:
        raise ValueError('Invalid child labels')
    if child.shape[0]%parent.shape[0] or child.shape[1]%parent.shape[1]:
        raise ValueError('Partition is not a refinement')
    fh,fw=(child.shape[i]//parent.shape[i] for i in range(2))
    changed=[]
    fill_only=True
    parent_real=parent<5
    for s in range(fh):
        for t in range(fw):
            part=child[s::fh,t::fw]
            different=part != parent
            if np.any(different & parent_real):
                fill_only=False
            ci,cj=np.nonzero(different & (part < 5))
            changed.append(np.column_stack((ci*fh+s,cj*fw+t)))
    coords=np.concatenate(changed)
    if native and len(coords):
        import hashlib
        import json
        import shutil
        import subprocess
        import tempfile
        grouped=fill_only and fh*fw<=64
        source=Path(__file__).with_name('check_fill.cpp' if grouped else 'check_extension.cpp')
        compiler=shutil.which('c++')
        if compiler:
            cache=source.parent/'.cache'
            cache.mkdir(exist_ok=True)
            binary=cache/('check_extension_'+hashlib.sha256(source.read_bytes()).hexdigest()[:16])
            if not binary.exists():
                subprocess.run([compiler,'-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
            with tempfile.TemporaryDirectory(prefix='almost5-integer-audit-') as tmp:
                tmp=Path(tmp)
                np.asarray(offsets,dtype=np.int32).tofile(tmp/'offsets.bin')
                if grouped:
                    np.asarray(parent,dtype=np.uint8).tofile(tmp/'grid.bin')
                    paint=np.column_stack((coords,child[coords[:,0],coords[:,1]]))
                    np.asarray(paint,dtype=np.int32).tofile(tmp/'changed.bin')
                    command=[str(binary),str(tmp/'grid.bin'),str(parent.shape[0]),str(parent.shape[1]),
                             str(fh),str(fw),str(tmp/'offsets.bin'),str(tmp/'changed.bin'),str(tmp/'report.json')]
                else:
                    np.asarray(child,dtype=np.uint8).tofile(tmp/'grid.bin')
                    np.asarray(coords,dtype=np.int32).tofile(tmp/'changed.bin')
                    command=[str(binary),str(tmp/'grid.bin'),str(child.shape[0]),str(child.shape[1]),
                             str(tmp/'offsets.bin'),str(tmp/'changed.bin'),str(tmp/'report.json')]
                subprocess.run(command,check=True)
                report=json.loads((tmp/'report.json').read_text())
                if report['checked']!=len(coords):
                    raise ValueError('Native checker did not check every changed real cell')
                return report['conflict_witnesses'],len(coords)
    total=0
    h,w=child.shape
    for index,(i,j) in enumerate(coords):
        neighbors=child[(i+offsets[:,0])%h,(j+offsets[:,1])%w]
        total+=int(np.count_nonzero(neighbors==child[i,j]))
        if index%10000==0:
            print('extension audit',index,'/',len(coords),'conflict witnesses',total,flush=True)
    return total,len(coords)
