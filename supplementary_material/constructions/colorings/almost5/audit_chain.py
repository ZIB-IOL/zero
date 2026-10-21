"""Independently audit a refinement chain from its root, using exact geometry."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import monotonic

import numpy as np
from geometry import read_lattice, unit_distance_offsets, direct_conflicts, extension_conflicts


def verify_chain(manifest,native=False):
    manifest=Path(manifest)
    spec=json.loads(manifest.read_text())
    lattice=manifest.parent/spec['lattice']
    parent=None
    reports=[]
    for entry in spec['steps']:
        start=monotonic()
        path=manifest.parent/entry['grid']
        raw=np.load(path)
        child=raw['coloring'] if 'coloring' in raw else raw['fixed_coloring']
        vectors,shape=read_lattice(lattice,entry['refinement'])
        if child.shape!=shape:
            raise ValueError('Grid shape disagrees with physical lattice subdivisions')
        offsets=unit_distance_offsets(vectors,shape)
        if parent is None:
            conflicts=direct_conflicts(child,offsets)
            changed=None
        else:
            conflicts,changed=extension_conflicts(parent,child,offsets,native=native)
        removed=int((child==5).sum())
        report={'grid':entry['grid'],'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
                'refinement':entry['refinement'],'shape':shape,'exact_offsets':len(offsets),
                'removed':removed,'removed_fraction':str(Fraction(removed,child.size)),
                'removed_percent':100*removed/child.size,'changed_real_cells_checked':changed,
                'conflict_witnesses':conflicts,'verified':conflicts==0,'seconds':monotonic()-start}
        reports.append(report)
        print(json.dumps(report),flush=True)
        if conflicts:
            raise ValueError('Chain step contains conflicts')
        parent=child
    return {'verified':True,'geometry':'exact rational closed cells',
            'lattice_sha256':hashlib.sha256(lattice.read_bytes()).hexdigest(),'steps':reports}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manifest',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--native',action='store_true',help='Use the independent C++ integer comparison loop')
    args=parser.parse_args()
    args.output.write_text(json.dumps(verify_chain(args.manifest,native=args.native),indent=2)+'\n')
