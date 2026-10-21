"""Independently verify a parent certificate plus sparse fine-cell additions."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from time import monotonic

import numpy as np

from audit_chain import verify_chain
from geometry import read_lattice,unit_distance_offsets
from sparse_geometry import sparse_fill_conflicts


def verify_sparse(manifest,native=False):
    manifest=Path(manifest)
    spec=json.loads(manifest.read_text())
    if spec['format']!='sparse_fill_v1':
        raise ValueError('Unsupported sparse certificate')
    chain_path=manifest.parent/spec['parent_chain']
    parent_audit=verify_chain(chain_path,native=native)
    chain=json.loads(chain_path.read_text())
    parent_entry=chain['steps'][-1]
    parent_path=chain_path.parent/parent_entry['grid']
    raw=np.load(parent_path)
    parent=raw['coloring'] if 'coloring' in raw else raw['fixed_coloring']
    if hashlib.sha256(parent_path.read_bytes()).hexdigest()!=parent_audit['steps'][-1]['sha256']:
        raise ValueError('Parent changed after ancestry audit')
    refinement=spec['refinement']
    if type(refinement) is not int or refinement<parent_entry['refinement'] or refinement%parent_entry['refinement']:
        raise ValueError('Invalid sparse refinement')
    factor=refinement//parent_entry['refinement']
    paint_path=manifest.parent/spec['paint']
    paint=np.load(paint_path)['paint']
    lattice=chain_path.parent/chain['lattice']
    vectors,shape=read_lattice(lattice,refinement)
    if shape!=tuple(n*factor for n in parent.shape):
        raise ValueError('Fine shape disagrees with physical lattice')
    if hashlib.sha256(lattice.read_bytes()).hexdigest()!=parent_audit['lattice_sha256']:
        raise ValueError('Lattice changed after ancestry audit')
    start=monotonic()
    offsets=unit_distance_offsets(vectors,shape)
    conflicts,checked=sparse_fill_conflicts(parent,paint,offsets,factor,native=native)
    total=shape[0]*shape[1]
    removed=int(np.count_nonzero(parent==5))*factor*factor-len(paint)
    extension={'paint_sha256':hashlib.sha256(paint_path.read_bytes()).hexdigest(),
               'parent_sha256':parent_audit['steps'][-1]['sha256'],'refinement':refinement,'shape':shape,
               'exact_offsets':len(offsets),'added_cells':len(paint),'changed_real_cells_checked':checked,
               'removed':removed,'total_fine_cells':total,'removed_fraction':str(Fraction(removed,total)),
               'removed_percent':100*removed/total,'conflict_witnesses':conflicts,'verified':conflicts==0,
               'seconds':monotonic()-start}
    print(json.dumps(extension),flush=True)
    if conflicts:
        raise ValueError('Sparse certificate contains conflicts')
    return {'verified':True,'format':'sparse_fill_v1','geometry':'exact rational closed cells',
            'parent_audit':parent_audit,'extension':extension}


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('manifest',type=Path)
    p.add_argument('--native',action='store_true')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    result=verify_sparse(args.manifest,native=args.native)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
