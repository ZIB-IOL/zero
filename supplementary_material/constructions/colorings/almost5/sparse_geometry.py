"""Exact fill checks for a dense parent plus sparse fine-cell additions.

Unchanged real parent cells are retained in their entirety. The remaining
parent cells are subdivided conceptually; only newly colored children are
stored. This defines the full coloring without allocating the fine grid.
"""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

import numpy as np


def validate_paint(parent,paint,factor):
    if parent.ndim!=2 or parent.dtype.kind not in 'iu' or parent.min()<0 or parent.max()>5:
        raise ValueError('Invalid parent labels')
    if not isinstance(factor,int) or factor<1:
        raise ValueError('Expected a positive integer refinement factor')
    h,w=parent.shape
    H,W=h*factor,w*factor
    if max(H,W)>np.iinfo(np.int32).max//2 or H*W>np.iinfo(np.int64).max:
        raise ValueError('Unsupported fine dimensions')
    if paint.ndim!=2 or paint.shape[1]!=3 or paint.dtype.kind not in 'iu':
        raise ValueError('Expected integer (row,column,color) triples')
    if len(paint):
        if (paint<0).any() or (paint[:,0]>=H).any() or (paint[:,1]>=W).any() or (paint[:,2]>=5).any():
            raise ValueError('Paint outside fine grid or invalid real color')
        ids=paint[:,0].astype(np.int64)*W+paint[:,1].astype(np.int64)
        if len(np.unique(ids))!=len(ids):
            raise ValueError('Duplicate paint coordinate')
        if (parent[paint[:,0]//factor,paint[:,1]//factor]!=5).any():
            raise ValueError('Sparse extension can only fill removed parent cells')
    return H,W


def sparse_fill_conflicts(parent,paint,offsets,factor,native=False):
    """Count all incident conflicts, exactly as a fully expanded grid would."""
    H,W=validate_paint(parent,paint,factor)
    if offsets.ndim!=2 or offsets.shape[1]!=2 or offsets.dtype.kind not in 'iu':
        raise ValueError('Invalid offset array')
    if len(offsets) and ((offsets<0).any() or (offsets[:,0]>=H).any() or (offsets[:,1]>=W).any()):
        raise ValueError('Offsets must be normalized to the fine torus')
    if len(np.unique(offsets,axis=0))!=len(offsets):
        raise ValueError('Expected distinct geometric offsets')
    if not len(paint):
        return 0,0
    if native and factor*factor<=64 and shutil.which('c++'):
        source=Path(__file__).with_name('check_fill.cpp')
        cache=source.parent/'.cache'
        cache.mkdir(exist_ok=True)
        binary=cache/('check_extension_'+hashlib.sha256(source.read_bytes()).hexdigest()[:16])
        if not binary.exists():
            subprocess.run([shutil.which('c++'),'-std=c++17','-O3',str(source),'-o',str(binary)],check=True)
        with tempfile.TemporaryDirectory(prefix='almost5-sparse-audit-') as tmp:
            tmp=Path(tmp)
            np.asarray(parent,dtype=np.uint8).tofile(tmp/'parent.bin')
            np.asarray(paint,dtype=np.int32).tofile(tmp/'paint.bin')
            np.asarray(offsets,dtype=np.int32).tofile(tmp/'offsets.bin')
            subprocess.run([str(binary),str(tmp/'parent.bin'),str(parent.shape[0]),str(parent.shape[1]),
                            str(factor),str(factor),str(tmp/'offsets.bin'),str(tmp/'paint.bin'),str(tmp/'report.json')],check=True)
            report=json.loads((tmp/'report.json').read_text())
            if report['checked']!=len(paint):
                raise ValueError('Incomplete native sparse audit')
            return report['conflict_witnesses'],len(paint)
    # Independent reference: query every fine displacement of every new cell.
    lookup={int(i)*W+int(j):int(c) for i,j,c in paint}
    conflicts=0
    for index,(i,j,c) in enumerate(paint):
        ni=(int(i)+offsets[:,0])%H
        nj=(int(j)+offsets[:,1])%W
        labels=parent[ni//factor,nj//factor]
        conflicts+=int(np.count_nonzero(labels==c))
        removed=labels==5
        conflicts+=sum(lookup.get(int(q)*W+int(r),5)==c for q,r in zip(ni[removed],nj[removed]))
        if index%10000==0:
            print('sparse reference audit',index,'/',len(paint),'conflicts',conflicts,flush=True)
    return conflicts,len(paint)
