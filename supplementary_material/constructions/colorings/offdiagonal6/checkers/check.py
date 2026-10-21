"""Standalone exact checker. Python standard library only; no solver needed.

Coverage: convex full-dimensional tiles have disjoint interiors on the torus
and areas sum to the lattice determinant. Local finiteness implies their
closures cover every point. The shared infinitesimal shift assigns boundaries.

Distances: the Minkowski difference of each convex pair has a connected norm
image. Vertex maxima and exact edge projections give its entire distance range.
A common infinitesimal shift cannot make a forbidden endpoint occur: it would
move both points into interiors while leaving their difference unchanged,
contradicting extremality. This permits non-strict endpoint inequalities.
"""
if not __debug__:
    raise RuntimeError("Run this proof checker without Python optimization (-O).")

import argparse,hashlib,json,time
from fractions import Fraction as F
from pathlib import Path

def sub(a,b):return a[0]-b[0],a[1]-b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def edges(p):return zip(p,p[1:]+p[:1])
def area2(p):return sum(cross(a,b) for a,b in edges(p))

def hull(points):
    pts=sorted(set(points));result=[]
    for seq in (pts,pts[::-1]):
        h=[]
        for p in seq:
            while len(h)>1 and cross(sub(h[-1],h[-2]),sub(p,h[-1]))<=0:h.pop()
            h.append(p)
        result+=h[:-1]
    return result

def range2(h):
    mx=max(dot(v,v) for v in h)
    if all(cross(sub(b,a),(-a[0],-a[1]))>=0 for a,b in edges(h)):return F(0),mx
    vals=[]
    for a,b in edges(h):
        e=sub(b,a);ee=dot(e,e);u=-dot(a,e)
        if u<=0:vals.append(dot(a,a))
        elif u>=ee:vals.append(dot(b,b))
        else:vals.append(cross(a,e)**2/ee)
    return min(vals),mx

def box(p):return min(v[0] for v in p),max(v[0] for v in p),min(v[1] for v in p),max(v[1] for v in p)
def ceil(x):return -(-x.numerator//x.denominator)

def verify(path):
    start=time.monotonic();data=json.loads(Path(path).read_text())
    L=[[F(x) for x in row] for row in data['lattice']]
    assert len(L)==2 and all(len(v)==2 for v in L)
    assert L[0][1]==0 and L[0][0]>0 and L[1][1]>0,'Expected triangular lattice with positive diagonal'
    ds=list(map(F,data['forbidden_distances']));assert len(ds)==6 and ds[:5]==[1]*5 and ds[5]>0
    low,high=map(F,data['sixth_forbidden_interval']);assert 0<low<=ds[5]<=high
    direction=tuple(map(F,data['boundary_direction']));assert len(direction)==2
    ps=[];cs=[]
    for row in data['polygons']:
        p=[tuple(map(F,v)) for v in row['vertices']];assert all(len(v)==2 for v in p) and len(p)>=3
        if area2(p)<0:p.reverse()
        assert area2(p)>0
        assert len(set(p))==len(p)
        assert all(cross(sub(b,a),sub(v,a))>=0 for a,b in edges(p) for v in p),'Non-convex polygon'
        assert all(cross(sub(b,a),direction)!=0 for a,b in edges(p)),'Boundary direction parallel to edge'
        assert 1<=row['color']<=6;ps.append(p);cs.append(row['color'])
    assert set(cs)==set(range(1,7))
    det=cross(*L);assert sum(area2(p) for p in ps)==2*det,'Area mismatch'
    bbox=box([v for p in ps for v in p]);R=max(ds+[high]);sx=bbox[1]-bbox[0];sy=bbox[3]-bbox[2]
    nb=ceil((sy+R)/L[1][1]);mb=ceil((sx+R+abs(L[1][0])*nb)/L[0][0])
    counts=dict(coverage_pairs=0,distance_pairs=0,bbox_skips=0);red_max=F(0);red_min=None
    for i,a in enumerate(ps):
        abox=box(a)
        for j in range(i,len(ps)):
            for m in range(-mb,mb+1):
                for n in range(-nb,nb+1):
                    sh=(m*L[0][0]+n*L[1][0],n*L[1][1]);b=[(v[0]+sh[0],v[1]+sh[1]) for v in ps[j]];bb=box(b)
                    if any(abox[2*k]-bb[2*k+1]>R or bb[2*k]-abox[2*k+1]>R for k in range(2)):
                        counts['bbox_skips']+=1;continue
                    h=hull([sub(y,x) for x in a for y in b]);same=i==j and m==n==0
                    if not same:
                        assert not all(cross(sub(v,u),(-u[0],-u[1]))>0 for u,v in edges(h)),('Interior overlap',i,j,m,n)
                        counts['coverage_pairs']+=1
                    if cs[i]!=cs[j]:continue
                    lo,hi=range2(h);d=ds[cs[i]-1]
                    assert hi<=d*d or lo>=d*d,('Forbidden distance',i,j,m,n,str(lo),str(hi))
                    if cs[i]==6:
                        assert hi<=low*low or lo>=high*high,('Forbidden interval',i,j,m,n)
                        if same:red_max=max(red_max,hi)
                        else:red_min=lo if red_min is None else min(red_min,lo)
                    counts['distance_pairs']+=1
    # All omitted pairs have distance>R, so cannot change this minimum if <=R².
    assert red_min is not None and red_min<=R*R
    return dict(verified=True,construction_sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),forbidden_distances=[str(d) for d in ds],sixth_forbidden_interval=[str(low),str(high)],boundary_direction=list(map(str,direction)),coverage='all plane points, including assigned boundaries',whole_regions=True,translation_bounds=[mb,nb],exact_area=str(det),red_diameter_squared=str(red_max),red_separation_squared=str(red_min),checks=counts,seconds=time.monotonic()-start)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('construction');ap.add_argument('--output');args=ap.parse_args()
    result=verify(args.construction)
    if args.output:Path(args.output).write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
