"""Exact algebraic whole-plane audit, including common-shift boundaries.

Coordinates are X=sqrt(3)*x,Y=y and belong to the specified real Q(b).
Coverage is checked on convex fan triangles; distance checks use entire
polygon components and all relevant periodic translations.
"""
if not __debug__:
    raise RuntimeError('Run without Python optimization (-O).')

import argparse
import hashlib
import json
import time
from functools import cmp_to_key,lru_cache
from pathlib import Path

import sympy as S
from sympy import QQ

ap=argparse.ArgumentParser();ap.add_argument('construction');ap.add_argument('--output',required=True);args=ap.parse_args()
raw=Path(args.construction).read_bytes();data=json.loads(raw);started=time.monotonic()
b=S.Symbol('b');coeff=list(map(S.Rational,data['field_polynomial']));poly=S.Poly.from_list(coeff,b)
assert poly.as_expr()==512*b**8-1536*b**7-1248*b**6+9072*b**5-4752*b**4-13608*b**3+12258*b**2-81*b-324
rootlo,roothi=map(S.Rational,data['field_root_interval'])
assert poly.count_roots(rootlo,roothi)==1
assert poly.count_roots(-S.oo,rootlo)==data['field_root_index']==3
assert poly.eval(rootlo)*poly.eval(roothi)<0
for _ in range(240):
    mid=(rootlo+roothi)/2
    if poly.eval(rootlo)*poly.eval(mid)<0:roothi=mid
    else:rootlo=mid
rootlo,roothi=QQ.convert(rootlo),QQ.convert(roothi)
K=QQ.algebraic_field(S.CRootOf(poly.as_expr(),3));zero,one,t=K.zero,K.one,K.unit
def parse(v):
    result=zero
    for c in v:result=result*t+K.convert(QQ(c))
    return result
def interval(x):
    lo=hi=QQ.zero
    for c in x.to_list():
        choices=[lo*rootlo,lo*roothi,hi*rootlo,hi*roothi]
        lo,hi=min(choices)+c,max(choices)+c
    return lo,hi
@lru_cache(maxsize=300000)
def sign(x):
    if not x:return 0
    lo,hi=interval(x)
    if lo>0:return 1
    if hi<0:return -1
    raise ArithmeticError('Unresolved algebraic sign')
def add(a,b):return a[0]+b[0],a[1]+b[1]
def sub(a,b):return a[0]-b[0],a[1]-b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def dot(a,b):return a[0]*b[0]/3+a[1]*b[1]
def edges(p):return list(zip(p,p[1:]+p[:1]))
def area(p):return sum((cross(a,b) for a,b in edges(p)),zero)
def shifted(p,v):return [add(x,v) for x in p]
@lru_cache(maxsize=50000)
def box(p):
    sort=lambda vals:sorted(vals,key=cmp_to_key(lambda x,y:sign(x-y)))
    x=sort([v[0] for v in p]);y=sort([v[1] for v in p])
    return x[0],x[-1],y[0],y[-1]
def separated(p,q):
    a,c=box(tuple(p)),box(tuple(q))
    if any(sign(a[2*k+1]-c[2*k])<=0 or sign(c[2*k+1]-a[2*k])<=0 for k in range(2)):return True
    for x,y in ((p,q),(q,p)):
        for v,w in edges(x):
            if all(sign(cross(sub(w,v),sub(z,v)))<=0 for z in y):return True
    return False
def far_box(p,q):
    a,c=box(tuple(p)),box(tuple(q))
    # Physical distance <=1 implies |delta X|<=sqrt(3)<2, |delta Y|<=1.
    return any(sign(a[2*k]-c[2*k+1]-R)>0 or sign(c[2*k]-a[2*k+1]-R)>0 for k,R in enumerate((2,1)))
def at_most(p,q,rr):return all(sign(dot(sub(v,w),sub(v,w))-rr)<=0 for v in p for w in q)
def at_least(p,q,rr):
    for x,y in ((p,q),(q,p)):
        for z in x:
            for v,w in edges(y):
                e,h=sub(w,v),sub(z,v);ee,projection=dot(e,e),dot(h,e)
                assert sign(ee)>0
                if sign(projection)<=0:residual=dot(h,h)-rr
                elif sign(projection-ee)>=0:residual=dot(sub(z,w),sub(z,w))-rr
                else:residual=dot(h,h)*ee-projection*projection-rr*ee
                if sign(residual)<0:return False
    return True
assert data['coordinate_metric']==['1/3','1']
L=[tuple(parse(x) for x in row) for row in data['lattice']]
assert L==[(2*t,zero),(t,t)]
low,high=map(parse,data['sixth_squared_interval'])
assert sign(low)>0 and sign(high-low)>0 and sign(one-high)>0
endpoint_poly=lambda D:16*D**8-112*D**7+355*D**6-541*D**5+322*D**4-106*D**3+55*D**2-25*D+4
assert endpoint_poly(high)==zero
components=[];triangles=[];colors=[]
direction=tuple(K.convert(QQ(x)) for x in data['boundary_direction'])
for row in data['polygons']:
    p=[tuple(parse(x) for x in v) for v in row['vertices']]
    if sign(area(p))<0:p.reverse()
    assert sign(area(p))>0 and len(set(p))==len(p)
    center=tuple(sum((v[k] for v in p),zero)/len(p) for k in range(2))
    pieces=[[center,v,w] for v,w in edges(p)]
    assert all(sign(area(q))>0 for q in pieces),'Nonpositive fan triangle'
    assert all(sign(cross(sub(w,v),direction))!=0 for v,w in edges(p))
    components.append(p);triangles.append(pieces);colors.append(row['color'])
assert len(components)==8 and sorted(colors)==[1,2,3,4,5,6,6,6]
assert sum((area(p) for p in components),zero)==2*cross(*L)
flat=box(tuple(v for p in components for v in p));sx,sy=flat[1]-flat[0],flat[3]-flat[2]
nb=next(n for n in range(30) if sign(n*t-sy-1)>0)
mb=next(m for m in range(30) if sign(2*m*t-nb*t-sx-2)>0)
counts=dict(triangle_coverage_pairs=0,component_distance_pairs=0,bbox_omissions=0)
endpoint_witness=False
for i,p in enumerate(components):
    for j in range(i,len(components)):
        for a in range(-mb,mb+1):
            for c in range(-nb,nb+1):
                shift=((2*a+c)*t,c*t);q=shifted(components[j],shift);same=i==j and a==c==0
                if far_box(p,q):counts['bbox_omissions']+=1;continue
                other=[shifted(tri,shift) for tri in triangles[j]]
                for k,first in enumerate(triangles[i]):
                    for ell,second in enumerate(other):
                        if same and ell<=k:continue
                        assert separated(first,second),('Coverage overlap',i,j,a,c,k,ell)
                        counts['triangle_coverage_pairs']+=1
                if colors[i]!=colors[j]:continue
                bound=low if colors[i]==6 else one
                if same:assert at_most(p,q,bound),('Component diameter',i)
                else:
                    bound=high if colors[i]==6 else one
                    assert at_least(p,q,bound),('Component separation',i,j,a,c)
                    if colors[i]==6 and any(dot(sub(v,w),sub(v,w))==high for v in p for w in q):endpoint_witness=True
                counts['component_distance_pairs']+=1
    print('exact component',i+1,'/ 8',flush=True)
assert endpoint_witness
result=dict(verified=True,scope='Whole-plane fixed coloring including algebraic upper endpoint',construction_sha256=hashlib.sha256(raw).hexdigest(),coordinate_metric=['1/3','1'],field_polynomial=str(poly.as_expr()),squared_interval_enclosures=[[str(x) for x in interval(v)] for v in (low,high)],translation_bounds=[mb,nb],checks=counts,upper_endpoint_attained_in_region_closures=True,boundary_rule='Color p by the unique component interior containing p+epsilon*(1,1/100) in scaled coordinates for all sufficiently small positive epsilon.',boundary_justification='If an assigned pair attained a forbidden distance, a common small shift would put both points in region interiors at that distance. Independent interior perturbations contradict the verified closed-region minimum or maximum.',seconds=time.monotonic()-started)
Path(args.output).write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k!='squared_interval_enclosures'},indent=2))
