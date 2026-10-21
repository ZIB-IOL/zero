"""Exact independent reconstruction/check of the published second construction.
Requires SymPy; all geometry is in Q(t), with t the specified isolated real root.
No optimizer, floating-point geometry, or imported family.py is used.
"""
if not __debug__:
    raise RuntimeError("Run this proof checker without Python optimization (-O).")

import argparse,json,time,hashlib
from functools import lru_cache
from fractions import Fraction
import sympy as S
from sympy import QQ

T=S.Symbol('t'); POLY=T**8-39*T**6+400*T**4-279*T**2+49
K=QQ.algebraic_field(S.CRootOf(POLY,3)); zero=K.zero; one=K.one; t=K.unit

# Rational root isolation and uniqueness (not a numerical root assertion).
rootlo=QQ(6570609747898614,10**16); roothi=QQ(6570609747898616,10**16)
assert S.Poly(POLY,T).count_roots(S.Rational(str(rootlo)),S.Rational(str(roothi)))==1
assert S.Poly(POLY,T).count_roots(-S.oo,S.Rational(str(rootlo)))==3
assert S.Poly(POLY,T).eval(S.Rational(str(rootlo)))<0
assert S.Poly(POLY,T).eval(S.Rational(str(roothi)))>0
for _ in range(220):
    mid=(rootlo+roothi)/2
    # Maintain a sign-changing isolating interval.
    if S.Poly(POLY,T).eval(S.Rational(str(mid)))<0: rootlo=mid
    else: roothi=mid

def interval(x):
    lo=hi=QQ.zero
    for c in x.to_list():
        vals=(lo*rootlo,lo*roothi,hi*rootlo,hi*roothi)
        lo,hi=min(vals)+c,max(vals)+c
    return lo,hi

@lru_cache(maxsize=300000)
def sign(x):
    if not x:return 0
    lo,hi=interval(x)
    if lo>0:return 1
    if hi<0:return -1
    raise ArithmeticError('Unresolved nonzero algebraic sign; increase root precision')

assert sign(8*t**7-234*t**5+1600*t**3-558*t)>0
s=-(t**4+18*t*t-7)/(5*t**3-3*t)
assert s*s==3*one and sign(s)>0
assert t**4+5*s*t**3+18*t*t-3*s*t-7==zero
r=s/2-t
c=(s+2*t)*(4-2*s*t-t*t)/(4*s)
assert c*c==(3-s*t-t*t)*(s+2*t)**2/(4*(1+(s+2*t)**2)) and sign(c)>0
u=one/4-c; v=3*s/4-t/2-c/(s+2*t); x=one/14;y=2*s/7

def pt(a,b):return (K.convert(a),K.convert(b))
def add(a,b):return (a[0]+b[0],a[1]+b[1])
def sub(a,b):return (a[0]-b[0],a[1]-b[1])
def dot(a,b):return a[0]*b[0]+a[1]*b[1]
def cross(a,b):return a[0]*b[1]-a[1]*b[0]
def edges(p):return list(zip(p,p[1:]+p[:1]))
def area2(p):return sum((cross(a,b) for a,b in edges(p)),zero)

# Direct transcription of the paper's labelled vertices and eight placed tiles.
p=dict(A=pt(0,0),AA=pt(1,0),PA=pt(one/2,s/2),PB=pt(-one/2,s/2),PC=pt(3*one/2,s/2),
 B=pt(0,s),C=pt(1,s),X=pt(one/2,r),XX=pt(one/2-r,0),XXX=pt(one/2+r,0),MAA=pt(one/2,-r),
 Y=pt(-one/2,s-r),YY=pt(-one/2+r,s),Z=pt(3*one/2,s-r),ZZ=pt(3*one/2-r,s),
 I1=pt(u,v),I2=pt(-u,s-v),I3=pt(1+u,s-v),I4=pt(1-u,v),
 M=pt(one/2-x,s-y),N=pt(one/2+x,s-y),NN=pt(one/2-x,s+y),MM=pt(one/2+x,s+y),
 MMM=pt(-one/2+x,y),NNN=pt(3*one/2-x,y),I5=pt(1+u,s+v),I6=pt(-u,s+v),PD=pt(one/2,3*s/2))
spec=[(1,'Z ZZ C N I3 I4 PC',pt(-1,-s)),(2,'Y YY B M I2 I1 PB',pt(1,-s)),
 (6,'PA I2 M N I3',pt(1,-s)),(5,'M N C MM NN B',pt(1,-s)),
 (6,'I5 MM NN I6 PD',pt(1,-s)),(3,'X XX A MMM I1 I2 PA',pt(0,0)),
 (4,'X XXX AA NNN I4 I3 PA',pt(0,0)),(6,'X XX MAA XXX',pt(0,0))]
polys=[]
for color,names,shift in spec:
    q=[add(p[n],shift) for n in names.split()]
    if sign(area2(q))<0:q.reverse()
    assert sign(area2(q))>0
    # Convexity, tested with ALL vertices against every oriented edge.
    assert all(sign(cross(sub(b,a),sub(w,a)))>=0 for a,b in edges(q) for w in q)
    polys.append(q)
assert sum((area2(q) for q in polys),zero)==4*s

@lru_cache(maxsize=None)
def bbox(p):
    from functools import cmp_to_key
    xs=sorted((v[0] for v in p),key=cmp_to_key(lambda a,b:sign(a-b)))
    ys=sorted((v[1] for v in p),key=cmp_to_key(lambda a,b:sign(a-b)))
    return xs[0],xs[-1],ys[0],ys[-1]

def separated(p,q):
    """True iff convex polygon interiors are disjoint (SAT)."""
    pbox=bbox(tuple(p));qbox=bbox(tuple(q))
    if any(sign(pbox[2*k+1]-qbox[2*k])<=0 or sign(qbox[2*k+1]-pbox[2*k])<=0 for k in range(2)):return True
    for a,b in edges(p):
        e=sub(b,a)
        if all(sign(cross(e,sub(w,a)))<=0 for w in q):return True
    for a,b in edges(q):
        e=sub(b,a)
        if all(sign(cross(e,sub(w,a)))<=0 for w in p):return True
    return False

def at_least(p,q,rr):
    if not separated(p,q):return False
    for a,b in ((p,q),(q,p)):
        for z in a:
            for v,w in edges(b):
                e=sub(w,v); h=sub(z,v); ee=dot(e,e); proj=dot(h,e)
                if sign(proj)<=0:val=dot(h,h)-rr
                elif sign(proj-ee)>=0:val=dot(sub(z,w),sub(z,w))-rr
                else:val=cross(e,h)**2-rr*ee
                if sign(val)<0:return False
    return True

def at_most(p,q,rr):return all(sign(dot(sub(a,b),sub(a,b))-rr)<=0 for a in p for b in q)

def far_bbox(p,q):
    a=bbox(tuple(p));b=bbox(tuple(q))
    return any(sign(a[2*k]-b[2*k+1]-1)>0 or sign(b[2*k]-a[2*k+1]-1)>0 for k in range(2))

def serial(x):return [str(c) for c in x.to_list()]

def run():
    start=time.monotonic();flat=[v for p in polys for v in p];box=bbox(tuple(flat))
    sx=box[1]-box[0];sy=box[3]-box[2]
    nb=next(n for n in range(20) if sign(n*s-sy-1)>0)
    mb=next(m for m in range(20) if sign(2*m-nb-sx-1)>0)
    # Outside these bounds: |n|s>sy+1 or |2m+n|>sx+1.
    counts=dict(coverage=0,unit=0,red=0,bbox_omissions=0)
    for i,a in enumerate(polys):
        for j in range(i,len(polys)):
            for m in range(-mb,mb+1):
                for n in range(-nb,nb+1):
                    shift=pt(2*m+n,n*s);b=[add(v,shift) for v in polys[j]]
                    same=i==j and m==n==0
                    if far_bbox(a,b):counts['bbox_omissions']+=1;continue
                    if not same:
                        assert separated(a,b),('coverage overlap',i,j,m,n)
                        counts['coverage']+=1
                    if spec[i][0]!=spec[j][0]:continue
                    if spec[i][0]<6:
                        assert (at_most(a,b,one) if same else at_least(a,b,one)),('unit',i,j,m,n)
                        counts['unit']+=1
                    else:
                        assert (at_most(a,b,(2*r)**2) if same else at_least(a,b,t*t)),('red',i,j,m,n)
                        counts['red']+=1
        print('exact tile',i+1,'/ 8',flush=True)
    # A nearest red vertex pair attains t in the closures. Their open
    # interiors attain distances immediately above t, proving this fixed
    # reference cannot extend its upper gap beyond t.
    nearest=[add(v,pt(-1,-s)) for v in polys[4]]
    assert any(dot(sub(a,b),sub(a,b))==t*t for a in polys[2] for b in nearest)
    # A common infinitesimal shift direction not parallel to any edge.
    direction=pt(1,one/100)
    assert all(sign(cross(sub(b,a),direction))!=0 for p in polys for a,b in edges(p))
    return dict(verified=True,scope='complete published fixed construction, including both distance endpoints',field_polynomial=str(POLY),root_interval=[str(rootlo),str(roothi)],dmin_interval=list(map(str,interval(2*r))),dmax_interval=list(map(str,interval(t))),lattice=[[["2"],["0"]],[["1"],serial(s)]],polygons=[dict(color=sp[0],vertices=[[serial(x),serial(y)] for x,y in p]) for sp,p in zip(spec,polys)],boundary='Assign p the interior tile containing p + epsilon*(1,1/100) for all sufficiently small epsilon>0.',translation_bounds=[mb,nb],checks=counts,seconds=time.monotonic()-start)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args()
    result=run()
    with open(args.output,'w') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in ('polygons','root_interval','dmin_interval','dmax_interval')},indent=2))
