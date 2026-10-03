#!/usr/bin/env python3
"""Independent finite diagnostics for Rudvalis V15 audit and V16 profile proof.
Floating-point tests check algebra and normalization, not asymptotic theorems.
Run with OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python verify_rudvalis_profile.py.
"""
from __future__ import annotations
import itertools as it
import json, math, pathlib
import numpy as np
import scipy.linalg as la
import mpmath as mp

OUT=pathlib.Path(__file__).with_name('verification_results_rudvalis_profile.json')
mp.mp.dps=65

def root(n: int):
    z=mp.findroot(lambda z:(n-2)*z+mp.log(2*mp.exp(z)-1)-2j*mp.pi,
                  (2j*mp.pi/n,2j*mp.pi/n-mp.mpf(40)/n**3))
    lam=mp.exp(z)
    return complex(lam),complex(z)

def onecard(n:int,Q:int=1):
    L=np.zeros((n,n))
    for i in range(n):
        j=(i+1)%n; c=0.5/(Q if i==0 else 1)
        L[i,i]+=c; L[j,j]+=c; L[i,j]-=c; L[j,i]-=c
    return L

def transposition_matrix(states, index, a,b):
    cols=[]
    for x in states:
        y=tuple(b if v==a else a if v==b else v for v in x)
        cols.append(index[y])
    return np.eye(len(states))[cols]

def tuple_checks(n,r,Q,tau):
    xs=list(it.permutations(range(n),r)); ind={x:i for i,x in enumerate(xs)}
    d=len(xs); I=np.eye(d); ell=onecard(n,Q); p=la.expm(-tau*ell)
    L=np.zeros((d,d)); S=np.zeros((d,d)); swap=np.zeros((d,d))
    kill=np.zeros(d)
    for a in range(n):
        for b in range(a+1,n):
            U=transposition_matrix(xs,ind,a,b)
            if (a+1)%n==b: c=0.5/(Q if a==0 else 1)
            elif (b+1)%n==a: c=.5
            else: c=0
            L+=c*(I-U); S+=p[a,b]*(I-U)
            for ix,x in enumerate(xs):
                if a in x and b in x:
                    swap[ix,:]+=p[a,b]*(I-U)[ix,:]
    for ix,x in enumerate(xs):
        kill[ix]=sum(p[a,b] for a in x for b in x if a!=b)
    grid=list(it.product(range(n),repeat=r)); gi={x:i for i,x in enumerate(grid)}
    J=np.zeros((len(grid),d))
    for i,x in enumerate(xs):J[gi[x],i]=1
    A=np.zeros((len(grid),len(grid)))
    for a in range(r):
        factors=[(np.eye(n)-p if k==a else np.eye(n)) for k in range(r)]
        B=factors[0]
        for f in factors[1:]:B=np.kron(B,f)
        A+=B
    C=[]
    for a in range(r):
        ys=list(it.permutations(range(n),r-1)); yi={y:i for i,y in enumerate(ys)}
        D=np.zeros((len(ys),d))
        for i,x in enumerate(xs):D[yi[x[:a]+x[a+1:]],i]=1
        C.append(D)
    H=la.null_space(np.concatenate(C,axis=0))
    ev=la.eigvalsh(H.T@L@H)
    ge=la.eigvalsh(ell)[1:]
    prod=np.sort([sum(1-np.exp(-tau*ge[k]) for k in ks)
                  for ks in it.product(range(n-1),repeat=r)])
    rhs=(prod[:len(ev)]-r*(r-1)*p.max())/tau
    identity=float(np.max(np.abs(J.T@A@J-S+swap-np.diag(kill))))
    dyn=float(la.eigvalsh(tau*L-S)[0])
    eigmargin=float(np.min(ev-rhs)) if len(ev) else None
    assert identity<1e-9 and dyn>-1e-9 and (eigmargin is None or eigmargin>-1e-9)
    return dict(n=n,r=r,Q=Q,tau=tau,harmonic_dim=H.shape[1],
                identity_residual=identity,dynamical_min_eigenvalue=dyn,
                ordered_eigenvalue_margin=eigmargin)

def regular_checks(n,Q):
    xs=list(it.permutations(range(n))); ind={x:i for i,x in enumerate(xs)}
    d=len(xs); I=np.eye(d)
    As=[(I-transposition_matrix(xs,ind,i,(i+1)%n))/2 for i in range(n)]
    K=I.copy()
    for A in As:K=(I-A)@K
    L=sum(As[1:],As[0]/Q); c=(1+5/Q)/4
    endpoint=la.eigvalsh(c*(I-K.T@K)-K.T@L@K)[0]
    resolvent=la.solve(I+L/c,I,assume_a='pos')
    rm=la.eigvalsh(resolvent-K@K.T)[0]
    central=np.zeros((d,d))
    for a in range(n):
        for b in range(a+1,n):central+=I-transposition_matrix(xs,ind,a,b)
    path=sum(As[1:]); pathmin=la.eigvalsh(path-central/n**3)[0]
    ell=onecard(n,Q); ps=la.expm(-.4*ell)
    smooth=np.zeros((d,d))
    for a in range(n):
        for b in range(a+1,n):smooth+=ps[a,b]*(I-transposition_matrix(xs,ind,a,b))
    octmin=la.eigvalsh(.4*L-smooth)[0]
    sv=la.svdvals(K); le=la.eigvalsh(L)
    traces=[]
    for m in (1,2,4):
        act=float(np.sum(np.abs(np.linalg.matrix_power(K,m))**2)-1)
        sm=float(np.sum(sv**(2*m))-1)
        rt=float(np.sum((1+le/c)**(-m))-1)
        assert act<=sm+1e-7 and sm<=rt+1e-7
        traces.append(dict(m=m,actual_chi_square=act,singular_sum=sm,resolvent_trace=rt))
    assert min(endpoint,rm,pathmin,octmin)>-1e-8
    return dict(n=n,Q=Q,endpoint_min=float(endpoint),resolvent_min=float(rm),
                path_min=float(pathmin),octopus_min=float(octmin),traces=traces)

def bracket_check(n):
    xs=list(it.permutations(range(n))); ind={x:i for i,x in enumerate(xs)}
    a0=np.array([ind[x[1:]+x[:1]] for x in xs])
    a1=np.array([ind[x[1:-1]+x[:1]+x[-1:]] for x in xs])
    lam,z=root(n); phi=np.array([lam**(-j) for j in range(n-1)]+[1],complex)
    phi*=math.sqrt(n/float(np.vdot(phi,phi).real)); aa=phi.conj()
    X=np.asarray(xs,dtype=int); Z=(aa[X]@phi)/math.sqrt(n)
    def T(f):return (f[a0]+f[a1])/2
    eigres=np.max(np.abs(T(Z)-lam*Z))
    f1=np.abs(Z)**2; f2=Z**2
    for _ in range(n):f1=T(f1); f2=T(f2)
    Lam=lam**n
    b=f1-abs(Lam)**2*np.abs(Z)**2; c=f2-Lam**2*Z**2
    vb=n/(n-1); q=abs(np.sum(phi**2))**2/(n*(n-1))
    assert eigres<1e-8
    assert np.min(b)>-1e-8
    assert abs(b.mean()-(1-abs(Lam)**2)*vb)<1e-8
    assert abs(c.mean()-(1-Lam**2)*q)<1e-8
    return dict(n=n,eigenfunction_residual=float(eigres),
                n2_mean_b=float(n*n*b.mean()),
                n5over2_b_sd=float(n**2.5*np.sqrt(np.mean(abs(b-b.mean())**2))),
                n5over2_c_sd=float(n**2.5*np.sqrt(np.mean(abs(c-c.mean())**2))),
                mean_c_abs=float(abs(c.mean())),phi_square_sum_abs=float(abs(np.sum(phi**2))),
                scaled_gamma=float(-math.log(abs(Lam))*n*n))

def root_check(n):
    lam,z=root(n)
    gamma=-n*z.real
    ph=np.array([np.exp(-(j)*z) for j in range(n-1)]+[1],complex)
    ph*=math.sqrt(n/float(np.vdot(ph,ph).real))
    return dict(n=n,gamma_times_n2=gamma*n*n,
                gamma_relative_error=gamma*n*n/(4*math.pi**2)-1,
                phi_max=float(np.max(abs(ph))),
                phi_sum_abs=float(abs(ph.sum())),phi_square_sum_abs=float(abs((ph**2).sum())))

def graph_checks(n):
    total=0.; maxdeg=0
    for bits in it.product((0,1),repeat=n-1):
        deg=np.zeros(n,int); bottom=n-1
        for j,bit in enumerate(bits):
            deg[j]+=1;deg[bottom]+=1
            if bit:bottom=j
        total+=float(deg@deg);maxdeg=max(maxdeg,int(deg.max()))
    return dict(n=n,expected_degree_square_sum=total/2**(n-1),ratio=(total/2**(n-1))/n,max_degree=maxdeg)

def main() -> None:
    res={'description':'Fresh independent finite algebra and normalization diagnostics; not an asymptotic proof.',
         'tuple':[],'regular':[],'bracket':[],'roots':[],'graphs':[]}
    for n,r in ((4,2),(5,2),(5,3),(6,2),(7,2)):
        for Q in (1,7,31):
            for tau in (.2,2.):
                res['tuple'].append(tuple_checks(n,r,Q,tau))
    for n in (4,5,6):
        for Q in (1,19):res['regular'].append(regular_checks(n,Q))
    for n in (4,5,6,7,8):res['bracket'].append(bracket_check(n))
    for n in (10,30,100,1000,10000):res['roots'].append(root_check(n))
    for n in (4,8,12,16):res['graphs'].append(graph_checks(n))
    OUT.write_text(json.dumps(res,indent=2))
    print(json.dumps({'tuple_cases':len(res['tuple']),'regular_cases':len(res['regular']),
                      'max_identity_residual':max(x['identity_residual'] for x in res['tuple']),
                      'bracket':res['bracket'],'roots':res['roots'],'graphs':res['graphs']},indent=2))

if __name__ == "__main__":
    main()
