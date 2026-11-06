import itertools, numpy as np, sympy as sp
from sympy import Poly, symbols, GF
x=symbols('x')
def factors(q,n,lam):
    P=Poly(x**n-lam,x,modulus=q)
    fl=P.factor_list()[1]
    out=[]
    for f,e in fl: out+= [f]*e
    return out
def coeffs(p,q,n):
    c=[int(a)%q for a in reversed(p.all_coeffs())]
    return c+[0]*(n-len(c))
def gens(q,n,lam):
    fs=factors(q,n,lam); seen={}
    for mask in itertools.product([0,1],repeat=len(fs)):
        g=Poly(1,x,modulus=q)
        for m,f in zip(mask,fs):
            if m: g=g*f
        key=tuple(coeffs(g,q,n+1))
        seen[key]=g
    return list(seen.values())
def mind(G,q):
    k=G.shape[0]
    if k==0: return None
    best=10**9
    M=np.array(list(itertools.product(range(q),repeat=k)))[1:]
    W=(M@G)%q
    return int((W!=0).sum(1).min())
def genmat(g,q,n):
    d=g.degree(); k=n-d; c=coeffs(g,q,n)
    G=np.zeros((k,n),int)
    for i in range(k): G[i,i:i+d+1]=c[:d+1]
    return G
def best(q,n,lam,maxk):
    res={}
    for g in gens(q,n,lam):
        k=n-g.degree()
        if 1<=k<=maxk and q**k<=3**11:
            d=mind(genmat(g,q,n),q); res[k]=max(res.get(k,0),d)
    return res
import sys
if __name__!="__main__": sys.argv=["x","5","1","2","1"]
q=int(sys.argv[1]); lams=[int(a) for a in sys.argv[2].split(',')]
for n in range(int(sys.argv[3]),int(sys.argv[4])+1):
    if n%q==0: continue
    rs={l:best(q,n,l,n-1) for l in lams}
    ks=sorted(set().union(*[set(r) for r in rs.values()]))
    diff=[(k,{l:rs[l].get(k,'-') for l in lams}) for k in ks if len(set(rs[l].get(k,-1) for l in lams))>1]
    better=[(k,v) for k,v in diff if all(isinstance(v[l],int) for l in lams) and max(v[l] for l in lams[1:])>v[lams[0]]]
    print(n, 'constacyclic strictly better (both exist):', better)
