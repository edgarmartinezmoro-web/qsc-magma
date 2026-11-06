import itertools, numpy as np, sympy as sp
from math import comb
x=sp.symbols('x')
def parse(s,p):
    d={'A':10,'B':11,'C':12,'D':13,'E':14,'F':15,'G':16}
    return [ (d[c] if c in d else int(c))%p for c in s]      # increasing degree
def pmod(a,p): 
    a=[c%p for c in a]
    while len(a)>1 and a[-1]==0: a.pop()
    return a
def pdivmod(a,b,p):
    a=pmod(a,p); b=pmod(b,p); inv=pow(b[-1],-1,p); q=[0]*max(1,len(a)-len(b)+1)
    a=a[:]
    for i in range(len(a)-len(b),-1,-1):
        c=a[i+len(b)-1]*inv%p; q[i]=c
        for j,bj in enumerate(b): a[i+j]=(a[i+j]-c*bj)%p
    return pmod(q,p), pmod(a[:len(b)-1] or [0],p)
def pmul(a,b,p):
    r=[0]*(len(a)+len(b)-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b): r[i+j]=(r[i+j]+ai*bj)%p
    return pmod(r,p)
def divides(g,f,p): return pdivmod(f,g,p)[1]==[0]
def genmat(g,n,p):
    d=len(g)-1;k=n-d;G=np.zeros((k,n),dtype=np.int64)
    for i in range(k): G[i,i:i+d+1]=g
    return G
def rref(M,p):
    M=M.copy()%p; r=0; rows,cols=M.shape; piv=[]
    for c in range(cols):
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        i=r+nz[0]; M[[r,i]]=M[[i,r]]; M[r]=(M[r]*pow(int(M[r,c]),-1,p))%p
        col=M[:,c].copy(); col[r]=0
        M=(M-np.outer(col,M[r]))%p
        piv.append(c); r+=1
        if r==rows: break
    return M[:r],piv
def rank(M,p): return len(rref(M,p)[1])
def nullspace(G,p):
    R,piv=rref(G,p); n=G.shape[1]; free=[c for c in range(n) if c not in piv]
    N=np.zeros((len(free),n),dtype=np.int64)
    for t,fc in enumerate(free):
        N[t,fc]=1
        for i,pc in enumerate(piv): N[t,pc]=(-R[i,fc])%p
    return N
def contains(A,B,p):  # rowspace(B) ⊆ rowspace(A)
    if B.shape[0]==0: return True
    return rank(np.vstack([A,B]),p)==rank(A,p)
def mindist(G,p,cap=None):
    """min distance = min #dependent columns of parity-check H"""
    n=G.shape[1]; k=G.shape[0]
    if k==0: return None
    if p**k<=2*10**5:
        M=np.array(list(itertools.product(range(p),repeat=k)),dtype=np.int64)[1:]
        return int(((M@G)%p!=0).sum(1).min())
    H=nullspace(G,p); r=H.shape[0]
    for w in range(1,n-k+2):
        if cap and w>cap: return None
        for S in itertools.combinations(range(n),w):
            if rank(H[:,S],p)<w: return w
    return n-k+1
def factors(f,p):
    P=sp.Poly(list(reversed(f)),x,modulus=p)
    out=[]
    for fac,e in P.factor_list()[1]:
        c=[int(a)%p for a in reversed(fac.all_coeffs())]
        out+=[c]*e
    return out
def divisors(f,p):
    fs=factors(f,p); seen={}
    for mask in itertools.product([0,1],repeat=len(fs)):
        g=[1]
        for m,fc in zip(mask,fs):
            if m: g=pmul(g,fc,p)
        inv=pow(g[-1],-1,p); g=[c*inv%p for c in g]
        seen[tuple(g)]=g
    return list(seen.values())
def order(k,p,maxe=10**6):
    # smallest e with k | x^e -1 ; k(0)!=0
    if len(k)==1: return 1
    r=[0,1]; e=1
    while e<=maxe:
        if pdivmod(r,k,p)[1]==[1]: return e
        r=pdivmod(pmul(r,[0,1],p),k,p)[1]; e+=1
    return None

def ppowmod(e,m,p):
    r=[1]; b=[0,1]
    while e:
        if e&1: r=pdivmod(pmul(r,b,p),m,p)[1]
        b=pdivmod(pmul(b,b,p),m,p)[1]; e>>=1
    return r
def order_irred(phi,p):
    m=len(phi)-1; N=p**m-1; o=N
    for pr,_ in sp.factorint(N).items():
        while o%pr==0 and ppowmod(o//pr,phi,p)==[1]: o//=pr
    return o
def order(k,p,maxe=None):
    if len(k)==1: return 1
    P=sp.Poly(list(reversed(k)),x,modulus=p)
    from math import lcm, ceil, log
    L=1
    for fac,e in P.factor_list()[1]:
        c=[int(a)%p for a in reversed(fac.all_coeffs())]
        if len(c)==2 and c[0]==0: return None
        o=order_irred(c,p); t=0
        while p**t<e: t+=1
        L=lcm(L,o*p**t)
    return L
