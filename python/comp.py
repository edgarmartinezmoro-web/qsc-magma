from lib import *
import sys, pickle, time
def ge(G,p,D):
    """True iff min distance of rowspace(G) >= D"""
    n=G.shape[1]; k=G.shape[0]
    if k==0: return True
    if D<=1: return True
    H=nullspace(G,p)
    if H.shape[0]==0: return D<=1
    for w in range(1,D):
        for S in itertools.combinations(range(n),w):
            if rank(H[:,S],p)<w: return False
    return True
def search(n,p,Kp,D1,D2,minord):
    found=[]; stats={}
    for lam in range(1,p):
        f=[(-lam)%p]+[0]*(n-1)+[1]
        divs=divisors(f,p)
        divs.sort(key=len)
        G={tuple(g):genmat(g,n,p) for g in divs}
        okU={}; okG={}
        for u in divs:
            ku=n-(len(u)-1)
            if len(u)-1 < Kp+1: continue
            key=tuple(u)
            okU[key]=ge(nullspace(G[key],p),p,D2) if ku<n else False
        cnt=0
        for u in divs:
            if not okU.get(tuple(u)): continue
            du=len(u)-1
            for g3 in divs:
                d3=len(g3)-1
                if du-d3!=Kp or not divides(g3,u,p): continue
                for g1 in divs:
                    if len(g1)>len(g3) or not divides(g1,g3,p): continue
                    if len(g1)==len(g3): continue
                    t,_=pdivmod(g3,g1,p); o=order(t,p)
                    if o is None or o<minord: continue
                    key=tuple(g1)
                    if key not in okG: okG[key]=ge(G[key],p,D1)
                    if okG[key]:
                        found.append((lam,g1,g3,u,o)); break
                if found: return found
    return found
if __name__=="__main__":
    n,p,Kp,D1,D2,mo=map(int,sys.argv[1:7])
    t=time.time(); r=search(n,p,Kp,D1,D2,mo)
    print(f"n={n} q={p} K'={Kp} target(d1,d2)=({D1},{D2}) ord>={mo}: ", "FOUND "+str(r[0]) if r else "NONE", f"[{time.time()-t:.0f}s]")
