"""synthetic folded meshes + fsaverage5 lh.pial for FlipOut reference answers"""
import json, sys, numpy as np, os
def icosphere(k):
    t=(1+5**.5)/2
    V=[[-1,t,0],[1,t,0],[-1,-t,0],[1,-t,0],[0,-1,t],[0,1,t],[0,-1,-t],[0,1,-t],[t,0,-1],[t,0,1],[-t,0,-1],[-t,0,1]]
    F=[[0,11,5],[0,5,1],[0,1,7],[0,7,10],[0,10,11],[1,5,9],[5,11,4],[11,10,2],[10,7,6],[7,1,8],[3,9,4],[3,4,2],[3,2,6],[3,6,8],[3,8,9],[4,9,5],[2,4,11],[6,2,10],[8,6,7],[9,8,1]]
    V=[np.array(v,float)/np.linalg.norm(v) for v in V]
    for _ in range(k):
        cache={};nf=[]
        def mid(a,b):
            key=(min(a,b),max(a,b))
            if key not in cache:
                m=V[a]+V[b];V.append(m/np.linalg.norm(m));cache[key]=len(V)-1
            return cache[key]
        for a,b,c in F:
            ab,bc,ca=mid(a,b),mid(b,c),mid(c,a);nf+=[[a,ab,ca],[b,bc,ab],[c,ca,bc],[ab,bc,ca]]
        F=nf
    return np.array(V),np.array(F)
out={}
for name,k,amp in (("bumpy642",3,0.18),("bumpy2562",4,0.22)):
    V,F=icosphere(k)
    th=np.arccos(V[:,2]);ph=np.arctan2(V[:,1],V[:,0])
    r=1+amp*np.sin(6*th)*np.cos(5*ph)
    out[name]={"vertices":(V*r[:,None]*0.05).round(9).tolist(),"faces":F.tolist()}
def fs(path):
    with open(path,"rb") as f:
        f.read(3);f.readline();f.readline()
        nv,nf=np.fromfile(f,">i4",2);V=np.fromfile(f,">f4",nv*3).reshape(nv,3).astype(float);F=np.fromfile(f,">i4",nf*3).reshape(nf,3)
    return V/1000,F
V,F=fs(os.path.join(os.environ["NXR_FSAVERAGE5"],"surf","lh.pial"))
out["fsaverage5_lh_pial"]={"vertices":V.tolist(),"faces":F.tolist(),"source":"$NXR_FSAVERAGE5/surf/lh.pial, mm/1000 (metres), FreeSurfer RAS"}
rng=np.random.default_rng(7)
for k,m in out.items():
    n=len(m["vertices"]);m["pairs"]=rng.integers(0,n,(40 if n<5000 else 60,2)).tolist()
json.dump(out,open("meshes.json","w"))
print({k:len(v["vertices"]) for k,v in out.items()})
