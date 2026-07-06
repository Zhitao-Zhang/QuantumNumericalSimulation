#!/usr/bin/env python3
# 验证方案B(位移二阶)波前是否为圆：对比 新FEM / 旧FEM备份 / FDM
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

isou, jsou = 100, 60
IT = "0800"
FDM = "/home/wanwb/zzt/Quantum/FDM/ImplicitV2"

def pressure(d):
    return 0.5*(np.loadtxt(f"{d}/snap_txx_it{IT}.dat") + np.loadtxt(f"{d}/snap_tzz_it{IT}.dat"))

def front(field, n=180, rmin=8, rmax=64):
    nz,nx = field.shape; amp=np.abs(field)
    ang=np.linspace(0,2*np.pi,n,endpoint=False); rs=np.arange(rmin,rmax,0.5); rad=[]
    for phi in ang:
        prof,rr=[],[]
        for r in rs:
            ix=int(round(isou+r*np.cos(phi))); iz=int(round(jsou+r*np.sin(phi)))
            if 0<=ix<nx and 0<=iz<nz: prof.append(amp[iz,ix]); rr.append(r)
        rad.append(rr[int(np.argmax(prof))] if (len(prof)>5 and max(prof)>0.08*amp.max()) else np.nan)
    return ang,np.array(rad)

def rep(name,ang,rad):
    ax=np.abs(np.cos(2*ang))>0.9; dg=np.abs(np.cos(2*ang))<0.15
    ra=np.nanmean(rad[ax]); rd=np.nanmean(rad[dg])
    an=(np.nanmax(rad)-np.nanmin(rad))/np.nanmean(rad)*100
    print(f"[{name:16s}] axis={ra:5.1f} diag={rd:5.1f} diag/axis={rd/ra:5.3f} aniso={an:5.1f}%")
    return ra,rd

# 新 FEM 就在当前目录
p_new = pressure(".")
p_fdm = pressure(FDM)
print(f"expect Vp*t=40 cells")
an,rn = front(p_new); rep("FEM-B (位移)",an,rn)
af,rf = front(p_fdm); rep("FDM (交错)",af,rf)

fig=plt.figure(figsize=(14,5))
for i,(p,t) in enumerate([(p_new,"FEM-B displacement (this)"),(p_fdm,"FDM staggered")]):
    ax=fig.add_subplot(1,3,i+1); vl=np.percentile(np.abs(p),99.5)
    ax.imshow(p,cmap="RdBu_r",origin="lower",vmin=-vl,vmax=vl)
    ax.plot(isou,jsou,"k*",ms=10); ax.set_title(f"{t}  pressure it={IT}")
    ax.set_xlabel("ix"); ax.set_ylabel("iz")
ax=fig.add_subplot(1,3,3,projection="polar")
gn=~np.isnan(rn); gf=~np.isnan(rf)
ax.plot(an[gn],rn[gn],"g.-",lw=1,ms=3,label="FEM-B (this)")
ax.plot(af[gf],rf[gf],"b.-",lw=1,ms=3,label="FDM")
ax.set_title("wavefront radius vs angle\n(circle=isotropic)")
ax.legend(loc="lower right",bbox_to_anchor=(1.3,-0.1))
plt.tight_layout(); plt.savefig("diag_planB_vs_fdm.png",dpi=120)
print("[saved] diag_planB_vs_fdm.png")
