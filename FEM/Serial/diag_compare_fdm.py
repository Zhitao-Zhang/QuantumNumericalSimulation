#!/usr/bin/env python3
# =============================================================================
#  diag_compare_fdm.py  —  FEM vs FDM 波前各向同性对照（决定性诊断）
#
#  用压强 p = (σxx+σzz)/2 作对照量：
#    - P 波压强场是旋转对称的（无 σxx 的 cos²θ 辐射方向图污染），
#      所以波前形状 = 真正的数值频散各向异性。
#    - 同时刻 it=800，FEM 与 FDM 用完全相同的参数/源/介质。
#  测量：沿 180 个角度提取波前前沿半径，比较 轴向 vs 对角 半径比。
# =============================================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FEM_DIR = "."
FDM_DIR = "/home/wanwb/zzt/Quantum/FDM/ImplicitV2"
IT = "0800"
isou, jsou = 100, 60   # 列=ix, 行=iz

def pressure(dir_):
    sxx = np.loadtxt(f"{dir_}/snap_txx_it{IT}.dat")
    szz = np.loadtxt(f"{dir_}/snap_tzz_it{IT}.dat")
    return 0.5 * (sxx + szz)

def front_radius(field, isou, jsou, n_ang=180, rmin=8, rmax=64):
    """沿每个角度提取主波前半径（剖面峰值位置）。"""
    nz, nx = field.shape
    amp = np.abs(field)
    ang = np.linspace(0, 2*np.pi, n_ang, endpoint=False)
    r_scan = np.arange(rmin, rmax, 0.5)
    rad = []
    for phi in ang:
        prof, rr = [], []
        for r in r_scan:
            ix = int(round(isou + r*np.cos(phi)))
            iz = int(round(jsou + r*np.sin(phi)))
            if 0 <= ix < nx and 0 <= iz < nz:
                prof.append(amp[iz, ix]); rr.append(r)
        if len(prof) > 5 and max(prof) > 0.08*amp.max():
            rad.append(rr[int(np.argmax(prof))])
        else:
            rad.append(np.nan)
    return ang, np.array(rad)

def aniso_report(name, ang, rad):
    axis_mask = np.abs(np.cos(2*ang)) > 0.9   # 0/90/180/270
    diag_mask = np.abs(np.cos(2*ang)) < 0.15  # 45/135/225/315
    r_ax = np.nanmean(rad[axis_mask])
    r_dg = np.nanmean(rad[diag_mask])
    r_all = rad[~np.isnan(rad)]
    aniso = (np.nanmax(rad) - np.nanmin(rad)) / np.nanmean(rad) * 100
    print(f"[{name}] 波前半径 轴向={r_ax:5.1f}  对角={r_dg:5.1f}  "
          f"diag/axis={r_dg/r_ax:5.3f}  总各向异性={aniso:4.1f}%  "
          f"(r范围 {r_all.min():.0f}-{r_all.max():.0f})")
    return r_ax, r_dg, aniso

p_fem = pressure(FEM_DIR)
p_fdm = pressure(FDM_DIR)
print(f"[load] FEM p shape={p_fem.shape}  FDM p shape={p_fdm.shape}")
print(f"[expect] Vp*t = 5000*{int(IT)}*1e-4 = {5000*int(IT)*1e-4:.0f} m = "
      f"{5000*int(IT)*1e-4/10:.0f} cells\n")

ang_f, rad_f = front_radius(p_fem, isou, jsou)
ang_d, rad_d = front_radius(p_fdm, isou, jsou)
aniso_report("FEM ", ang_f, rad_f)
aniso_report("FDM ", ang_d, rad_d)

# ---- 图：并排波场 + 波前极坐标 ----
fig = plt.figure(figsize=(14, 5))
vlim_f = np.percentile(np.abs(p_fem), 99.5)
vlim_d = np.percentile(np.abs(p_fdm), 99.5)

ax1 = fig.add_subplot(1, 3, 1)
ax1.imshow(p_fem, cmap="RdBu_r", origin="lower", vmin=-vlim_f, vmax=vlim_f)
ax1.plot(isou, jsou, "k*", ms=10); ax1.set_title(f"FEM pressure it={IT}")
ax1.set_xlabel("ix"); ax1.set_ylabel("iz")

ax2 = fig.add_subplot(1, 3, 2)
ax2.imshow(p_fdm, cmap="RdBu_r", origin="lower", vmin=-vlim_d, vmax=vlim_d)
ax2.plot(isou, jsou, "k*", ms=10); ax2.set_title(f"FDM pressure it={IT}")
ax2.set_xlabel("ix"); ax2.set_ylabel("iz")

ax3 = fig.add_subplot(1, 3, 3, projection="polar")
gf = ~np.isnan(rad_f); gd = ~np.isnan(rad_d)
ax3.plot(ang_f[gf], rad_f[gf], "r.-", lw=1, ms=3, label="FEM front")
ax3.plot(ang_d[gd], rad_d[gd], "b.-", lw=1, ms=3, label="FDM front")
ax3.set_title("Wavefront radius vs angle\n(circle=isotropic, lobed=anisotropic)")
ax3.legend(loc="lower right", bbox_to_anchor=(1.25, -0.12))

plt.tight_layout()
plt.savefig("diag_fem_vs_fdm_wavefront.png", dpi=120)
print("\n[saved] diag_fem_vs_fdm_wavefront.png")
