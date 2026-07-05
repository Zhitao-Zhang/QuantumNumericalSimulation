#!/usr/bin/env python3
# =============================================================================
#  diag_dispersion.py  —  诊断"菱形波前"来自 collocated 离散
#
#  两部分：
#   (1) 解析数值频散：collocated 中心差分  vs  交错(FDM)  的相速度极坐标图。
#       证明 collocated 沿轴向/对角线相速度不一致 → 波前变菱形。
#   (2) 从真实快照测量波前各向异性（沿角度提取波前半径），与解析预测对照。
# =============================================================================
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- 物理/网格参数（与 config.h 对齐）----
Vp   = 5000.0
dx   = 10.0
dt   = 1.0e-4
f0   = 30.0
lam  = Vp / f0                 # 主频波长 ~167 m
G    = lam / dx                # 每波长网格数 ~16.7

print(f"[params] Vp={Vp}  dx={dx}  lambda={lam:.1f} m  points/wavelength G={G:.2f}")

# =============================================================================
# (1) 解析数值频散
#
#   连续:      ω = Vp·k                     (完全各向同性)
#
#   交错(FDM)：一阶导用 (f_{i+1/2}-f_{i-1/2})/dx  → 符号 (2/dx)·sin(kx·dx/2)
#       ω_num = (Vp/dx)·sqrt( (2 sin(kx dx/2))² + (2 sin(kz dx/2))² )
#             各向异性小（半角 sin）。
#
#   collocated：一阶导用 (f_{i+1}-f_{i-1})/(2dx) → 符号 (1/dx)·sin(kx·dx)
#       等效间距 2dx；ω_num = (Vp/dx)·sqrt( sin²(kx dx) + sin²(kz dx) )
#       在 kx dx = π（轴向 Nyquist）处 sin=0 → 该方向完全传不动(奇偶解耦)。
#       各向异性大 → 波前菱形。
# =============================================================================
def vph_over_vp(scheme, G_ppw, angles):
    """给定每波长网格数 G，返回各角度的 归一化相速度 v_ph/Vp。"""
    k = 2.0 * np.pi / (G_ppw * dx)      # 目标波数 (对应主频波长)
    out = []
    for phi in angles:
        kx = k * np.cos(phi)
        kz = k * np.sin(phi)
        if scheme == "staggered":
            sx = (2.0 / dx) * np.sin(kx * dx / 2.0)
            sz = (2.0 / dx) * np.sin(kz * dx / 2.0)
        elif scheme == "collocated":
            sx = (1.0 / dx) * np.sin(kx * dx)      # 中心差分：跳中心 → 等效 2dx
            sz = (1.0 / dx) * np.sin(kz * dx)
        else:
            raise ValueError
        omega_num = Vp * np.sqrt(sx * sx + sz * sz)
        omega_exact = Vp * k
        out.append(omega_num / omega_exact)        # = v_ph_num / Vp
    return np.array(out)

angles = np.linspace(0, 2 * np.pi, 361)
# 用一个偏高频分量（比主频短，更能暴露各向异性）：G≈6 网格/波长
for Gtest in [G, 6.0]:
    vs = vph_over_vp("staggered",  Gtest, angles)
    vc = vph_over_vp("collocated", Gtest, angles)
    aniso_s = (vs.max() - vs.min()) / vs.mean() * 100
    aniso_c = (vc.max() - vc.min()) / vc.mean() * 100
    print(f"[G={Gtest:4.1f}] staggered  vph/Vp in [{vs.min():.3f},{vs.max():.3f}] "
          f"anisotropy={aniso_s:5.1f}%")
    print(f"[G={Gtest:4.1f}] collocated vph/Vp in [{vc.min():.3f},{vc.max():.3f}] "
          f"anisotropy={aniso_c:5.1f}%")

# 极坐标对比图（用 G=6 高频分量，各向异性最直观）
Gp = 6.0
vs = vph_over_vp("staggered",  Gp, angles)
vc = vph_over_vp("collocated", Gp, angles)

fig = plt.figure(figsize=(12, 5.5))
ax1 = fig.add_subplot(1, 2, 1, projection="polar")
ax1.plot(angles, vs, "b-", lw=2, label="staggered (FDM)")
ax1.plot(angles, vc, "r-", lw=2, label="collocated (this FEM)")
ax1.plot(angles, np.ones_like(angles), "k--", lw=1, label="exact (isotropic)")
ax1.set_title(f"Phase velocity vs angle  (G={Gp:.0f} pts/wavelength)\n"
              "collocated 沿对角线更快 → 波前菱形", fontsize=11)
ax1.legend(loc="lower right", bbox_to_anchor=(1.15, -0.15))

# 用相速度反推"等时波前形状"：波前半径 ∝ v_ph(φ)
ax2 = fig.add_subplot(1, 2, 2, projection="polar")
ax2.plot(angles, vs, "b-", lw=2, label="staggered → 近圆")
ax2.plot(angles, vc, "r-", lw=2, label="collocated → 菱形")
ax2.fill(angles, vc, "r", alpha=0.15)
ax2.set_title("等时波前形状 (r ∝ v_ph)\n红色明显是 45°方向凸出的菱形", fontsize=11)
ax2.legend(loc="lower right", bbox_to_anchor=(1.15, -0.15))

plt.tight_layout()
plt.savefig("diag_wavefront_anisotropy.png", dpi=120)
print("[saved] diag_wavefront_anisotropy.png")

# =============================================================================
# (2) 从真实快照测量波前各向异性
# =============================================================================
import glob, os

def load_snap(fn):
    return np.loadtxt(fn)   # 200 x 200

# 找一个波前已经传开、但还没到 PML 的时刻
# 选一个波前已明显传开、且未进 PML 的时刻。Vp*t = 5000*it*1e-4 m。
# it=800 → t=0.08s → 400 m → 40 cells（PML 30 cells，中心距边 ~70，安全）。
cand = "snap_txx_it0800.dat"
if os.path.exists(cand):
    field = load_snap(cand)
    nz, nx = field.shape
    print(f"\n[snap] {cand}  shape={field.shape}")
    # 源位置（config: isou=100, jsou=60）。注意文件行=iz, 列=ix
    isou, jsou = 100, 60
    amp = np.abs(field)
    # 屏蔽源附近 8 cells（避免源自身大振幅干扰），从外向内找**最外侧**显著波前
    ang_list = np.linspace(0, 2 * np.pi, 180, endpoint=False)
    r_scan = np.arange(8, 65, 0.5)
    # 每个角度上取该射线剖面的峰值位置作为主波前半径（比阈值更稳）
    radius = []
    for phi in ang_list:
        prof = []
        rr   = []
        for r in r_scan:
            ix = int(round(isou + r * np.cos(phi)))
            iz = int(round(jsou + r * np.sin(phi)))
            if 0 <= ix < nx and 0 <= iz < nz:
                prof.append(amp[iz, ix]); rr.append(r)
        if len(prof) > 5 and max(prof) > 0.05 * amp.max():
            radius.append(rr[int(np.argmax(prof))])   # 峰值处 = 主波前
        else:
            radius.append(np.nan)
    radius = np.array(radius)
    good = ~np.isnan(radius)
    if good.sum() > 10:
        r_axis = np.nanmean(radius[(np.abs(np.cos(2*ang_list)) > 0.9)])  # 轴向 0/90/180/270
        r_diag = np.nanmean(radius[(np.abs(np.cos(2*ang_list)) < 0.2)])  # 对角 45/135...
        print(f"[measure] 波前半径 轴向≈{r_axis:.1f} cells  对角≈{r_diag:.1f} cells  "
              f"diag/axis={r_diag/r_axis:.3f}")
        print("          若 diag/axis 明显 >1 → 对角线传得快 → 菱形（与解析预测一致）")

        figm = plt.figure(figsize=(11, 5))
        axa = figm.add_subplot(1, 2, 1)
        im = axa.imshow(field, cmap="RdBu_r", origin="lower",
                        vmin=-np.percentile(amp, 99.5), vmax=np.percentile(amp, 99.5))
        axa.plot(isou, jsou, "k*", ms=10)
        axa.set_title(f"{cand}  σxx 波场\n(菱形波前肉眼可见)")
        axa.set_xlabel("ix"); axa.set_ylabel("iz")
        plt.colorbar(im, ax=axa, fraction=0.046)

        axb = figm.add_subplot(1, 2, 2, projection="polar")
        axb.plot(ang_list[good], radius[good], "r.-", lw=1)
        axb.set_title("测得波前半径 vs 角度\n(凸向 45° = 菱形)")
        plt.tight_layout()
        plt.savefig("diag_measured_wavefront.png", dpi=120)
        print("[saved] diag_measured_wavefront.png")
    else:
        print("[measure] 波前点太少，可能该时刻波还没传开或已进 PML")
else:
    print(f"\n[snap] {cand} 不存在，跳过实测部分")

print("\n[done]")
