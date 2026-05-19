/* =============================================================================
 *  2D Velocity–Stress Elastic Wave Equation
 *  Crank–Nicolson Implicit Time Scheme + ADE-PML  (Serial CPU C++)
 *
 *  在 v1（无 PML）基础上加入 PML，按 Quantum/Implicit/main.cu 的 PML 设置：
 *    四周 NP 个网格作为 PML 吸收区，d(α) 在 PML 内呈二次衰减分布，
 *    d0 = 3 V0 ln(1/RC) / (2 DP)，DP = NP·h。
 *
 *  关键设计（与 Implicit/main.cu 显式 PML 的不同）：
 *      Implicit 版本里 ψ_α 是显式更新 + Picard 迭代逼近 CN；
 *      本版本把 ψ_α 真正当作未知量，与 (vx, vz, σxx, σzz, σxz) 一起进入
 *      Ax=b 线性系统，每步一次直接求解（无 Picard 内迭代）。
 *
 *  ===========================================================================
 *  PML 的 ADE (Auxiliary Differential Equation) 形式
 *  ===========================================================================
 *  把每个空间导数替换为「拉伸坐标」下的导数：
 *
 *      ∂̃_α  =  ∂_α  +  ψ_α                                  (ADE 形式)
 *
 *  其中 ψ_α 满足一阶 ODE
 *
 *      ∂ψ_α/∂t  +  d_α ψ_α  =  -d_α (∂_α field)             (*)
 *
 *  使用指数积分 + 梯形规则（与 Implicit/main.cu 完全一致）对 (*) 时间离散：
 *
 *      ψ_α^{n+1} = e^{-d_α Δt} ψ_α^n
 *                + (-d_α Δt / 2)·[ e^{-d_α Δt} (∂_α field)^n + (∂_α field)^{n+1} ]
 *
 *  整理成「左端=新时层未知量、右端=旧时层已知量」的 CN 矩阵行形式：
 *
 *      ψ_α^{n+1}  +  (d_α Δt/2)·(∂_α field)^{n+1}
 *           =  E_α · [ ψ_α^n  -  (d_α Δt/2)·(∂_α field)^n ],
 *      其中 E_α = e^{-d_α Δt}.
 *
 *  在 d_α = 0 区（PML 外）：E_α = 1，方程退化为  ψ^{n+1} = ψ^n，
 *  若初值 ψ ≡ 0 则 ψ 永远是 0，矩阵自动「关闭」PML 行为，结果与 v1 完全一致。
 *
 *  ===========================================================================
 *  CN 主方程的修改（把 ∂α 替换为 ∂α + ψ_α）
 *  ===========================================================================
 *  令 C=λ+2μ, L=λ, M=μ, cax=Δt/(2 dx), caz=Δt/(2 dz). 五个主方程变为：
 *
 *      [vx eq]   ρx vx^{n+1}
 *                - (Δt/2)(∂xσxx + ψ_xSxx)^{n+1}
 *                - (Δt/2)(∂zσxz + ψ_zSxz)^{n+1}
 *              = ρx vx^{n} + (Δt/2)[same terms]^{n}
 *
 *      [vz eq]   ρz vz^{n+1}
 *                - (Δt/2)(∂xσxz + ψ_xSxz)^{n+1}
 *                - (Δt/2)(∂zσzz + ψ_zSzz)^{n+1}
 *              = ρz vz^{n} + (Δt/2)[same terms]^{n}
 *
 *      [σxx eq]  σxx^{n+1}
 *                - (Δt/2) C (∂xvx + ψ_xVx)^{n+1}
 *                - (Δt/2) L (∂zvz + ψ_zVz)^{n+1}
 *              = σxx^{n} + (Δt/2)[same terms]^{n}
 *
 *      [σzz eq]  σzz^{n+1}
 *                - (Δt/2) L (∂xvx + ψ_xVx)^{n+1}
 *                - (Δt/2) C (∂zvz + ψ_zVz)^{n+1}
 *              = σzz^{n} + (Δt/2)[same terms]^{n}
 *
 *      [σxz eq]  σxz^{n+1}
 *                - (Δt/2) M (∂xvz + ψ_xVz)^{n+1}
 *                - (Δt/2) M (∂zvx + ψ_zVx)^{n+1}
 *              = σxz^{n} + (Δt/2)[same terms]^{n}
 *
 *  ===========================================================================
 *  8 个 ψ 辅助场（按物理位置分组）
 *  ===========================================================================
 *      位置在 vx     (Nz × (Nx+1)) ：ψ_xSxx [对 ∂xσxx], ψ_zSxz [对 ∂zσxz]
 *      位置在 vz     ((Nz+1) × Nx) ：ψ_xSxz [对 ∂xσxz], ψ_zSzz [对 ∂zσzz]
 *      位置在 σxx/σzz (Nz × Nx)    ：ψ_xVx  [对 ∂xvx ], ψ_zVz  [对 ∂zvz ]  (σxx 和 σzz 共用)
 *      位置在 σxz   ((Nz+1)×(Nx+1)):ψ_xVz  [对 ∂xvz ], ψ_zVx  [对 ∂zvx ]
 *
 *  ψ_α 的 d_α 取决于：①导数方向(x 或 z)、② ψ 所在的物理位置(整点 / 半点)。
 *  本程序预先在 5 条 1D 阻尼曲线上采样：
 *      damp_x_int [ix]      x = ix·dx,         ix=0..Nx-1
 *      damp_x_half[ix]      x = (ix+1/2)·dx,   ix=0..Nx
 *      damp_z_int [iz]      z = iz·dz,         iz=0..Nz-1
 *      damp_z_hv  [iz]      z = (iz+1/2)·dz,   iz=0..Nz   (vz 与 ψ_zSzz)
 *      damp_z_hxz [iz]      z = (iz-1/2)·dz,   iz=0..Nz   (σxz 与 ψ_zVx)
 *
 *  对应 8 个 ψ：
 *      ψ_xSxx :  damp_x_half[ix]           ψ_zSxz :  damp_z_int [iz]
 *      ψ_xSxz :  damp_x_int [ix]           ψ_zSzz :  damp_z_hv  [iz]
 *      ψ_xVx  :  damp_x_int [ix]           ψ_zVz  :  damp_z_int [iz]
 *      ψ_xVz  :  damp_x_half[ix]           ψ_zVx  :  damp_z_hxz [iz]
 *
 *  ===========================================================================
 *  未知向量 x 的 13 块 (按行块顺序 = 块偏移顺序)
 *  ===========================================================================
 *      0 : vx     [Nz × (Nx+1)]
 *      1 : vz     [(Nz+1) × Nx]
 *      2 : σxx    [Nz × Nx]
 *      3 : σzz    [Nz × Nx]
 *      4 : σxz    [(Nz+1) × (Nx+1)]
 *      5 : ψ_xSxx [Nz × (Nx+1)]    same shape as vx
 *      6 : ψ_zSxz [Nz × (Nx+1)]    same shape as vx
 *      7 : ψ_xSxz [(Nz+1) × Nx]    same shape as vz
 *      8 : ψ_zSzz [(Nz+1) × Nx]    same shape as vz
 *      9 : ψ_xVx  [Nz × Nx]        same shape as σxx
 *     10 : ψ_zVz  [Nz × Nx]        same shape as σxx
 *     11 : ψ_xVz  [(Nz+1)×(Nx+1)]  same shape as σxz
 *     12 : ψ_zVx  [(Nz+1)×(Nx+1)]  same shape as σxz
 *
 *  ===========================================================================
 *  边界规则（与 v1 完全一致：缺邻居 -> 跳过该差分项）
 *  ===========================================================================
 *  每条 ψ 行使用的 ∂α 必须与对应主方程里那条 ∂α 一字不差，否则方程不闭合。
 *  例如 vz 方程里 ∂xσxz 的合法条件是 ix≥1 且 iz≤Nz-1，那么 ψ_xSxz 行里
 *  的 σxz 系数也必须遵守同样条件（不合法时只剩 ψ^{n+1} = E·ψ^n，没有 σ 耦合）。
 *
 *  ===========================================================================
 *  完整 PML d(α) 阻尼曲线（与 Implicit/main.cu 同 form）
 *  ===========================================================================
 *      d0   = 3 V0 ln(1/RC) / (2 DP)        ; V0 = Vp_val, RC = 1e-6, DP = NP·h
 *      x < xleft     :  d(x) = d0 · ((xleft - x)/DP)^2
 *      x > xright    :  d(x) = d0 · ((x - xright)/DP)^2
 *      otherwise     :  d(x) = 0
 *      xleft  =  NP   · dx
 *      xright = (Nx-1-NP)·dx
 *      z 方向同构。
 *
 *  作者：CN-SGFD+PML 串行原型
 * ===========================================================================*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <map>
#include <algorithm>

// ----------------------------------------------------------------------------
// 1. 网格与时间参数（可在编译时通过 -DNZ=.. 等覆盖）
// ----------------------------------------------------------------------------
#ifndef NZ
#define NZ 200         // 应力网格在 z 方向的大小 Nz
#endif
#ifndef NX
#define NX 200         // 应力网格在 x 方向的大小 Nx
#endif
#ifndef NT
#define NT 4500         // 时间步数
#endif
#ifndef NPMLDEF
#define NPMLDEF 30     // 四周 PML 厚度（单位：网格数）
#endif
#ifndef DX_VAL
#define DX_VAL 10.0    // x 方向网格步长 (m)
#endif
#ifndef DZ_VAL
#define DZ_VAL 10.0    // z 方向网格步长 (m)
#endif
#ifndef DT_VAL
#define DT_VAL 1.0e-4  // 时间步长 (s)，CN 无条件稳定
#endif


static const int    Nz = NZ;
static const int    Nx = NX;
static const int    NP = NPMLDEF;  // PML thickness in cells
static const double dx = DX_VAL;   // x 方向网格步长 (m)
static const double dz = DZ_VAL;   // z 方向网格步长 (m)
static const double dt = DT_VAL;   // 时间步长 (s)，CN 无条件稳定

// ----------------------------------------------------------------------------
// 2. 介质参数：固体背景 + 含水井孔（非均质模型）
// ----------------------------------------------------------------------------
//   井孔几何：贯穿 z 方向、宽度 = 2·NBH 个格点，中心在源 ix=isou
//             即井孔占据 ix ∈ [isou-NBH, isou+NBH-1] 共 2·NBH 格
//   井内：水（Vp=1500, Vs=0, ρ=1000）
//   井外：固体（Vp=5000, Vs=3300, ρ=2450 — 可通过宏覆盖）
//   NBH=0 表示退化为均匀介质（与 v2 行为相同）
#ifndef NBH
#define NBH 40        // 井孔半宽（cells），整宽 = 2·NBH
#endif

// ---- 固体（井外）----
static const double Vp_val   = 5000.0;
static const double Vs_val   = 3300.0;
static const double rho_val  = 2450.0;
static const double mu_val   = rho_val * Vs_val * Vs_val;                // μ
static const double lam2u_val= rho_val * Vp_val * Vp_val;                // λ + 2μ
static const double lam_val  = lam2u_val - 2.0 * mu_val;                 // λ


// ---- 固体2（井外）----
static const double Vp_val2   = 6000.0;
static const double Vs_val2   = 3400.0;
static const double rho_val2  = 2770.0;
static const double mu_val2   = rho_val2 * Vs_val2 * Vs_val2;                // μ
static const double lam2u_val2= rho_val2 * Vp_val2 * Vp_val2;                // λ + 2μ
static const double lam_val2  = lam2u_val2 - 2.0 * mu_val2;                 // λ

// ---- 流体（井内）----
static const double Vp_fluid   = 1500.0;
static const double Vs_fluid   = 0.0;
static const double rho_fluid  = 1000.0;
static const double mu_fluid   = rho_fluid * Vs_fluid * Vs_fluid;        // = 0
static const double lam2u_fluid= rho_fluid * Vp_fluid * Vp_fluid;        // = ρ·Vp²
static const double lam_fluid  = lam2u_fluid - 2.0 * mu_fluid;           // = λ (= ρ·Vp²)

// ----------------------------------------------------------------------------
// 3. PML 阻尼曲线参数（与 Quantum/Implicit/main.cu 一致）
// ----------------------------------------------------------------------------
static const double RC_pml   = 1.0e-6;            // 反射系数
static const double V0_pml   = Vp_val;            // 参考速度
static const double DPx      = (double)NP * dx;   // PML 厚度（物理单位，x）
static const double DPz      = (double)NP * dz;   // PML 厚度（物理单位，z）
static const double d0_x     = (NP > 0)
                             ? 3.0 * V0_pml * std::log(1.0 / RC_pml) / (2.0 * DPx)
                             : 0.0;
static const double d0_z     = (NP > 0)
                             ? 3.0 * V0_pml * std::log(1.0 / RC_pml) / (2.0 * DPz)
                             : 0.0;
static const double xleft_p  = (double)NP * dx;
static const double xright_p = (double)(NX - 1 - NP) * dx;
static const double zleft_p  = (double)NP * dz;
static const double zright_p = (double)(NZ - 1 - NP) * dz;

static double dampx_at(double x)
{
    if (NP <= 0) return 0.0;
    if (x < xleft_p)  return d0_x * std::pow((xleft_p  - x) / DPx, 2);
    if (x > xright_p) return d0_x * std::pow((x - xright_p) / DPx, 2);
    return 0.0;
}
static double dampz_at(double z)
{
    if (NP <= 0) return 0.0;
    if (z < zleft_p)  return d0_z * std::pow((zleft_p  - z) / DPz, 2);
    if (z > zright_p) return d0_z * std::pow((z - zright_p) / DPz, 2);
    return 0.0;
}

// 5 条采样曲线 + 各自的 E = exp(-d·dt)
static std::vector<double> g_dx_int,  g_ex_int;     // damp & exp at x = ix·dx,        ix=0..Nx-1
static std::vector<double> g_dx_half, g_ex_half;    // damp & exp at x = (ix+0.5)·dx,  ix=0..Nx
static std::vector<double> g_dz_int,  g_ez_int;     // damp & exp at z = iz·dz,        iz=0..Nz-1
static std::vector<double> g_dz_hv,   g_ez_hv;      // damp & exp at z = (iz+0.5)·dz,  iz=0..Nz
static std::vector<double> g_dz_hxz,  g_ez_hxz;     // damp & exp at z = (iz-0.5)·dz,  iz=0..Nz

static void build_damping_profiles()
{
    g_dx_int .resize(Nx);     g_ex_int .resize(Nx);
    g_dx_half.resize(Nx + 1); g_ex_half.resize(Nx + 1);
    g_dz_int .resize(Nz);     g_ez_int .resize(Nz);
    g_dz_hv  .resize(Nz + 1); g_ez_hv  .resize(Nz + 1);
    g_dz_hxz .resize(Nz + 1); g_ez_hxz .resize(Nz + 1);

    for (int ix = 0; ix < Nx; ix++) {
        g_dx_int[ix] = dampx_at(ix * dx);
        g_ex_int[ix] = std::exp(-g_dx_int[ix] * dt);
    }
    for (int ix = 0; ix <= Nx; ix++) {
        g_dx_half[ix] = dampx_at((ix + 0.5) * dx);
        g_ex_half[ix] = std::exp(-g_dx_half[ix] * dt);
    }
    for (int iz = 0; iz < Nz; iz++) {
        g_dz_int[iz] = dampz_at(iz * dz);
        g_ez_int[iz] = std::exp(-g_dz_int[iz] * dt);
    }
    for (int iz = 0; iz <= Nz; iz++) {
        g_dz_hv [iz] = dampz_at((iz + 0.5) * dz);
        g_ez_hv [iz] = std::exp(-g_dz_hv [iz] * dt);
        g_dz_hxz[iz] = dampz_at((iz - 0.5) * dz);
        g_ez_hxz[iz] = std::exp(-g_dz_hxz[iz] * dt);
    }
}

// ----------------------------------------------------------------------------
// 4. 声源参数（Ricker 子波，应力点源，加到 σxx 和 σzz 的 RHS）
// ----------------------------------------------------------------------------
static const double PI = 3.14159265358979323846;
#ifndef F0_HZ
#define F0_HZ 30.0     // Ricker 子波主频 (Hz)
#endif
static const double F0 = F0_HZ;
static const double T0 = 1.2 / F0;
static const int    jsou = NZ / 2 - 40;
static const int    isou = NX / 2;
static int          jrec = NZ / 2 + 10;
static int          irec = NX / 2;

// ----------------------------------------------------------------------------
// 5. 未知向量布局：5 块主场 + 8 块 ψ 辅助场
// ----------------------------------------------------------------------------
// 主块尺寸
static int Nvx, Nvz, Nxx, Nzz, Nxz;
// 各自的块偏移
static int off_vx, off_vz, off_xx, off_zz, off_xz;
// PML ψ 块的偏移
static int off_pxSxx, off_pzSxz;   // 位置=vx
static int off_pxSxz, off_pzSzz;   // 位置=vz
static int off_pxVx,  off_pzVz;    // 位置=σxx/σzz
static int off_pxVz,  off_pzVx;    // 位置=σxz
static int Ntot;

// 主块的下标
static inline int id_vx(int iz, int ix) { return off_vx + ix + iz * (Nx + 1); }
static inline int id_vz(int iz, int ix) { return off_vz + ix + iz * Nx;       }
static inline int id_xx(int iz, int ix) { return off_xx + ix + iz * Nx;       }
static inline int id_zz(int iz, int ix) { return off_zz + ix + iz * Nx;       }
static inline int id_xz(int iz, int ix) { return off_xz + ix + iz * (Nx + 1); }

// 8 个 ψ 的下标
static inline int id_pxSxx(int iz, int ix) { return off_pxSxx + ix + iz * (Nx + 1); }
static inline int id_pzSxz(int iz, int ix) { return off_pzSxz + ix + iz * (Nx + 1); }
static inline int id_pxSxz(int iz, int ix) { return off_pxSxz + ix + iz * Nx;       }
static inline int id_pzSzz(int iz, int ix) { return off_pzSzz + ix + iz * Nx;       }
static inline int id_pxVx (int iz, int ix) { return off_pxVx  + ix + iz * Nx;       }
static inline int id_pzVz (int iz, int ix) { return off_pzVz  + ix + iz * Nx;       }
static inline int id_pxVz (int iz, int ix) { return off_pxVz  + ix + iz * (Nx + 1); }
static inline int id_pzVx (int iz, int ix) { return off_pzVx  + ix + iz * (Nx + 1); }

// 单一波场的本地下标（不带 block offset），用于裸数组寻址
static inline int lx_vx(int iz, int ix) { return ix + iz * (Nx + 1); }
static inline int lx_vz(int iz, int ix) { return ix + iz * Nx;       }
static inline int lx_xx(int iz, int ix) { return ix + iz * Nx;       }
static inline int lx_zz(int iz, int ix) { return ix + iz * Nx;       }
static inline int lx_xz(int iz, int ix) { return ix + iz * (Nx + 1); }

// ----------------------------------------------------------------------------
// 6. 逐格点介质参数（非均匀介质：含水井孔 + 固体背景）
// ----------------------------------------------------------------------------
//   定义在整点格的「单元」属性（Nz × Nx）：rho_cell, lam_cell, mu_cell；
//   再插值到 5 个交错位置上：
//      g_rho_x  : 在 vx 位置 (x=ix+0.5, z=iz)，大小 Nz × (Nx+1)，沿 x 算术平均
//      g_rho_z  : 在 vz 位置 (x=ix,     z=iz+0.5)，大小 (Nz+1) × Nx，沿 z 算术平均
//      g_C      : 在 σxx/σzz 位置 (x=ix, z=iz)，大小 Nz × Nx，= λ+2μ 直接取
//      g_L      : 同上位置，= λ 直接取
//      g_M      : 在 σxz 位置 (x=ix+0.5, z=iz-0.5)，大小 (Nz+1) × (Nx+1)，4 邻平均
//   边界外的"虚邻居"采用最邻近 clamp（不影响物理结果，因 README §6 的差分跳过策略）。
//   井外格点取固体 (Vp_val, Vs_val, rho_val)；井内取流体 (Vp_fluid, Vs_fluid, rho_fluid)。

static std::vector<double> g_rho_x, g_rho_z;
static std::vector<double> g_C, g_L, g_M;

static void build_medium()
{
    // (1) 在整点 (iz, ix) 上设置 cell-wise 属性
    std::vector<double> rho_cell((size_t)Nz * Nx);
    std::vector<double> lam_cell((size_t)Nz * Nx);
    std::vector<double> mu_cell ((size_t)Nz * Nx);

    int ix_lo = 100;
    int ix_hi = 120;

    int n_fluid_cells = 0;
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            bool in_borehole = (NBH > 0) && (iz >= ix_lo) && (iz < ix_hi);
            double r, vp, vs;
            if (in_borehole) {
                r = rho_val2; vp = Vp_val2; vs = Vs_val2;
                n_fluid_cells++;
            } else {
                r = rho_val;   vp = Vp_val;   vs = Vs_val;
            } 
            const double mu_c    = r * vs * vs;
            const double lam2u_c = r * vp * vp;
            const double lam_c   = lam2u_c - 2.0 * mu_c;
            rho_cell[lx_xx(iz, ix)] = r;
            lam_cell[lx_xx(iz, ix)] = lam_c;
            mu_cell [lx_xx(iz, ix)] = mu_c;
        }
    }

    auto clampi = [](int v, int lo, int hi){ return v < lo ? lo : (v > hi ? hi : v); };

    // (2) g_rho_x: vx[iz, ix] 位置 (x=ix+0.5)，介于 cell ix 和 cell ix+1 之间
    g_rho_x.assign((size_t)Nz * (Nx + 1), 0.0);
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int ixL = clampi(ix,     0, Nx - 1);
            int ixR = clampi(ix + 1, 0, Nx - 1);
            g_rho_x[lx_vx(iz, ix)] = 0.5 * (rho_cell[lx_xx(iz, ixL)] + rho_cell[lx_xx(iz, ixR)]);
        }
    }

    // (3) g_rho_z: vz[iz, ix] 位置 (z=iz+0.5)，介于 cell iz 和 cell iz+1 之间
    g_rho_z.assign((size_t)(Nz + 1) * Nx, 0.0);
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int izL = clampi(iz,     0, Nz - 1);
            int izR = clampi(iz + 1, 0, Nz - 1);
            g_rho_z[lx_vz(iz, ix)] = 0.5 * (rho_cell[lx_xx(izL, ix)] + rho_cell[lx_xx(izR, ix)]);
        }
    }

    // (4) g_C, g_L: σxx/σzz 在 (ix, iz) 整点上，直接取 cell 值
    g_C.assign((size_t)Nz * Nx, 0.0);
    g_L.assign((size_t)Nz * Nx, 0.0);
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double mu_c  = mu_cell [lx_xx(iz, ix)];
            const double lam_c = lam_cell[lx_xx(iz, ix)];
            g_C[lx_xx(iz, ix)] = lam_c + 2.0 * mu_c;     // λ+2μ
            g_L[lx_xx(iz, ix)] = lam_c;                  // λ
        }
    }

    // (5) g_M: σxz 在 (ix+0.5, iz-0.5)，4 个相邻整点 cell 算术平均
    //   注：流-固界面用算术平均，文献中有时也用调和平均；这里选算术保持 μ 在
    //   界面附近平滑过渡（不会被流体的 μ=0 完全压制），数值更稳健。
    g_M.assign((size_t)(Nz + 1) * (Nx + 1), 0.0);
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int izs[2] = { clampi(iz - 1, 0, Nz - 1), clampi(iz, 0, Nz - 1) };
            int ixs[2] = { clampi(ix,     0, Nx - 1), clampi(ix + 1, 0, Nx - 1) };
            double s = 0.0;
            for (int a = 0; a < 2; a++)
                for (int b = 0; b < 2; b++)
                    s += mu_cell[lx_xx(izs[a], ixs[b])];
            g_M[lx_xz(iz, ix)] = 0.25 * s;
        }
    }

    // 打印介质摘要
    printf("Medium: solid (Vp=%g, Vs=%g, rho=%g) + fluid borehole (Vp=%g, Vs=%g, rho=%g)\n",
           Vp_val, Vs_val, rho_val, Vp_fluid, Vs_fluid, rho_fluid);
    if (NBH > 0) {
        printf("  Borehole: vertical column, ix ∈ [%d, %d) (%d cells wide), centered at isou=%d\n",
               ix_lo, ix_hi, 2 * NBH, isou);
        printf("  Fluid cells: %d / %d (%.2f%%)\n",
               n_fluid_cells, Nz * Nx, 100.0 * n_fluid_cells / (double)(Nz * Nx));
    } else {
        printf("  NBH=0  -> homogeneous solid medium\n");
    }
}

// ----------------------------------------------------------------------------
// 7. CSR 矩阵 + Jacobi 预处理 BiCGSTAB（与 v1 相同）
// ----------------------------------------------------------------------------
struct CSRMatrix {
    int n;
    std::vector<int>    rowptr;
    std::vector<int>    colind;
    std::vector<double> values;
    std::vector<double> diag;

    void build_from_rows(int n_, std::vector<std::map<int,double>> &rows)
    {
        n = n_;
        rowptr.assign(n + 1, 0);
        for (int i = 0; i < n; i++)
            rowptr[i + 1] = rowptr[i] + (int)rows[i].size();
        int nnz = rowptr[n];
        colind.resize(nnz);
        values.resize(nnz);
        diag.assign(n, 0.0);
        int k = 0;
        for (int i = 0; i < n; i++) {
            for (auto &kv : rows[i]) {
                colind[k] = kv.first;
                values[k] = kv.second;
                if (kv.first == i) diag[i] = kv.second;
                k++;
            }
        }
        for (int i = 0; i < n; i++)
            if (diag[i] == 0.0) diag[i] = 1.0;
    }

    void multiply(const std::vector<double> &x, std::vector<double> &y) const
    {
        if ((int)y.size() != n) y.assign(n, 0.0);
        for (int i = 0; i < n; i++) {
            double s = 0.0;
            for (int k = rowptr[i]; k < rowptr[i + 1]; k++)
                s += values[k] * x[colind[k]];
            y[i] = s;
        }
    }
};

static int bicgstab(const CSRMatrix &A,
                    const std::vector<double> &b,
                    std::vector<double> &x,
                    int maxiter, double tol, double *final_res)
{
    const int n = A.n;
    std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0),
                        s(n), t(n), ph(n), sh(n);

    A.multiply(x, r);
    for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
    r0 = r;

    double bnorm = 0.0;
    for (int i = 0; i < n; i++) bnorm += b[i] * b[i];
    bnorm = std::sqrt(bnorm);
    if (bnorm == 0.0) bnorm = 1.0;

    double rho_prev = 1.0, alpha = 1.0, omega = 1.0;

    for (int iter = 0; iter < maxiter; iter++) {
        double rho = 0.0;
        for (int i = 0; i < n; i++) rho += r0[i] * r[i];
        if (std::fabs(rho) < 1e-300) {
            r0 = r;
            rho = 0.0;
            for (int i = 0; i < n; i++) rho += r0[i] * r[i];
            rho_prev = 1.0; alpha = 1.0; omega = 1.0;
            std::fill(p.begin(), p.end(), 0.0);
            std::fill(v.begin(), v.end(), 0.0);
        }
        double beta = (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

        for (int i = 0; i < n; i++) ph[i] = p[i] / A.diag[i];
        A.multiply(ph, v);

        double r0v = 0.0;
        for (int i = 0; i < n; i++) r0v += r0[i] * v[i];
        if (std::fabs(r0v) < 1e-300) { if (final_res) *final_res = 1.0; return iter + 1; }
        alpha = rho / r0v;

        for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
        double snorm = 0.0;
        for (int i = 0; i < n; i++) snorm += s[i] * s[i];
        snorm = std::sqrt(snorm);
        if (snorm / bnorm < tol) {
            for (int i = 0; i < n; i++) x[i] += alpha * ph[i];
            if (final_res) *final_res = snorm / bnorm;
            return iter + 1;
        }

        for (int i = 0; i < n; i++) sh[i] = s[i] / A.diag[i];
        A.multiply(sh, t);
        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < n; i++) { ts += t[i] * s[i]; tt += t[i] * t[i]; }
        if (tt < 1e-300) { if (final_res) *final_res = 1.0; return iter + 1; }
        omega = ts / tt;

        for (int i = 0; i < n; i++)
            x[i] += alpha * ph[i] + omega * sh[i];
        for (int i = 0; i < n; i++)
            r[i] = s[i] - omega * t[i];

        double rnorm = 0.0;
        for (int i = 0; i < n; i++) rnorm += r[i] * r[i];
        rnorm = std::sqrt(rnorm);
        if (rnorm / bnorm < tol) {
            if (final_res) *final_res = rnorm / bnorm;
            return iter + 1;
        }
        if (std::fabs(omega) < 1e-300) {
            if (final_res) *final_res = rnorm / bnorm;
            return iter + 1;
        }
        rho_prev = rho;
    }
    {
        std::vector<double> tmp(n);
        A.multiply(x, tmp);
        double rn = 0.0;
        for (int i = 0; i < n; i++) { double d = tmp[i] - b[i]; rn += d * d; }
        if (final_res) *final_res = std::sqrt(rn) / bnorm;
    }
    return maxiter;
}

// ----------------------------------------------------------------------------
// 7. 矩阵 A 组装：13 行块
//
//    主方程 5 块加入 -Δt/2 (或 -Δt·C/2 等) 对 ψ 的耦合；
//    ψ 8 块每行 = 1·ψ + (d·Δt/2)·(∂α applied to σ or v).
//
//    边界规则：和对应主方程的 ∂α 跳过条件保持一致。
// ----------------------------------------------------------------------------
static void assemble_A(std::vector<std::map<int,double>> &rows)
{
    // 介质参数现在是逐格点的：在每个循环点上从 g_rho_x / g_rho_z / g_C / g_L / g_M 取值。
    const double cax = dt / (2.0 * dx);
    const double caz = dt / (2.0 * dz);
    const double hd  = 0.5 * dt;          // 主方程里 ψ 的系数前缀（再乘 C/L/M）

    // ===========================================================
    // (A) 5 个主方程
    // ===========================================================

    // ---- vx 方程 ----  iz∈[0, Nz-1], ix∈[0, Nx]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_vx(iz, ix);
            const double rho_x = g_rho_x[lx_vx(iz, ix)];
            rows[row][id_vx(iz, ix)] += rho_x;

            // ∂x σxx 项：ix ≤ Nx-2
            if (ix <= Nx - 2) {
                rows[row][id_xx(iz, ix    )] += +cax;
                rows[row][id_xx(iz, ix + 1)] += -cax;
            }
            // ∂z σxz 项：iz ∈ [0, Nz-1]  (loop 内全合法)
            rows[row][id_xz(iz    , ix)] += +caz;
            rows[row][id_xz(iz + 1, ix)] += -caz;

            // PML ψ 耦合：-Δt/2 · ψ_xSxx, -Δt/2 · ψ_zSxz
            rows[row][id_pxSxx(iz, ix)] += -hd;
            rows[row][id_pzSxz(iz, ix)] += -hd;
        }
    }

    // ---- vz 方程 ----  iz∈[0, Nz], ix∈[0, Nx-1]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_vz(iz, ix);
            const double rho_z = g_rho_z[lx_vz(iz, ix)];
            rows[row][id_vz(iz, ix)] += rho_z;

            // ∂x σxz 项：ix ≥ 1, iz ≤ Nz-1
            if (ix >= 1 && iz <= Nz - 1) {
                rows[row][id_xz(iz + 1, ix - 1)] += +cax;
                rows[row][id_xz(iz + 1, ix    )] += -cax;
            }
            // ∂z σzz 项：iz ≤ Nz-2
            if (iz <= Nz - 2) {
                rows[row][id_zz(iz    , ix)] += +caz;
                rows[row][id_zz(iz + 1, ix)] += -caz;
            }

            rows[row][id_pxSxz(iz, ix)] += -hd;
            rows[row][id_pzSzz(iz, ix)] += -hd;
        }
    }

    // ---- σxx 方程 ----  iz∈[0, Nz-1], ix∈[0, Nx-1]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_xx(iz, ix);
            const double C = g_C[lx_xx(iz, ix)];
            const double L = g_L[lx_xx(iz, ix)];
            rows[row][id_xx(iz, ix)] += 1.0;

            if (ix >= 1) {
                rows[row][id_vx(iz, ix - 1)] += +cax * C;
                rows[row][id_vx(iz, ix    )] += -cax * C;
            }
            if (iz >= 1) {
                rows[row][id_vz(iz - 1, ix)] += +caz * L;
                rows[row][id_vz(iz    , ix)] += -caz * L;
            }

            rows[row][id_pxVx(iz, ix)] += -hd * C;
            rows[row][id_pzVz(iz, ix)] += -hd * L;
        }
    }

    // ---- σzz 方程 ----  iz∈[0, Nz-1], ix∈[0, Nx-1]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_zz(iz, ix);
            const double C = g_C[lx_xx(iz, ix)];
            const double L = g_L[lx_xx(iz, ix)];
            rows[row][id_zz(iz, ix)] += 1.0;

            if (ix >= 1) {
                rows[row][id_vx(iz, ix - 1)] += +cax * L;
                rows[row][id_vx(iz, ix    )] += -cax * L;
            }
            if (iz >= 1) {
                rows[row][id_vz(iz - 1, ix)] += +caz * C;
                rows[row][id_vz(iz    , ix)] += -caz * C;
            }

            rows[row][id_pxVx(iz, ix)] += -hd * L;
            rows[row][id_pzVz(iz, ix)] += -hd * C;
        }
    }

    // ---- σxz 方程 ----  iz∈[0, Nz], ix∈[0, Nx]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_xz(iz, ix);
            const double M = g_M[lx_xz(iz, ix)];
            rows[row][id_xz(iz, ix)] += 1.0;

            if (iz >= 1 && ix <= Nx - 2) {
                rows[row][id_vz(iz - 1, ix    )] += +cax * M;
                rows[row][id_vz(iz - 1, ix + 1)] += -cax * M;
            }
            if (iz >= 1 && iz <= Nz - 1) {
                rows[row][id_vx(iz - 1, ix)] += +caz * M;
                rows[row][id_vx(iz    , ix)] += -caz * M;
            }

            rows[row][id_pxVz(iz, ix)] += -hd * M;
            rows[row][id_pzVx(iz, ix)] += -hd * M;
        }
    }

    // ===========================================================
    // (B) 8 个 ψ 行块
    //
    //   通式：  1·ψ^{n+1}  +  (d·Δt/2)·(∂α field)^{n+1}  =  RHS
    //   其中 (∂α field)^{n+1} 用 1/h 的中心差分系数展开到对应 σ 或 v 的两个未知量上。
    //   边界跳过条件必须与对应主方程的 ∂α 完全相同。
    // ===========================================================
    const double half_dt = 0.5 * dt;

    // ---- ψ_xSxx ----  位置 vx，∂α = ∂x σxx，stencil (σxx[iz,ix], σxx[iz,ix+1])
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_pxSxx(iz, ix);
            rows[row][id_pxSxx(iz, ix)] += 1.0;

            const double d = g_dx_half[ix];
            if (d > 0.0 && ix <= Nx - 2) {
                // (d Δt/2) · (σxx[ix+1] - σxx[ix]) / dx
                const double coef = (half_dt * d) / dx;
                rows[row][id_xx(iz, ix    )] += -coef;
                rows[row][id_xx(iz, ix + 1)] += +coef;
            }
        }
    }

    // ---- ψ_zSxz ----  位置 vx，∂α = ∂z σxz，stencil (σxz[iz,ix], σxz[iz+1,ix])
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_pzSxz(iz, ix);
            rows[row][id_pzSxz(iz, ix)] += 1.0;

            const double d = g_dz_int[iz];
            if (d > 0.0) {
                const double coef = (half_dt * d) / dz;
                rows[row][id_xz(iz    , ix)] += -coef;
                rows[row][id_xz(iz + 1, ix)] += +coef;
            }
        }
    }

    // ---- ψ_xSxz ----  位置 vz，∂α = ∂x σxz，stencil (σxz[iz+1,ix-1], σxz[iz+1,ix])
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_pxSxz(iz, ix);
            rows[row][id_pxSxz(iz, ix)] += 1.0;

            const double d = g_dx_int[ix];
            if (d > 0.0 && ix >= 1 && iz <= Nz - 1) {
                const double coef = (half_dt * d) / dx;
                rows[row][id_xz(iz + 1, ix - 1)] += -coef;
                rows[row][id_xz(iz + 1, ix    )] += +coef;
            }
        }
    }

    // ---- ψ_zSzz ----  位置 vz，∂α = ∂z σzz，stencil (σzz[iz,ix], σzz[iz+1,ix])
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_pzSzz(iz, ix);
            rows[row][id_pzSzz(iz, ix)] += 1.0;

            const double d = g_dz_hv[iz];
            if (d > 0.0 && iz <= Nz - 2) {
                const double coef = (half_dt * d) / dz;
                rows[row][id_zz(iz    , ix)] += -coef;
                rows[row][id_zz(iz + 1, ix)] += +coef;
            }
        }
    }

    // ---- ψ_xVx ----  位置 σxx/σzz，∂α = ∂x vx，stencil (vx[iz,ix-1], vx[iz,ix])
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_pxVx(iz, ix);
            rows[row][id_pxVx(iz, ix)] += 1.0;

            const double d = g_dx_int[ix];
            if (d > 0.0 && ix >= 1) {
                const double coef = (half_dt * d) / dx;
                rows[row][id_vx(iz, ix - 1)] += -coef;
                rows[row][id_vx(iz, ix    )] += +coef;
            }
        }
    }

    // ---- ψ_zVz ----  位置 σxx/σzz，∂α = ∂z vz，stencil (vz[iz-1,ix], vz[iz,ix])
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            int row = id_pzVz(iz, ix);
            rows[row][id_pzVz(iz, ix)] += 1.0;

            const double d = g_dz_int[iz];
            if (d > 0.0 && iz >= 1) {
                const double coef = (half_dt * d) / dz;
                rows[row][id_vz(iz - 1, ix)] += -coef;
                rows[row][id_vz(iz    , ix)] += +coef;
            }
        }
    }

    // ---- ψ_xVz ----  位置 σxz，∂α = ∂x vz，stencil (vz[iz-1,ix], vz[iz-1,ix+1])
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_pxVz(iz, ix);
            rows[row][id_pxVz(iz, ix)] += 1.0;

            const double d = g_dx_half[ix];
            if (d > 0.0 && iz >= 1 && ix <= Nx - 2) {
                const double coef = (half_dt * d) / dx;
                rows[row][id_vz(iz - 1, ix    )] += -coef;
                rows[row][id_vz(iz - 1, ix + 1)] += +coef;
            }
        }
    }

    // ---- ψ_zVx ----  位置 σxz，∂α = ∂z vx，stencil (vx[iz-1,ix], vx[iz,ix])
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            int row = id_pzVx(iz, ix);
            rows[row][id_pzVx(iz, ix)] += 1.0;

            const double d = g_dz_hxz[iz];
            if (d > 0.0 && iz >= 1 && iz <= Nz - 1) {
                const double coef = (half_dt * d) / dz;
                rows[row][id_vx(iz - 1, ix)] += -coef;
                rows[row][id_vx(iz    , ix)] += +coef;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// 8. 右端向量 b 的组装：13 行块
//
//   主方程 RHS：在 v1 基础上多了 +(Δt/2)·ψ_old (或 +(Δt·C/2)·ψ_old 等)
//   ψ 行 RHS  ：E·(ψ_old - (d Δt/2)·(∂α)_old)
//
//   边界规则与 assemble_A 完全一致。声源加到 σxx, σzz 的 RHS。
// ----------------------------------------------------------------------------
static void assemble_b(std::vector<double> &b,
                       const std::vector<double> &vx_old,
                       const std::vector<double> &vz_old,
                       const std::vector<double> &txx_old,
                       const std::vector<double> &tzz_old,
                       const std::vector<double> &txz_old,
                       const std::vector<double> &pxSxx_old,
                       const std::vector<double> &pzSxz_old,
                       const std::vector<double> &pxSxz_old,
                       const std::vector<double> &pzSzz_old,
                       const std::vector<double> &pxVx_old,
                       const std::vector<double> &pzVz_old,
                       const std::vector<double> &pxVz_old,
                       const std::vector<double> &pzVx_old,
                       double sou_t)
{
    // 介质参数现在是逐格点的，从 g_rho_x / g_rho_z / g_C / g_L / g_M 取值
    const double cax   = dt / (2.0 * dx);
    const double caz   = dt / (2.0 * dz);
    const double hd    = 0.5 * dt;
    const double half_dt = 0.5 * dt;

    std::fill(b.begin(), b.end(), 0.0);

    // ---- (A) 主方程 RHS ----

    // vx 方程
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            const double rho_x = g_rho_x[lx_vx(iz, ix)];
            double val = rho_x * vx_old[lx_vx(iz, ix)];
            if (ix <= Nx - 2)
                val += cax * (txx_old[lx_xx(iz, ix + 1)] - txx_old[lx_xx(iz, ix)]);
            val += caz * (txz_old[lx_xz(iz + 1, ix)] - txz_old[lx_xz(iz, ix)]);
            val += hd * pxSxx_old[lx_vx(iz, ix)];
            val += hd * pzSxz_old[lx_vx(iz, ix)];
            b[id_vx(iz, ix)] = val;
        }
    }

    // vz 方程
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double rho_z = g_rho_z[lx_vz(iz, ix)];
            double val = rho_z * vz_old[lx_vz(iz, ix)];
            if (ix >= 1 && iz <= Nz - 1)
                val += cax * (txz_old[lx_xz(iz + 1, ix    )]
                            - txz_old[lx_xz(iz + 1, ix - 1)]);
            if (iz <= Nz - 2)
                val += caz * (tzz_old[lx_zz(iz + 1, ix)]
                            - tzz_old[lx_zz(iz    , ix)]);
            val += hd * pxSxz_old[lx_vz(iz, ix)];
            val += hd * pzSzz_old[lx_vz(iz, ix)];
            b[id_vz(iz, ix)] = val;
        }
    }

    // σxx 方程
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double C = g_C[lx_xx(iz, ix)];
            const double L = g_L[lx_xx(iz, ix)];
            double val = txx_old[lx_xx(iz, ix)];
            if (ix >= 1)
                val += cax * C * (vx_old[lx_vx(iz, ix    )]
                                - vx_old[lx_vx(iz, ix - 1)]);
            if (iz >= 1)
                val += caz * L * (vz_old[lx_vz(iz    , ix)]
                                - vz_old[lx_vz(iz - 1, ix)]);
            val += hd * C * pxVx_old[lx_xx(iz, ix)];
            val += hd * L * pzVz_old[lx_xx(iz, ix)];
            b[id_xx(iz, ix)] = val;
        }
    }

    // σzz 方程
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const double C = g_C[lx_xx(iz, ix)];
            const double L = g_L[lx_xx(iz, ix)];
            double val = tzz_old[lx_zz(iz, ix)];
            if (ix >= 1)
                val += cax * L * (vx_old[lx_vx(iz, ix    )]
                                - vx_old[lx_vx(iz, ix - 1)]);
            if (iz >= 1)
                val += caz * C * (vz_old[lx_vz(iz    , ix)]
                                - vz_old[lx_vz(iz - 1, ix)]);
            val += hd * L * pxVx_old[lx_zz(iz, ix)];
            val += hd * C * pzVz_old[lx_zz(iz, ix)];
            b[id_zz(iz, ix)] = val;
        }
    }

    // σxz 方程
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            const double M = g_M[lx_xz(iz, ix)];
            double val = txz_old[lx_xz(iz, ix)];
            if (iz >= 1 && ix <= Nx - 2)
                val += cax * M * (vz_old[lx_vz(iz - 1, ix + 1)]
                                - vz_old[lx_vz(iz - 1, ix    )]);
            if (iz >= 1 && iz <= Nz - 1)
                val += caz * M * (vx_old[lx_vx(iz    , ix)]
                                - vx_old[lx_vx(iz - 1, ix)]);
            val += hd * M * pxVz_old[lx_xz(iz, ix)];
            val += hd * M * pzVx_old[lx_xz(iz, ix)];
            b[id_xz(iz, ix)] = val;
        }
    }

    // 声源（应力点源，加到 σxx 与 σzz 的 RHS）
    b[id_xx(jsou, isou)] += dt * sou_t;
    b[id_zz(jsou, isou)] += dt * sou_t;

    // ---- (B) 8 个 ψ 行的 RHS ----
    //   通式：  RHS = E·(ψ_old  -  (d·Δt/2)·(∂α)_old)

    // ψ_xSxx：∂α = ∂x σxx,  位置 vx,  d = g_dx_half[ix],  E = g_ex_half[ix]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            double E = g_ex_half[ix];
            double d = g_dx_half[ix];
            double dax = 0.0;
            if (d > 0.0 && ix <= Nx - 2) {
                dax = (txx_old[lx_xx(iz, ix + 1)] - txx_old[lx_xx(iz, ix)]) / dx;
            }
            b[id_pxSxx(iz, ix)] = E * (pxSxx_old[lx_vx(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_zSxz：∂α = ∂z σxz, 位置 vx, d = g_dz_int[iz]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            double E = g_ez_int[iz];
            double d = g_dz_int[iz];
            double dax = 0.0;
            if (d > 0.0) {
                dax = (txz_old[lx_xz(iz + 1, ix)] - txz_old[lx_xz(iz, ix)]) / dz;
            }
            b[id_pzSxz(iz, ix)] = E * (pzSxz_old[lx_vx(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_xSxz：∂α = ∂x σxz, 位置 vz, d = g_dx_int[ix]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            double E = g_ex_int[ix];
            double d = g_dx_int[ix];
            double dax = 0.0;
            if (d > 0.0 && ix >= 1 && iz <= Nz - 1) {
                dax = (txz_old[lx_xz(iz + 1, ix    )]
                     - txz_old[lx_xz(iz + 1, ix - 1)]) / dx;
            }
            b[id_pxSxz(iz, ix)] = E * (pxSxz_old[lx_vz(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_zSzz：∂α = ∂z σzz, 位置 vz, d = g_dz_hv[iz]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            double E = g_ez_hv[iz];
            double d = g_dz_hv[iz];
            double dax = 0.0;
            if (d > 0.0 && iz <= Nz - 2) {
                dax = (tzz_old[lx_zz(iz + 1, ix)] - tzz_old[lx_zz(iz, ix)]) / dz;
            }
            b[id_pzSzz(iz, ix)] = E * (pzSzz_old[lx_vz(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_xVx：∂α = ∂x vx, 位置 σxx, d = g_dx_int[ix]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            double E = g_ex_int[ix];
            double d = g_dx_int[ix];
            double dax = 0.0;
            if (d > 0.0 && ix >= 1) {
                dax = (vx_old[lx_vx(iz, ix    )]
                     - vx_old[lx_vx(iz, ix - 1)]) / dx;
            }
            b[id_pxVx(iz, ix)] = E * (pxVx_old[lx_xx(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_zVz：∂α = ∂z vz, 位置 σxx, d = g_dz_int[iz]
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            double E = g_ez_int[iz];
            double d = g_dz_int[iz];
            double dax = 0.0;
            if (d > 0.0 && iz >= 1) {
                dax = (vz_old[lx_vz(iz    , ix)]
                     - vz_old[lx_vz(iz - 1, ix)]) / dz;
            }
            b[id_pzVz(iz, ix)] = E * (pzVz_old[lx_xx(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_xVz：∂α = ∂x vz, 位置 σxz, d = g_dx_half[ix]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            double E = g_ex_half[ix];
            double d = g_dx_half[ix];
            double dax = 0.0;
            if (d > 0.0 && iz >= 1 && ix <= Nx - 2) {
                dax = (vz_old[lx_vz(iz - 1, ix + 1)]
                     - vz_old[lx_vz(iz - 1, ix    )]) / dx;
            }
            b[id_pxVz(iz, ix)] = E * (pxVz_old[lx_xz(iz, ix)] - half_dt * d * dax);
        }
    }

    // ψ_zVx：∂α = ∂z vx, 位置 σxz, d = g_dz_hxz[iz]
    for (int iz = 0; iz <= Nz; iz++) {
        for (int ix = 0; ix <= Nx; ix++) {
            double E = g_ez_hxz[iz];
            double d = g_dz_hxz[iz];
            double dax = 0.0;
            if (d > 0.0 && iz >= 1 && iz <= Nz - 1) {
                dax = (vx_old[lx_vx(iz    , ix)]
                     - vx_old[lx_vx(iz - 1, ix)]) / dz;
            }
            b[id_pzVx(iz, ix)] = E * (pzVx_old[lx_xz(iz, ix)] - half_dt * d * dax);
        }
    }
}

// ----------------------------------------------------------------------------
// 9. 把解向量 x 拆回 13 个二维波场数组
// ----------------------------------------------------------------------------
static void unpack_solution(const std::vector<double> &x,
                            std::vector<double> &vx,
                            std::vector<double> &vz,
                            std::vector<double> &txx,
                            std::vector<double> &tzz,
                            std::vector<double> &txz,
                            std::vector<double> &pxSxx,
                            std::vector<double> &pzSxz,
                            std::vector<double> &pxSxz,
                            std::vector<double> &pzSzz,
                            std::vector<double> &pxVx,
                            std::vector<double> &pzVz,
                            std::vector<double> &pxVz,
                            std::vector<double> &pzVx)
{
    for (int i = 0; i < Nvx; i++) vx   [i] = x[off_vx    + i];
    for (int i = 0; i < Nvz; i++) vz   [i] = x[off_vz    + i];
    for (int i = 0; i < Nxx; i++) txx  [i] = x[off_xx    + i];
    for (int i = 0; i < Nzz; i++) tzz  [i] = x[off_zz    + i];
    for (int i = 0; i < Nxz; i++) txz  [i] = x[off_xz    + i];
    for (int i = 0; i < Nvx; i++) pxSxx[i] = x[off_pxSxx + i];
    for (int i = 0; i < Nvx; i++) pzSxz[i] = x[off_pzSxz + i];
    for (int i = 0; i < Nvz; i++) pxSxz[i] = x[off_pxSxz + i];
    for (int i = 0; i < Nvz; i++) pzSzz[i] = x[off_pzSzz + i];
    for (int i = 0; i < Nxx; i++) pxVx [i] = x[off_pxVx  + i];
    for (int i = 0; i < Nxx; i++) pzVz [i] = x[off_pzVz  + i];
    for (int i = 0; i < Nxz; i++) pxVz [i] = x[off_pxVz  + i];
    for (int i = 0; i < Nxz; i++) pzVx [i] = x[off_pzVx  + i];
}

// ----------------------------------------------------------------------------
// 10. wfile2d：与 Quantum/Implicit/main.cu 中签名一致，ASCII 输出
// ----------------------------------------------------------------------------
void wfile2d(const char* fn, const float* data, int nrows, int ncols)
{
    FILE* fp = fopen(fn, "wt");
    if (!fp) return;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++)
            fprintf(fp, "%e ", data[(size_t)i * (size_t)ncols + (size_t)j]);
        fprintf(fp, "\n");
    }
    fclose(fp);
}

static void wfile2d_d(const char* fn, const std::vector<double> &v, int nrows, int ncols)
{
    std::vector<float> tmp((size_t)nrows * ncols);
    for (size_t k = 0; k < tmp.size(); k++) tmp[k] = (float)v[k];
    wfile2d(fn, tmp.data(), nrows, ncols);
}

// ----------------------------------------------------------------------------
// 11. main
// ----------------------------------------------------------------------------
int main(void)
{
    // ---- 11.1 块尺寸 + 偏移 ----
    Nvx = Nz * (Nx + 1);
    Nvz = (Nz + 1) * Nx;
    Nxx = Nz * Nx;
    Nzz = Nz * Nx;
    Nxz = (Nz + 1) * (Nx + 1);

    off_vx    = 0;
    off_vz    = off_vx + Nvx;
    off_xx    = off_vz + Nvz;
    off_zz    = off_xx + Nxx;
    off_xz    = off_zz + Nzz;
    off_pxSxx = off_xz    + Nxz;
    off_pzSxz = off_pxSxx + Nvx;
    off_pxSxz = off_pzSxz + Nvx;
    off_pzSzz = off_pxSxz + Nvz;
    off_pxVx  = off_pzSzz + Nvz;
    off_pzVz  = off_pxVx  + Nxx;
    off_pxVz  = off_pzVz  + Nxx;
    off_pzVx  = off_pxVz  + Nxz;
    Ntot      = off_pzVx  + Nxz;

    printf("==== 2D CN-SGFD with ADE-PML (serial) ====\n");
    printf("Grid : Nz=%d, Nx=%d,  dx=%g, dz=%g,  dt=%g,  NT=%d,  NP=%d\n",
           Nz, Nx, dx, dz, dt, NT, NP);
    printf("Main block offsets : vx=%d, vz=%d, xx=%d, zz=%d, xz=%d\n",
           off_vx, off_vz, off_xx, off_zz, off_xz);
    printf("PML  block offsets : pxSxx=%d pzSxz=%d pxSxz=%d pzSzz=%d pxVx=%d pzVz=%d pxVz=%d pzVx=%d\n",
           off_pxSxx, off_pzSxz, off_pxSxz, off_pzSzz, off_pxVx, off_pzVz, off_pxVz, off_pzVx);
    printf("Total unknowns Ntot=%d  (main=%d, psi=%d)\n",
           Ntot, off_xz + Nxz, Ntot - (off_xz + Nxz));

    if (jrec >= Nz) jrec = Nz - 1;
    if (irec >= Nx) irec = Nx - 1;
    if (jrec < 0) jrec = 0;
    if (irec < 0) irec = 0;
    printf("Source (jsou,isou)=(%d,%d)   Recv (jrec,irec)=(%d,%d)\n",
           jsou, isou, jrec, irec);
    printf("Solid: Vp=%g  Vs=%g  rho=%g   (mu=%.4e, lambda=%.4e, lambda+2mu=%.4e)\n",
           Vp_val, Vs_val, rho_val, mu_val, lam_val, lam2u_val);
    printf("Fluid: Vp=%g  Vs=%g  rho=%g   (mu=%.4e, lambda=%.4e, lambda+2mu=%.4e)\n",
           Vp_fluid, Vs_fluid, rho_fluid, mu_fluid, lam_fluid, lam2u_fluid);

    // ---- 11.2 介质 + PML 阻尼曲线 ----
    build_medium();
    build_damping_profiles();
    printf("PML d0_x=%g  d0_z=%g   xleft=%g xright=%g  zleft=%g zright=%g\n",
           d0_x, d0_z, xleft_p, xright_p, zleft_p, zright_p);

    // ---- 11.3 分配 5 + 8 = 13 个旧时刻波场数组 ----
    std::vector<double> vx (Nvx, 0.0), vz (Nvz, 0.0);
    std::vector<double> txx(Nxx, 0.0), tzz(Nzz, 0.0), txz(Nxz, 0.0);

    std::vector<double> pxSxx(Nvx, 0.0), pzSxz(Nvx, 0.0);
    std::vector<double> pxSxz(Nvz, 0.0), pzSzz(Nvz, 0.0);
    std::vector<double> pxVx (Nxx, 0.0), pzVz (Nxx, 0.0);
    std::vector<double> pxVz (Nxz, 0.0), pzVx (Nxz, 0.0);

    std::vector<double> b(Ntot, 0.0);
    std::vector<double> xsol(Ntot, 0.0);

    // ---- 11.4 组装 A（均匀介质 + 常 PML 阻尼，A 不随时间变化）----
    printf("Assembling A ...\n");
    clock_t tA0 = clock();
    std::vector<std::map<int,double>> rows(Ntot);
    assemble_A(rows);
    CSRMatrix A;
    A.build_from_rows(Ntot, rows);
    rows.clear();
    rows.shrink_to_fit();
    printf("A: n=%d  nnz=%d  (assemble %.2f s)\n",
           A.n, (int)A.values.size(), (double)(clock()-tA0)/CLOCKS_PER_SEC);

    // README §16.2 验证（4×4 网格时检查 row 0 的 5 个主块 + 2 个 ψ 项）
    if (Nz == 4 && Nx == 4) {
        const double cax_exp = dt / (2.0 * dx);
        const double caz_exp = dt / (2.0 * dz);
        printf("[verify] Row 0 (vx eq at (0,0)) of A, expected per README §9:\n");
        for (int k = A.rowptr[0]; k < A.rowptr[1]; k++) {
            const char *hint = "";
            int c = A.colind[k];
            if      (c == 0)            hint = "  <- vx(0,0)    expect +rho_x";
            else if (c == 40)           hint = "  <- xx(0,0)    expect +cax";
            else if (c == 41)           hint = "  <- xx(0,1)    expect -cax";
            else if (c == 72)           hint = "  <- xz(0,0)    expect +caz";
            else if (c == 77)           hint = "  <- xz(1,0)    expect -caz";
            else if (c == off_pxSxx)    hint = "  <- pxSxx(0,0) PML, expect -dt/2";
            else if (c == off_pzSxz)    hint = "  <- pzSxz(0,0) PML, expect -dt/2";
            printf("  A[0,%4d] = % .6e %s\n", c, A.values[k], hint);
        }
        printf("  reference: +cax=%.6e  +caz=%.6e  -dt/2=%.6e\n",
               cax_exp, caz_exp, -0.5 * dt);
    }

    // ---- 11.5 输出容器 ----
    std::vector<float> sou_record(NT, 0.0f);
    std::vector<float> trace_txx (NT, 0.0f);
    std::vector<float> trace_tzz (NT, 0.0f);
    std::vector<float> trace_p   (NT, 0.0f);

    // 沿震源 z 轴的 σxx 时空记录：
    //   行优先 row-major，第 iz 行 = (iz, isou) 处 σxx 的所有时刻值；
    //   形状 Nz × NT，文件存储为 ASCII，每行 NT 个值。
    //   行编号方向：iz=0 是顶部，iz=Nz-1 是底部。
    std::vector<float> trace_zaxis_txx((size_t)Nz * (size_t)NT, 0.0f);

    // ---- 11.6 时间循环 ----
    const int    bicg_max = 5000;
#ifndef BICG_TOL
#define BICG_TOL 1.0e-9
#endif
    const double bicg_tol = BICG_TOL;

    clock_t t0 = clock();
    for (int it = 0; it < NT; it++) {
        double tt = (double)it * dt;
        double a = PI * F0 * (tt - T0);
        double sou_t = -(1.0 - 2.0 * a * a) * std::exp(-a * a);
        sou_record[it] = (float)sou_t;

        // warm start: 把上一步的全 13 块解作为 BiCGSTAB 的初值
        for (int i = 0; i < Nvx; i++) xsol[off_vx    + i] = vx   [i];
        for (int i = 0; i < Nvz; i++) xsol[off_vz    + i] = vz   [i];
        for (int i = 0; i < Nxx; i++) xsol[off_xx    + i] = txx  [i];
        for (int i = 0; i < Nzz; i++) xsol[off_zz    + i] = tzz  [i];
        for (int i = 0; i < Nxz; i++) xsol[off_xz    + i] = txz  [i];
        for (int i = 0; i < Nvx; i++) xsol[off_pxSxx + i] = pxSxx[i];
        for (int i = 0; i < Nvx; i++) xsol[off_pzSxz + i] = pzSxz[i];
        for (int i = 0; i < Nvz; i++) xsol[off_pxSxz + i] = pxSxz[i];
        for (int i = 0; i < Nvz; i++) xsol[off_pzSzz + i] = pzSzz[i];
        for (int i = 0; i < Nxx; i++) xsol[off_pxVx  + i] = pxVx [i];
        for (int i = 0; i < Nxx; i++) xsol[off_pzVz  + i] = pzVz [i];
        for (int i = 0; i < Nxz; i++) xsol[off_pxVz  + i] = pxVz [i];
        for (int i = 0; i < Nxz; i++) xsol[off_pzVx  + i] = pzVx [i];

        assemble_b(b, vx, vz, txx, tzz, txz,
                   pxSxx, pzSxz, pxSxz, pzSzz, pxVx, pzVz, pxVz, pzVx,
                   sou_t);

        double res = 0.0;
        int iters = bicgstab(A, b, xsol, bicg_max, bicg_tol, &res);

        unpack_solution(xsol, vx, vz, txx, tzz, txz,
                        pxSxx, pzSxz, pxSxz, pzSzz, pxVx, pzVz, pxVz, pzVx);

        int rec_idx_xx = lx_xx(jrec, irec);
        int rec_idx_zz = lx_zz(jrec, irec);
        trace_txx[it] = (float)txx[rec_idx_xx];
        trace_tzz[it] = (float)tzz[rec_idx_zz];
        trace_p  [it] = 0.5f * (trace_txx[it] + trace_tzz[it]);

        // 记录沿震源 z 轴 (x = isou) 整列 σxx：第 iz 行的第 it 列 = txx[iz, isou] @ time it
        for (int iz = 0; iz < Nz; iz++)
            trace_zaxis_txx[(size_t)iz * (size_t)NT + (size_t)it]
                = (float)txx[lx_xx(iz, isou)];

        if (it % 20 == 0) {
            printf("it=%4d  t=%.4f s  iters=%4d  res=%.3e  sou=%.3e  txx(rec)=%.3e\n",
                   it, tt, iters, res, sou_t, trace_txx[it]);
        }

        if (it % 50 == 0) {
            char fn[128];
            snprintf(fn, sizeof(fn), "snap_vx_it%04d.dat", it);
            wfile2d_d(fn, vx,  Nz,     Nx + 1);
            snprintf(fn, sizeof(fn), "snap_vz_it%04d.dat", it);
            wfile2d_d(fn, vz,  Nz + 1, Nx);
            snprintf(fn, sizeof(fn), "snap_txx_it%04d.dat", it);
            wfile2d_d(fn, txx, Nz,     Nx);
            snprintf(fn, sizeof(fn), "snap_tzz_it%04d.dat", it);
            wfile2d_d(fn, tzz, Nz,     Nx);
            snprintf(fn, sizeof(fn), "snap_txz_it%04d.dat", it);
            wfile2d_d(fn, txz, Nz + 1, Nx + 1);
        }
    }
    printf("Time loop done. Elapsed %.2f s\n", (double)(clock()-t0)/CLOCKS_PER_SEC);

    wfile2d("source_time.dat",   sou_record.data(), 1, NT);
    wfile2d("trace_txx.dat",     trace_txx.data(),  1, NT);
    wfile2d("trace_tzz.dat",     trace_tzz.data(),  1, NT);
    wfile2d("trace_pressure.dat",trace_p.data(),    1, NT);

    // 沿震源 z 轴的 σxx 时空记录：Nz 行 × NT 列
    //   行编号 iz=0..Nz-1 表示从顶到底；列编号 it=0..NT-1 表示时间
    //   读取后即可绘制成 z-t 图（VSP 风格）
    wfile2d("trace_zaxis_txx.dat", trace_zaxis_txx.data(), Nz, NT);
    printf("Saved trace_zaxis_txx.dat  shape = Nz x NT = %d x %d  (col x_src = %d)\n",
           Nz, NT, isou);

    printf("Done. snap_<field>_itXXXX.dat ; trace_*.dat ; trace_zaxis_txx.dat ; source_time.dat\n");
    return 0;
}
