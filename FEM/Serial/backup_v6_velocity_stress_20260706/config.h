/* =============================================================================
 *  config.h
 *  ---------------------------------------------------------------------------
 *  全局仿真参数。所有几何 / 声源 / 材料参数保持与
 *  Quantum/FDM/ImplicitV2/main.cpp 一致，以便 FEM 输出与 FDM 逐点对照。
 *
 *  与 FDM 的关键区别：
 *    - FEM 使用非交错的 P1 三角形连续单元，所有场（主 5 场 + 8 个 ψ 辅助场）
 *      全部定义在**同一套** Nx * Nz 网格节点上（node-centered）。
 *    - 因此 FEM 里所有 13 个块尺寸都是 Nnodes = Nx * Nz。
 *
 *  ---------------------------------------------------------------------------
 *  未知向量 q 的 13 块（按块偏移顺序）
 *  ---------------------------------------------------------------------------
 *      0  : vx          (Nnodes)
 *      1  : vz          (Nnodes)
 *      2  : σxx         (Nnodes)
 *      3  : σzz         (Nnodes)
 *      4  : σxz         (Nnodes)
 *      5  : ψ_xSxx      (Nnodes)   memory of ∂x σxx   (in vx eq)   d = d_x(node)
 *      6  : ψ_zSxz      (Nnodes)   memory of ∂z σxz   (in vx eq)   d = d_z(node)
 *      7  : ψ_xSxz      (Nnodes)   memory of ∂x σxz   (in vz eq)   d = d_x(node)
 *      8  : ψ_zSzz      (Nnodes)   memory of ∂z σzz   (in vz eq)   d = d_z(node)
 *      9  : ψ_xVx       (Nnodes)   memory of ∂x vx    (in σxx/σzz) d = d_x(node)
 *      10 : ψ_zVz       (Nnodes)   memory of ∂z vz    (in σxx/σzz) d = d_z(node)
 *      11 : ψ_xVz       (Nnodes)   memory of ∂x vz    (in σxz)     d = d_x(node)
 *      12 : ψ_zVx       (Nnodes)   memory of ∂z vx    (in σxz)     d = d_z(node)
 *
 *  可通过编译宏覆盖：
 *    -DNZ=200 -DNX=200 -DNT=4500
 *    -DDX_VAL=10.0 -DDZ_VAL=10.0 -DDT_VAL=1.0e-4
 *    -DF0_HZ=30.0 -DNBH=40
 *    -DNPMLDEF=30           PML 厚度（cells）；设为 0 完全关闭 PML，退化为旧行为
 *    -DRC_PML=1.0e-6        PML 目标反射系数
 * ===========================================================================*/
#ifndef FEM_CONFIG_H
#define FEM_CONFIG_H

#include <cmath>

// ----------------------------------------------------------------------------
// 1. 网格与时间参数
// ----------------------------------------------------------------------------
#ifndef NZ
#define NZ 200
#endif
#ifndef NX
#define NX 200
#endif
#ifndef NT
#define NT 4100
#endif
#ifndef DX_VAL
#define DX_VAL 10.0
#endif
#ifndef DZ_VAL
#define DZ_VAL 10.0
#endif
#ifndef DT_VAL
#define DT_VAL 1.0e-4
#endif

namespace femcfg {

constexpr int    Nx    = NX;
constexpr int    Nz    = NZ;
constexpr int    Ntime = NT;
constexpr double dx    = DX_VAL;
constexpr double dz    = DZ_VAL;
constexpr double dt    = DT_VAL;

constexpr int NcellX  = Nx - 1;
constexpr int NcellZ  = Nz - 1;
constexpr int Ncells  = NcellX * NcellZ;
constexpr int Ntri    = Ncells * 2;
constexpr int Nnodes  = Nx * Nz;

// 5 主场 + 8 个 ψ 辅助场 = 13 个块
constexpr int NfieldsPerNode = 13;
constexpr int Ntot           = Nnodes * NfieldsPerNode;

// 主场偏移
constexpr int off_vx     = 0  * Nnodes;
constexpr int off_vz     = 1  * Nnodes;
constexpr int off_xx     = 2  * Nnodes;
constexpr int off_zz     = 3  * Nnodes;
constexpr int off_xz     = 4  * Nnodes;

// PML ψ 偏移（8 个）
constexpr int off_pxSxx  = 5  * Nnodes;   // in vx eq, ∂x σxx
constexpr int off_pzSxz  = 6  * Nnodes;   // in vx eq, ∂z σxz
constexpr int off_pxSxz  = 7  * Nnodes;   // in vz eq, ∂x σxz
constexpr int off_pzSzz  = 8  * Nnodes;   // in vz eq, ∂z σzz
constexpr int off_pxVx   = 9  * Nnodes;   // in σxx/σzz, ∂x vx
constexpr int off_pzVz   = 10 * Nnodes;   // in σxx/σzz, ∂z vz
constexpr int off_pxVz   = 11 * Nnodes;   // in σxz,     ∂x vz
constexpr int off_pzVx   = 12 * Nnodes;   // in σxz,     ∂z vx

// ----------------------------------------------------------------------------
// 2. 介质参数（固体背景 + 中间高速层，与 FDM 一致）
// ----------------------------------------------------------------------------
#ifndef NBH
#define NBH 40
#endif
constexpr int NBH_val = NBH;
constexpr int borehole_iz_lo = 100;
constexpr int borehole_iz_hi = 120;

constexpr double Vp1  = 5000.0;
constexpr double Vs1  = 3300.0;
constexpr double rho1 = 2450.0;

constexpr double Vp2  = 6000.0;
constexpr double Vs2  = 3400.0;
constexpr double rho2 = 2770.0;

// ----------------------------------------------------------------------------
// 3. 声源（Ricker，与 FDM 相同）
// ----------------------------------------------------------------------------
#ifndef F0_HZ
#define F0_HZ 30.0
#endif
constexpr double PI = 3.14159265358979323846;
constexpr double F0 = F0_HZ;
constexpr double T0 = 1.2 / F0;

// 小网格下自动 clamp 声源 / 接收点位置
constexpr int _clamp01(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
constexpr int isou = _clamp01(NX / 2,          1, NX - 2);
constexpr int jsou = _clamp01(NZ / 2 - 40,     1, NZ - 2);
constexpr int irec = _clamp01(NX / 2,          1, NX - 2);
constexpr int jrec = _clamp01(NZ / 2 + 10,     1, NZ - 2);

inline double ricker(int it)
{
    const double tt = (double)it * dt;
    const double a  = PI * F0 * (tt - T0);
    return -(1.0 - 2.0 * a * a) * std::exp(-a * a);
}

// ----------------------------------------------------------------------------
// 4. PML 参数（与 FDM/ImplicitV2 完全一致的 d(α) 剖面）
// ----------------------------------------------------------------------------
#ifndef NPMLDEF
#define NPMLDEF 30
#endif
#ifndef RC_PML
#define RC_PML 1.0e-6
#endif

constexpr int    NP      = NPMLDEF;                    // PML 厚度（cells）
constexpr double RC_pml  = RC_PML;                     // 目标反射系数
constexpr double V0_pml  = Vp1;                        // 参考速度（用背景 P 波）
constexpr double DPx     = (double)NP * dx;
constexpr double DPz     = (double)NP * dz;

// d0 系数（NP=0 时关闭 PML）
constexpr double d0_x = (NP > 0)
                      ? (3.0 * V0_pml * 6.907755278982137 / (2.0 * DPx))  // ln(1/1e-6)≈6.9078
                      : 0.0;
constexpr double d0_z = (NP > 0)
                      ? (3.0 * V0_pml * 6.907755278982137 / (2.0 * DPz))
                      : 0.0;

// 更精确的 d0（构造函数里再算一遍确保和真实 -log(RC) 一致，见 assembler.cpp）
inline double compute_d0_x() { return 3.0 * V0_pml * std::log(1.0 / RC_pml) / (2.0 * DPx); }
inline double compute_d0_z() { return 3.0 * V0_pml * std::log(1.0 / RC_pml) / (2.0 * DPz); }

constexpr double xleft_p  = (double)NP           * dx;
constexpr double xright_p = (double)(Nx - 1 - NP) * dx;
constexpr double zleft_p  = (double)NP           * dz;
constexpr double zright_p = (double)(Nz - 1 - NP) * dz;

inline double dampx_at(double x)
{
    if (NP <= 0) return 0.0;
    const double d0 = compute_d0_x();
    if (x < xleft_p)  return d0 * std::pow((xleft_p  - x) / DPx, 2);
    if (x > xright_p) return d0 * std::pow((x - xright_p) / DPx, 2);
    return 0.0;
}
inline double dampz_at(double z)
{
    if (NP <= 0) return 0.0;
    const double d0 = compute_d0_z();
    if (z < zleft_p)  return d0 * std::pow((zleft_p  - z) / DPz, 2);
    if (z > zright_p) return d0 * std::pow((z - zright_p) / DPz, 2);
    return 0.0;
}

// ----------------------------------------------------------------------------
// 5. BiCGSTAB 求解器缺省参数
// ----------------------------------------------------------------------------
#ifndef BICG_MAX_IT
#define BICG_MAX_IT 5000
#endif
#ifndef BICG_TOL
#define BICG_TOL 1.0e-9
#endif
constexpr int    bicg_max_it = BICG_MAX_IT;
constexpr double bicg_tol    = BICG_TOL;

// ----------------------------------------------------------------------------
// 6. 输出控制
// ----------------------------------------------------------------------------
#ifndef SNAP_STRIDE
#define SNAP_STRIDE 50
#endif
constexpr int snap_stride = SNAP_STRIDE;

// ----------------------------------------------------------------------------
// 7. 常用索引：把 (ix,iz) 二维节点索引转换为 1D global node id
// ----------------------------------------------------------------------------
inline int node_id(int ix, int iz)   { return ix + iz * Nx; }
inline int dof_vx    (int ix, int iz){ return off_vx    + node_id(ix, iz); }
inline int dof_vz    (int ix, int iz){ return off_vz    + node_id(ix, iz); }
inline int dof_xx    (int ix, int iz){ return off_xx    + node_id(ix, iz); }
inline int dof_zz    (int ix, int iz){ return off_zz    + node_id(ix, iz); }
inline int dof_xz    (int ix, int iz){ return off_xz    + node_id(ix, iz); }
inline int dof_pxSxx (int ix, int iz){ return off_pxSxx + node_id(ix, iz); }
inline int dof_pzSxz (int ix, int iz){ return off_pzSxz + node_id(ix, iz); }
inline int dof_pxSxz (int ix, int iz){ return off_pxSxz + node_id(ix, iz); }
inline int dof_pzSzz (int ix, int iz){ return off_pzSzz + node_id(ix, iz); }
inline int dof_pxVx  (int ix, int iz){ return off_pxVx  + node_id(ix, iz); }
inline int dof_pzVz  (int ix, int iz){ return off_pzVz  + node_id(ix, iz); }
inline int dof_pxVz  (int ix, int iz){ return off_pxVz  + node_id(ix, iz); }
inline int dof_pzVx  (int ix, int iz){ return off_pzVx  + node_id(ix, iz); }

inline int cell_id(int icx, int icz) { return icx + icz * NcellX; }

} // namespace femcfg

#endif // FEM_CONFIG_H
