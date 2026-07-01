/* =============================================================================
 *  config.h
 *  ---------------------------------------------------------------------------
 *  全局仿真参数。所有几何 / 声源 / 材料参数刻意保持与
 *  Quantum/FDM/ImplicitV2/main.cpp 一致，以便 FEM 输出与 FDM 逐点对照。
 *
 *  与 FDM 的关键区别：
 *    - FEM 使用非交错的 P1 三角形连续单元，五个未知量 (vx, vz, sxx, szz, sxz)
 *      全部定义在**同一套** Nx * Nz 网格节点上（node-centered）。
 *    - 因此 FEM 不需要 (Nx+1)/(Nz+1) 半点网格；所有场都是 Nz × Nx。
 *    - PML 暂不启用（README §14），只保留后续扩展的接口。
 *
 *  可通过编译宏覆盖：
 *    -DNZ=200 -DNX=200 -DNT=4500
 *    -DDX_VAL=10.0 -DDZ_VAL=10.0 -DDT_VAL=1.0e-4
 *    -DF0_HZ=30.0 -DNBH=40
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
#define NT 1000
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

// 节点网格尺寸：FEM 中 Nx * Nz 个节点（与 FDM 应力网格同尺寸）。
constexpr int    Nx    = NX;
constexpr int    Nz    = NZ;
constexpr int    Ntime = NT;
constexpr double dx    = DX_VAL;
constexpr double dz    = DZ_VAL;
constexpr double dt    = DT_VAL;

// 单元(cell)数量与三角形数量：结构化 (Nx-1)*(Nz-1) 个矩形，每个切 2 个三角
constexpr int NcellX  = Nx - 1;
constexpr int NcellZ  = Nz - 1;
constexpr int Ncells  = NcellX * NcellZ;
constexpr int Ntri    = Ncells * 2;
constexpr int Nnodes  = Nx * Nz;

// 每个节点 5 个 DOF：vx, vz, sxx, szz, sxz
constexpr int NfieldsPerNode = 5;
constexpr int Ntot           = Nnodes * NfieldsPerNode;

// 全局 DOF 块偏移（按块堆叠布局）
constexpr int off_vx  = 0 * Nnodes;
constexpr int off_vz  = 1 * Nnodes;
constexpr int off_xx  = 2 * Nnodes;
constexpr int off_zz  = 3 * Nnodes;
constexpr int off_xz  = 4 * Nnodes;

// ----------------------------------------------------------------------------
// 2. 介质参数（固体背景 + 中间高速层，与 FDM 一致）
//    FDM 里 “borehole” 变量名虽写 borehole，但条件是 iz∈[100,120)，
//    实际是一条水平高速夹层。此处严格复刻。
// ----------------------------------------------------------------------------
#ifndef NBH
#define NBH 40
#endif
constexpr int NBH_val = NBH;
constexpr int borehole_iz_lo = 100;
constexpr int borehole_iz_hi = 120;

// 固体 1（背景）
constexpr double Vp1  = 5000.0;
constexpr double Vs1  = 3300.0;
constexpr double rho1 = 2450.0;

// 固体 2（水平夹层 iz∈[100,120)）
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

// 声源节点位置（FDM: jsou = Nz/2-40, isou = Nx/2；FEM 用同一 (ix,iz)）
//   对小网格自动 clamp，避免 make small 时 NZ/2-40 < 0
constexpr int _clamp01(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
constexpr int isou = _clamp01(NX / 2,          1, NX - 2);
constexpr int jsou = _clamp01(NZ / 2 - 40,     1, NZ - 2);

// 接收器节点位置
constexpr int irec = _clamp01(NX / 2,          1, NX - 2);
constexpr int jrec = _clamp01(NZ / 2 + 10,     1, NZ - 2);

inline double ricker(int it)
{
    const double tt = (double)it * dt;
    const double a  = PI * F0 * (tt - T0);
    return -(1.0 - 2.0 * a * a) * std::exp(-a * a);
}

// ----------------------------------------------------------------------------
// 4. BiCGSTAB 求解器缺省参数
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
// 5. 输出控制
// ----------------------------------------------------------------------------
#ifndef SNAP_STRIDE
#define SNAP_STRIDE 50
#endif
constexpr int snap_stride = SNAP_STRIDE;

// ----------------------------------------------------------------------------
// 6. 常用索引：把 (ix,iz) 二维节点索引转换为 1D global node id
//     节点 (ix, iz) 位置 = (ix*dx, iz*dz)
// ----------------------------------------------------------------------------
inline int node_id(int ix, int iz)   { return ix + iz * Nx; }
inline int dof_vx (int ix, int iz)   { return off_vx + node_id(ix, iz); }
inline int dof_vz (int ix, int iz)   { return off_vz + node_id(ix, iz); }
inline int dof_xx (int ix, int iz)   { return off_xx + node_id(ix, iz); }
inline int dof_zz (int ix, int iz)   { return off_zz + node_id(ix, iz); }
inline int dof_xz (int ix, int iz)   { return off_xz + node_id(ix, iz); }

// cell 索引：(icx, icz) → cell id
inline int cell_id(int icx, int icz) { return icx + icz * NcellX; }

} // namespace femcfg

#endif // FEM_CONFIG_H
