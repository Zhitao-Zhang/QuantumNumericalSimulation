/* =============================================================================
 *  config.h   —  位移二阶格式（SESSION 7 起，方案 B）
 *  ---------------------------------------------------------------------------
 *  【重大变更】丢弃速度–应力一阶 collocated 格式，改成**位移二阶格式**：
 *
 *      ρ ∂²u/∂t² = ∇·σ(u) + f,      σ = C : ε(u),   ε = ½(∇u + ∇uᵀ)
 *
 *  离散后（P1 三角，标准刚度阵）：
 *
 *      M ü + C u̇ + K u = f
 *          M = ∫ ρ N_i N_j        （consistent mass，2 dof/node）
 *          K = ∫ Bᵀ D B           （标准刚度阵，连中心节点 → 无奇偶解耦）
 *          C = 吸收层的质量比例阻尼（Rayleigh mass-proportional，仅边界层非零）
 *
 *  时间推进：Newmark 平均加速度（β=1/4, γ=1/2，无条件稳定、二阶、无人工阻尼），
 *      每步解  A a^{n+1} = b，A = M + γ·dt·C + β·dt²·K（与时间无关，装配一次）。
 *      保留 ILinearSolver / BiCGSTAB / CIM 落位口不变。
 *
 *  为什么换：velocity-stress 5 场全 collocated 在同一节点 → 一阶导等效 2Δx
 *  中心差分 → 奇偶解耦 + 各向异性频散 → 波前菱形（已用 diag_*.py 确诊）。
 *  位移格式里梯度以 ∫∇N·∇N 出现（连中心节点），各向异性天然小，波前是圆。
 *
 *  ---------------------------------------------------------------------------
 *  未知向量 q 的 2 块布局（按块偏移顺序）
 *  ---------------------------------------------------------------------------
 *      0 : ux   (Nnodes)     水平位移
 *      1 : uz   (Nnodes)     垂直位移
 *
 *  可通过编译宏覆盖：
 *    -DNZ=200 -DNX=200 -DNT=4100
 *    -DDX_VAL=10.0 -DDZ_VAL=10.0 -DDT_VAL=1.0e-4
 *    -DF0_HZ=30.0 -DNBH=40
 *    -DNABSDEF=30           吸收层厚度（cells）；设为 0 关闭吸收层
 *    -DRC_ABS=1.0e-3        吸收层目标残余
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

// ----------------------------------------------------------------------------
// 单元阶数：默认 P2（6 节点二次三角，抗频散）；-DUSE_P1_MESH 退回 P1（3 节点）
// ----------------------------------------------------------------------------
// 【P2 节点布局的关键设计】P2 的 6 个节点（3 顶点 + 3 边中点）恰好都落在一张
// 半间距 (dx/2) 的规则细网格上：
//     顶点     → 细网格偶数格点 (2*ix, 2*iz)
//     边中点   → 细网格半格点   (两端点索引平均)
// 细网格尺寸 MX×MZ = (2Nx-1)×(2Nz-1)，节点 id = IX + IZ*MX。
// 好处：PML 剖面按坐标算照常有效（IX=id%MX, IZ=id/MX，物理坐标 = IX*(dx/2)）；
//       边中点天然唯一编号，无需边表；快照只导出顶点子网格 (stride 2) = 原 Nx×Nz，
//       画图脚本与 FDM 逐点对照完全不用改。
#ifdef USE_P1_MESH
constexpr int    MX     = Nx;          // P1：节点网格 = 顶点网格
constexpr int    MZ     = Nz;
constexpr double hgrid  = dx;          // 节点网格间距 = 单元边长
constexpr double hgridz = dz;
constexpr int    node_stride = 1;      // 顶点在节点网格里的步长
#else
constexpr int    MX     = 2 * Nx - 1;  // P2：半间距细网格
constexpr int    MZ     = 2 * Nz - 1;
constexpr double hgrid  = 0.5 * dx;    // 细网格间距 = dx/2
constexpr double hgridz = 0.5 * dz;
constexpr int    node_stride = 2;      // 顶点落在细网格偶数格点
#endif

constexpr int Nnodes  = MX * MZ;       // 全局节点数（P2 时约为 P1 的 4 倍）
constexpr int Nvert   = Nx * Nz;       // 顶点子网格节点数（快照/对照用）

// 2 个位移场 = 2 个块
constexpr int NfieldsPerNode = 2;
constexpr int Ntot           = Nnodes * NfieldsPerNode;

// 块偏移
constexpr int off_ux = 0 * Nnodes;
constexpr int off_uz = 1 * Nnodes;

// 细网格 (IX,IZ) -> 全局 node id；物理坐标 = IX*hgrid, IZ*hgridz
inline int fine_node_id(int IX, int IZ) { return IX + IZ * MX; }
// 顶点 (ix,iz) -> 细网格 node id（顶点在偶数格点）
inline int vert_node_id(int ix, int iz) { return (ix * node_stride) + (iz * node_stride) * MX; }

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

// 声源为**各向同性矩张量**（爆炸源）：等效节点力 F_a = M0(t)·∇N_a。
//   M0(t) = SRC_SCALE · ricker(it)。SRC_SCALE 用 dx·dz（把 FDM 的
//   stress-glut 密度换算成矩），只影响整体振幅，不影响波前形状。
#ifndef SRC_SCALE
#define SRC_SCALE (DX_VAL * DZ_VAL)
#endif
constexpr double src_scale = SRC_SCALE;

// ----------------------------------------------------------------------------
// 4. Newmark 时间积分参数（平均加速度：无条件稳定、二阶、无人工阻尼）
// ----------------------------------------------------------------------------
#ifndef NEWMARK_BETA
#define NEWMARK_BETA 0.25
#endif
#ifndef NEWMARK_GAMMA
#define NEWMARK_GAMMA 0.5
#endif
constexpr double nm_beta  = NEWMARK_BETA;
constexpr double nm_gamma = NEWMARK_GAMMA;

// 混合质量系数 α：M = (1-α)·consistent + α·lumped，抗 P1 频散尾。
//   α=0 纯 consistent；理论最优抵消 P1 主导频散约 α=0.5（1/3-2/3 规则）。
#ifndef MASS_BLEND
#define MASS_BLEND 0.5
#endif
constexpr double mass_blend = MASS_BLEND;

// ----------------------------------------------------------------------------
// 5. ADE-PML 吸收边界（unsplit，复坐标拉伸；SESSION 8）
//   拉伸算子 s_α = 1 + d_α/(iω)，d_α(α) 为方向性阻尼曲线（1/s）。
//   位移二阶弱形式展开后（详见 assembler.cpp 顶部推导）：
//       M ü + C u̇ + K_tot u + F_ψ = F_src
//   其中：
//       C     = ∫ ρ(d_x + d_z) N_a N_b         （质量比例，对称）
//       K_tot = 标准刚度 K + 角区弹簧 P         （P = ∫ ρ d_x d_z N_a N_b，对称）
//       F_ψ   = ψ 记忆变量的 RHS 修正（显式更新，不进 A）
//   NABS=0 时 d_α≡0 → C=P=ψ=0，精确退化回无 PML 的纯位移格式。
//
//   d_α 剖面：多项式（二次）分级，d0 由目标反射系数 RC 反推（Collino-Tsogka）。
// ----------------------------------------------------------------------------
#ifndef NABSDEF
#define NABSDEF 30
#endif
#ifndef RC_ABS
#define RC_ABS 1.0e-4
#endif

constexpr int    NABS    = NABSDEF;                // PML 厚度（cells）
constexpr double RC_abs  = RC_ABS;                 // 目标理论反射系数
constexpr double V0_abs  = Vp1;                    // 参考速度（背景 P 波）
constexpr double DAx     = (double)NABS * dx;      // PML 物理厚度 x
constexpr double DAz     = (double)NABS * dz;

// d0 = -(p+1)·V0·ln(RC)/(2·L)，p=2（二次剖面），标准 PML 取值。
inline double compute_d0_x() { return (NABS > 0) ? (-3.0 * V0_abs * std::log(RC_abs) / (2.0 * DAx)) : 0.0; }
inline double compute_d0_z() { return (NABS > 0) ? (-3.0 * V0_abs * std::log(RC_abs) / (2.0 * DAz)) : 0.0; }

constexpr double xleft_a  = (double)NABS            * dx;
constexpr double xright_a = (double)(Nx - 1 - NABS) * dx;
constexpr double zleft_a  = (double)NABS            * dz;
constexpr double zright_a = (double)(Nz - 1 - NABS) * dz;

// 方向性阻尼曲线（PML 核心）：x 方向只看 x 坐标，z 方向只看 z 坐标。
inline double dampx_at(double x)
{
    if (NABS <= 0) return 0.0;
    const double d0 = compute_d0_x();
    if (x < xleft_a)  return d0 * std::pow((xleft_a  - x) / DAx, 2);
    if (x > xright_a) return d0 * std::pow((x - xright_a) / DAx, 2);
    return 0.0;
}
inline double dampz_at(double z)
{
    if (NABS <= 0) return 0.0;
    const double d0 = compute_d0_z();
    if (z < zleft_a)  return d0 * std::pow((zleft_a  - z) / DAz, 2);
    if (z > zright_a) return d0 * std::pow((z - zright_a) / DAz, 2);
    return 0.0;
}

// ----------------------------------------------------------------------------
// 6. BiCGSTAB 求解器缺省参数
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
// 7. 输出控制
// ----------------------------------------------------------------------------
#ifndef SNAP_STRIDE
#define SNAP_STRIDE 50
#endif
constexpr int snap_stride = SNAP_STRIDE;

// ----------------------------------------------------------------------------
// 8. 常用索引
//   注意：node_id 现在按顶点 (ix,iz) 映射到细网格全局 id（P2 时含 stride）。
//   dof 直接给全局 node id 加块偏移。
// ----------------------------------------------------------------------------
inline int node_id(int ix, int iz) { return vert_node_id(ix, iz); }
inline int dof_ux (int ix, int iz) { return off_ux + node_id(ix, iz); }
inline int dof_uz (int ix, int iz) { return off_uz + node_id(ix, iz); }
inline int cell_id(int icx, int icz) { return icx + icz * NcellX; }

} // namespace femcfg

#endif // FEM_CONFIG_H
