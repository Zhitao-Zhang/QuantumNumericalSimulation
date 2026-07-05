/* =============================================================================
 *  assembler.cpp   (支持 ADE-PML 的 13 块装配)
 *  ---------------------------------------------------------------------------
 *  质量矩阵：默认使用 **consistent mass**  (SESSION 5 关键修正)
 *
 *  背景：先前默认 lumped mass 的判断是错的 —— 对 velocity-stress **一阶** 系统
 *  用 lumped mass 会导致灾难性问题：
 *      σ 方程   :  M · ∂σ/∂t = D · v
 *      lumped 后：M 变成对角，即 σ_k 的更新是**纯 pointwise**；
 *      但 D · v 是有空间 stencil 的  → 两边空间尺度不匹配，
 *      σ 只在源 stencil 覆盖的少数节点上响应，波无法向外传播（诊断已确认）。
 *
 *  正确用法：
 *      - 对 velocity-stress **一阶** 系统：必须用 consistent mass（默认）
 *      - lumped mass 仅对 **displacement** 二阶系统 (M ü = -K u) 是安全的
 *
 *  编译期开关：
 *      -DUSE_LUMPED_MASS   诊断用（等价旧默认行为，波传不动）
 *      -DUSE_CONSISTENT_MASS  显式选择（本项已是默认）
 *
 *  声源施加：见 assemble_rhs 底部。FEM 弱形式下 `∫ N_k · δ(x-x_k)·f = f`，
 *  而 FDM 直接给节点/格点加 `dt·f`（等价于给出总量 f·(dx·dz) 的源）。为了让
 *  两者的**物理振幅**一致，FEM 里要把源乘以节点的“面积权重” M_L[k]。
 *
 *  P1 三角单元的闭式积分公式（|A_e| 为单元面积）：
 *
 *      ∫_e N_i N_j dA          = (|A_e|/12) * (1 + δ_ij)     (consistent)
 *                              = (|A_e|/3) * δ_ij            (lumped)
 *      ∫_e N_j dA              = |A_e|/3
 *      ∂N_i/∂x                 = 常数 (dN_dx[i])
 *      ∫_e (∂N_i/∂x) N_j dA    = dN_dx[i] * |A_e|/3
 *      ∫_e N_i (∂N_j/∂x) dA    = dN_dx[j] * |A_e|/3          (= 上式转置)
 *
 *  由此得到 3x3 局部矩阵：
 *      Mρ^e_ij  = ρ_e · M0^e_ij
 *      M_C^e_ij = C_e · M0^e_ij         (C = λ+2μ)
 *      M_L^e_ij = λ_e · M0^e_ij         (材料 λ；不要与节点 lumped mass 混淆)
 *      M_M^e_ij = μ_e · M0^e_ij
 *      Gx^e_ij  = dN_dx[i]·|A_e|/3     Gz^e_ij  = dN_dz[i]·|A_e|/3
 *      DxC^e_ij = C_e·dN_dx[j]·|A_e|/3  DzL^e_ij = λ_e·dN_dz[j]·|A_e|/3
 *      G_S,α^e_ij = dN_dα[j]·|A_e|/3    (用于 ψ 行的 stretched-coord 耦合)
 *
 *  ψ 方程（每节点 lumped mass 对角 + 与 field 的 G_S 耦合）：
 *      Row_k :  [M_L]_kk · ψ^{n+1}_k
 *              + (dt/2)·d_α(x_k)·Σ_j G_S_kj · field^{n+1}_j
 *              = E_α(x_k) · [ [M_L]_kk · ψ^n_k
 *                            - (dt/2)·d_α(x_k)·Σ_j G_S_kj · field^n_j ]
 *
 *  d_α = 0 处（PML 外）：E_α = 1，方程退化为 ψ^{n+1} = ψ^n，若 ψ^0=0 恒零。
 * ===========================================================================*/
#include "assembler.h"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>

namespace femasm {

using femcfg::Nx;
using femcfg::Nz;
using femcfg::Nnodes;
using femcfg::Ntot;
using femcfg::off_vx;   using femcfg::off_vz;
using femcfg::off_xx;   using femcfg::off_zz;   using femcfg::off_xz;
using femcfg::off_pxSxx; using femcfg::off_pzSxz;
using femcfg::off_pxSxz; using femcfg::off_pzSzz;
using femcfg::off_pxVx;  using femcfg::off_pzVz;
using femcfg::off_pxVz;  using femcfg::off_pzVx;
using femcfg::isou;      using femcfg::jsou;
using femcfg::dx;        using femcfg::dz;
using femcfg::NP;

using femla::RowMap;
using femla::SparseMatrixCSR;
using femmesh::Triangle;

static inline int gid_vx    (int node) { return off_vx    + node; }
static inline int gid_vz    (int node) { return off_vz    + node; }
static inline int gid_xx    (int node) { return off_xx    + node; }
static inline int gid_zz    (int node) { return off_zz    + node; }
static inline int gid_xz    (int node) { return off_xz    + node; }
static inline int gid_pxSxx (int node) { return off_pxSxx + node; }
static inline int gid_pzSxz (int node) { return off_pzSxz + node; }
static inline int gid_pxSxz (int node) { return off_pxSxz + node; }
static inline int gid_pzSzz (int node) { return off_pzSzz + node; }
static inline int gid_pxVx  (int node) { return off_pxVx  + node; }
static inline int gid_pzVz  (int node) { return off_pzVz  + node; }
static inline int gid_pxVz  (int node) { return off_pxVz  + node; }
static inline int gid_pzVx  (int node) { return off_pzVx  + node; }

// 从 node 全局 id 反解出 ix, iz（用于查节点级 d_α）
static inline void node_ij(int node, int &ix, int &iz)
{
    iz = node / Nx;
    ix = node - iz * Nx;
}

// -----------------------------------------------------------------------------
// 预计算 PML 阻尼曲线（节点级采样）
// -----------------------------------------------------------------------------
static void build_pml_profiles(double dt, FemSystem &sys)
{
    sys.dampx_at_ix.assign(Nx, 0.0);
    sys.dampz_at_iz.assign(Nz, 0.0);
    sys.Ex_at_ix   .assign(Nx, 1.0);
    sys.Ez_at_iz   .assign(Nz, 1.0);

    if (NP <= 0) {
        printf("[pml] NP=0  -> PML disabled (all d=0, E=1)\n");
        sys.pml_active_x = sys.pml_active_z = 0;
        return;
    }

    int nx_active = 0, nz_active = 0;
    for (int ix = 0; ix < Nx; ix++) {
        sys.dampx_at_ix[ix] = femcfg::dampx_at(ix * dx);
        sys.Ex_at_ix   [ix] = std::exp(-sys.dampx_at_ix[ix] * dt);
        if (sys.dampx_at_ix[ix] > 0.0) nx_active++;
    }
    for (int iz = 0; iz < Nz; iz++) {
        sys.dampz_at_iz[iz] = femcfg::dampz_at(iz * dz);
        sys.Ez_at_iz   [iz] = std::exp(-sys.dampz_at_iz[iz] * dt);
        if (sys.dampz_at_iz[iz] > 0.0) nz_active++;
    }
    sys.pml_active_x = nx_active;
    sys.pml_active_z = nz_active;

    double dxmax = 0.0, dzmax = 0.0;
    for (double v : sys.dampx_at_ix) if (v > dxmax) dxmax = v;
    for (double v : sys.dampz_at_iz) if (v > dzmax) dzmax = v;
    printf("[pml] NP=%d  d0_x=%g  d0_z=%g\n",
           NP, femcfg::compute_d0_x(), femcfg::compute_d0_z());
    printf("[pml] xleft=%g xright=%g  zleft=%g zright=%g\n",
           femcfg::xleft_p, femcfg::xright_p, femcfg::zleft_p, femcfg::zright_p);
    printf("[pml] active x-columns=%d/%d  z-rows=%d/%d   max(dampx)=%g  max(dampz)=%g\n",
           nx_active, Nx, nz_active, Nz, dxmax, dzmax);
}

// -----------------------------------------------------------------------------
// 装配 A 与 M
// -----------------------------------------------------------------------------
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  FemSystem &out)
{
#ifdef USE_Q1_MESH
    printf("[assemble] Ntot=%d  Nquad=%d  (Q1 bilinear) ...\n", Ntot, mesh.n_quads());
#else
    printf("[assemble] Ntot=%d  Ntri=%d  (P1 tri Union Jack) ...\n", Ntot, mesh.n_tris());
#endif
    fflush(stdout);
    clock_t tA0 = clock();

    // ---- PML 曲线 ----
    build_pml_profiles(dt, out);

    // 便捷指针
    const std::vector<double> &dampx = out.dampx_at_ix;
    const std::vector<double> &dampz = out.dampz_at_iz;

    std::vector<RowMap> rowsA(Ntot);
    std::vector<RowMap> rowsM(Ntot);

    // 时间耦合系数 hc：
    //   默认 Crank–Nicolson  hc = dt/2  （与 FDM/ImplicitV2 一致，A = M + (dt/2)K）
    //   -DUSE_BACKWARD_EULER  hc = dt    （旧版 Backward Euler，A = M + dt·K）
#ifdef USE_BACKWARD_EULER
    const double hc = dt;
    printf("[time] Backward Euler  (hc=dt=%g)\n", dt);
#else
    const double hc = 0.5 * dt;
    printf("[time] Crank-Nicolson  (hc=dt/2=%g)  — matches FDM/ImplicitV2\n", hc);
#endif

#ifdef USE_Q1_MESH
    // =============================================================
    //   Q1 四边形装配（每个 cell 一个 4-node 双线性单元）
    //
    //   本地节点编号：0=SW, 1=SE, 2=NE, 3=NW  (CCW)
    //   在规则矩形 dx × dz 单元上，闭式解析元素矩阵：
    //
    //     M0[i][j]  (consistent)  = (dx·dz / 36) · mMc[i][j]
    //     mMc = |4 2 1 2|         每行和 = 9 → row-sum diag = dx·dz/4
    //           |2 4 2 1|
    //           |1 2 4 2|
    //           |2 1 2 4|
    //
    //     Dx[i][j]  = ∫∫ N_i ∂N_j/∂x  = (dz/12) · mDx[i][j]
    //     mDx = |-2 +2 +1 -1|         (行 0,1 共享；行 2,3 共享)
    //           |-2 +2 +1 -1|
    //           |-1 +1 +2 -2|
    //           |-1 +1 +2 -2|
    //
    //     Dz[i][j]  = ∫∫ N_i ∂N_j/∂z  = (dx/12) · mDz[i][j]
    //     mDz = |-2 -1 +1 +2|
    //           |-1 -2 +2 +1|
    //           |-1 -2 +2 +1|
    //           |-2 -1 +1 +2|
    //
    //   注意 Dx ≠ Dx^T：Q1 元里 ∂N/∂x 不是常数（依赖 z），因此
    //   `∫ N_i ∂N_j/∂x` 与 `∫ (∂N_i/∂x) N_j` 是不同矩阵（互为转置）。
    //   动量方程走分部积分用 Dx^T = Dx[j][i]，本构方程不 IBP 用 Dx[i][j]。
    // =============================================================
    // Q1 element matrices (analytical, dx·dz rectangle)
    //   consistent M0 = (|A|/36) · mMc     row sum = 9 → row-sum diag = |A|/4
    static const int mMc [4][4] = { {4,2,1,2},{2,4,2,1},{1,2,4,2},{2,1,2,4} };
    static const int mDx [4][4] = { {-2,2,1,-1},{-2,2,1,-1},{-1,1,2,-2},{-1,1,2,-2} };
    static const int mDz [4][4] = { {-2,-1,1,2},{-1,-2,2,1},{-1,-2,2,1},{-2,-1,1,2} };

    for (const auto &qe : mesh.quads) {
        const femmat::CellProps &mp = mat.at_cell(qe.cell_id);
        const double rho_e = mp.rho;
        const double lam_e = mp.lambda;
        const double mu_e  = mp.mu;
        const double C_e   = lam_e + 2.0 * mu_e;

        const double dxq = qe.dx_val;
        const double dzq = qe.dz_val;
        const double A_e = dxq * dzq;
        const double A4  = A_e / 4.0;
        const double kDx = dzq / 12.0;
        const double kDz = dxq / 12.0;
        (void)A4;   // may be unused in consistent branch

        // Q1 4x4 mass 模板 —— **默认 consistent**（SESSION 5 修正）
#ifdef USE_LUMPED_MASS
        auto m0_ = [&](int i, int j) -> double { return (i == j) ? A4 : 0.0; };
#else
        const double A36 = A_e / 36.0;
        auto m0_ = [&](int i, int j) -> double { return A36 * (double)mMc[i][j]; };
#endif
        auto Dxq = [&](int i, int j) -> double { return kDx * (double)mDx[i][j]; };
        auto Dzq = [&](int i, int j) -> double { return kDz * (double)mDz[i][j]; };

        const int gnode[4] = { qe.nodes[0], qe.nodes[1], qe.nodes[2], qe.nodes[3] };

        // ---- (A) 主 5 方程的 mass + K 耦合 + 主-ψ 耦合 ----
        for (int i = 0; i < 4; i++) {
            const int gi = gnode[i];
            for (int j = 0; j < 4; j++) {
                const int gj = gnode[j];
                const double mij   = m0_(i, j);
                const double mij_r = rho_e * mij;

                rowsA[gid_vx(gi)][gid_vx(gj)] += mij_r;
                rowsA[gid_vz(gi)][gid_vz(gj)] += mij_r;
                rowsA[gid_xx(gi)][gid_xx(gj)] += mij;
                rowsA[gid_zz(gi)][gid_zz(gj)] += mij;
                rowsA[gid_xz(gi)][gid_xz(gj)] += mij;

                // 动量方程（IBP）→ 系数用 Dx^T[i][j] = Dx[j][i]
                rowsA[gid_vx(gi)][gid_xx(gj)] += +hc * Dxq(j, i);
                rowsA[gid_vx(gi)][gid_xz(gj)] += +hc * Dzq(j, i);
                rowsA[gid_vz(gi)][gid_xz(gj)] += +hc * Dxq(j, i);
                rowsA[gid_vz(gi)][gid_zz(gj)] += +hc * Dzq(j, i);

                // 本构方程（不 IBP）→ 系数用 Dx[i][j]
                rowsA[gid_xx(gi)][gid_vx(gj)] += -hc * C_e   * Dxq(i, j);
                rowsA[gid_xx(gi)][gid_vz(gj)] += -hc * lam_e * Dzq(i, j);
                rowsA[gid_zz(gi)][gid_vx(gj)] += -hc * lam_e * Dxq(i, j);
                rowsA[gid_zz(gi)][gid_vz(gj)] += -hc * C_e   * Dzq(i, j);
                rowsA[gid_xz(gi)][gid_vx(gj)] += -hc * mu_e  * Dzq(i, j);
                rowsA[gid_xz(gi)][gid_vz(gj)] += -hc * mu_e  * Dxq(i, j);

                // 主-ψ 耦合：-hc·M
                rowsA[gid_vx(gi)][gid_pxSxx(gj)] += -hc * mij;
                rowsA[gid_vx(gi)][gid_pzSxz(gj)] += -hc * mij;
                rowsA[gid_vz(gi)][gid_pxSxz(gj)] += -hc * mij;
                rowsA[gid_vz(gi)][gid_pzSzz(gj)] += -hc * mij;
                rowsA[gid_xx(gi)][gid_pxVx (gj)] += -hc * C_e   * mij;
                rowsA[gid_xx(gi)][gid_pzVz (gj)] += -hc * lam_e * mij;
                rowsA[gid_zz(gi)][gid_pxVx (gj)] += -hc * lam_e * mij;
                rowsA[gid_zz(gi)][gid_pzVz (gj)] += -hc * C_e   * mij;
                rowsA[gid_xz(gi)][gid_pzVx (gj)] += -hc * mu_e  * mij;
                rowsA[gid_xz(gi)][gid_pxVz (gj)] += -hc * mu_e  * mij;

                // M 矩阵（主 5 场）
                rowsM[gid_vx(gi)][gid_vx(gj)] += mij_r;
                rowsM[gid_vz(gi)][gid_vz(gj)] += mij_r;
                rowsM[gid_xx(gi)][gid_xx(gj)] += mij;
                rowsM[gid_zz(gi)][gid_zz(gj)] += mij;
                rowsM[gid_xz(gi)][gid_xz(gj)] += mij;
            }
        }

        // ---- (B) 8 个 ψ 行：lumped mass 对角 + 与 field 的耦合 ----
        for (int i = 0; i < 4; i++) {
            const int gi = gnode[i];
            int ix_i, iz_i;
            node_ij(gi, ix_i, iz_i);
            const double dx_i = dampx[ix_i];
            const double dz_i = dampz[iz_i];

            rowsA[gid_pxSxx(gi)][gid_pxSxx(gi)] += A4;
            rowsA[gid_pzSxz(gi)][gid_pzSxz(gi)] += A4;
            rowsA[gid_pxSxz(gi)][gid_pxSxz(gi)] += A4;
            rowsA[gid_pzSzz(gi)][gid_pzSzz(gi)] += A4;
            rowsA[gid_pxVx (gi)][gid_pxVx (gi)] += A4;
            rowsA[gid_pzVz (gi)][gid_pzVz (gi)] += A4;
            rowsA[gid_pxVz (gi)][gid_pxVz (gi)] += A4;
            rowsA[gid_pzVx (gi)][gid_pzVx (gi)] += A4;

            rowsM[gid_pxSxx(gi)][gid_pxSxx(gi)] += A4;
            rowsM[gid_pzSxz(gi)][gid_pzSxz(gi)] += A4;
            rowsM[gid_pxSxz(gi)][gid_pxSxz(gi)] += A4;
            rowsM[gid_pzSzz(gi)][gid_pzSzz(gi)] += A4;
            rowsM[gid_pxVx (gi)][gid_pxVx (gi)] += A4;
            rowsM[gid_pzVz (gi)][gid_pzVz (gi)] += A4;
            rowsM[gid_pxVz (gi)][gid_pxVz (gi)] += A4;
            rowsM[gid_pzVx (gi)][gid_pzVx (gi)] += A4;

            const double hd_x = 0.5 * dt * dx_i;
            const double hd_z = 0.5 * dt * dz_i;
            if (hd_x != 0.0 || hd_z != 0.0) {
                for (int j = 0; j < 4; j++) {
                    const int gj = gnode[j];
                    const double gxj = Dxq(i, j);   // ∫ N_i ∂N_j/∂x   (Q1: 依赖 i)
                    const double gzj = Dzq(i, j);
                    if (hd_x != 0.0) {
                        rowsA[gid_pxSxx(gi)][gid_xx(gj)] += hd_x * gxj;
                        rowsA[gid_pxSxz(gi)][gid_xz(gj)] += hd_x * gxj;
                        rowsA[gid_pxVx (gi)][gid_vx(gj)] += hd_x * gxj;
                        rowsA[gid_pxVz (gi)][gid_vz(gj)] += hd_x * gxj;
                    }
                    if (hd_z != 0.0) {
                        rowsA[gid_pzSxz(gi)][gid_xz(gj)] += hd_z * gzj;
                        rowsA[gid_pzSzz(gi)][gid_zz(gj)] += hd_z * gzj;
                        rowsA[gid_pzVz (gi)][gid_vz(gj)] += hd_z * gzj;
                        rowsA[gid_pzVx (gi)][gid_vx(gj)] += hd_z * gzj;
                    }
                }
            }
        }
    }
#else
    // (节点级 lumped mass 在 element loop 中通过累加 A_e/3 得到，无需显式引用 mesh.M_L)

    // ---- 遍历每个三角形，装配全部 3x3 局部贡献 ----
    for (const auto &t : mesh.triangles) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double rho_e = mp.rho;
        const double lam_e = mp.lambda;
        const double mu_e  = mp.mu;
        const double C_e   = lam_e + 2.0 * mu_e;

        const double A_e     = t.A;
        const double A3      = A_e / 3.0;
        const double dNdx[3] = { t.dN_dx[0], t.dN_dx[1], t.dN_dx[2] };
        const double dNdz[3] = { t.dN_dz[0], t.dN_dz[1], t.dN_dz[2] };
        const int gnode[3] = { t.nodes[0], t.nodes[1], t.nodes[2] };

        // P1 tri 3x3 mass 模板 —— **默认 consistent**（SESSION 5 修正）
        //     consistent : M0_ij = (|A|/12)·(1 + δ_ij)
        //     lumped     : M0_ij = (|A|/3)·δ_ij  ← 仅诊断用，对 velocity-stress 会破坏波传播
        (void)A3;   // suppress unused warning in consistent branch (A3 in ψ block)
#ifdef USE_LUMPED_MASS
        auto m0_ = [&](int i, int j) -> double { return (i == j) ? A3 : 0.0; };
#else
        const double A12 = A_e / 12.0;
        auto m0_ = [&](int i, int j) -> double { return A12 * (i == j ? 2.0 : 1.0); };
#endif

        // ---------------------------------------------------------------
        // (A) 主 5 方程的主 mass 与 K 耦合（同旧版）+  PML ψ 耦合项
        // ---------------------------------------------------------------
        for (int i = 0; i < 3; i++) {
            const int gi = gnode[i];

            for (int j = 0; j < 3; j++) {
                const int gj = gnode[j];
                const double mij = m0_(i, j);
                const double mij_r = rho_e * mij;

                // 主方程质量对角块
                rowsA[gid_vx(gi)][gid_vx(gj)] += mij_r;
                rowsA[gid_vz(gi)][gid_vz(gj)] += mij_r;
                rowsA[gid_xx(gi)][gid_xx(gj)] += mij;
                rowsA[gid_zz(gi)][gid_zz(gj)] += mij;
                rowsA[gid_xz(gi)][gid_xz(gj)] += mij;

                // 主方程 K 耦合 (dt · G / D 型)
                //   vx 行： +dt·Gx σxx + dt·Gz σxz
                //         Gx_ij = dN_dx[i] · A/3
                rowsA[gid_vx(gi)][gid_xx(gj)] += +hc * dNdx[i] * A3;
                rowsA[gid_vx(gi)][gid_xz(gj)] += +hc * dNdz[i] * A3;
                //   vz 行： +hc·Gx σxz + hc·Gz σzz
                rowsA[gid_vz(gi)][gid_xz(gj)] += +hc * dNdx[i] * A3;
                rowsA[gid_vz(gi)][gid_zz(gj)] += +hc * dNdz[i] * A3;
                //   σxx 行： -hc·DxC vx - hc·DzL vz
                rowsA[gid_xx(gi)][gid_vx(gj)] += -hc * C_e   * dNdx[j] * A3;
                rowsA[gid_xx(gi)][gid_vz(gj)] += -hc * lam_e * dNdz[j] * A3;
                //   σzz 行： -hc·DxL vx - hc·DzC vz
                rowsA[gid_zz(gi)][gid_vx(gj)] += -hc * lam_e * dNdx[j] * A3;
                rowsA[gid_zz(gi)][gid_vz(gj)] += -hc * C_e   * dNdz[j] * A3;
                //   σxz 行： -hc·DzM vx - hc·DxM vz
                rowsA[gid_xz(gi)][gid_vx(gj)] += -hc * mu_e  * dNdz[j] * A3;
                rowsA[gid_xz(gi)][gid_vz(gj)] += -hc * mu_e  * dNdx[j] * A3;

                // 主方程与 ψ 的耦合（stretched-coord 修正）
                rowsA[gid_vx(gi)][gid_pxSxx(gj)] += -hc * mij;
                rowsA[gid_vx(gi)][gid_pzSxz(gj)] += -hc * mij;

                rowsA[gid_vz(gi)][gid_pxSxz(gj)] += -hc * mij;
                rowsA[gid_vz(gi)][gid_pzSzz(gj)] += -hc * mij;

                rowsA[gid_xx(gi)][gid_pxVx (gj)] += -hc * C_e   * mij;
                rowsA[gid_xx(gi)][gid_pzVz (gj)] += -hc * lam_e * mij;

                rowsA[gid_zz(gi)][gid_pxVx (gj)] += -hc * lam_e * mij;
                rowsA[gid_zz(gi)][gid_pzVz (gj)] += -hc * C_e   * mij;

                rowsA[gid_xz(gi)][gid_pzVx (gj)] += -hc * mu_e  * mij;
                rowsA[gid_xz(gi)][gid_pxVz (gj)] += -hc * mu_e  * mij;

                // ----- M 矩阵（用于每步 RHS 的 M·q^n） -----
                // 主 5 场：与旧版相同
                rowsM[gid_vx(gi)][gid_vx(gj)] += mij_r;
                rowsM[gid_vz(gi)][gid_vz(gj)] += mij_r;
                rowsM[gid_xx(gi)][gid_xx(gj)] += mij;
                rowsM[gid_zz(gi)][gid_zz(gj)] += mij;
                rowsM[gid_xz(gi)][gid_xz(gj)] += mij;
            }
        }

        // ---------------------------------------------------------------
        // (B) 8 个 ψ 行：
        //     A[ψ_k, ψ_k]   += A_e/3     (lumped mass, 累积 element 贡献)
        //     A[ψ_k, field_j] += (dt/2)·d_α(x_k)·dNdα^e[j]·A_e/3
        //
        //     M 中 ψ 行只有 M_L 对角（一样通过 A_e/3 累积）
        //     RHS 每步再乘 E_α(x_k) 并加梯度修正
        // ---------------------------------------------------------------
        for (int i = 0; i < 3; i++) {
            const int gi = gnode[i];
            int ix_i, iz_i;
            node_ij(gi, ix_i, iz_i);
            const double dx_i = dampx[ix_i];
            const double dz_i = dampz[iz_i];

            // ψ 对角（M_L 的 element 贡献）
            rowsA[gid_pxSxx(gi)][gid_pxSxx(gi)] += A3;
            rowsA[gid_pzSxz(gi)][gid_pzSxz(gi)] += A3;
            rowsA[gid_pxSxz(gi)][gid_pxSxz(gi)] += A3;
            rowsA[gid_pzSzz(gi)][gid_pzSzz(gi)] += A3;
            rowsA[gid_pxVx (gi)][gid_pxVx (gi)] += A3;
            rowsA[gid_pzVz (gi)][gid_pzVz (gi)] += A3;
            rowsA[gid_pxVz (gi)][gid_pxVz (gi)] += A3;
            rowsA[gid_pzVx (gi)][gid_pzVx (gi)] += A3;

            rowsM[gid_pxSxx(gi)][gid_pxSxx(gi)] += A3;
            rowsM[gid_pzSxz(gi)][gid_pzSxz(gi)] += A3;
            rowsM[gid_pxSxz(gi)][gid_pxSxz(gi)] += A3;
            rowsM[gid_pzSzz(gi)][gid_pzSzz(gi)] += A3;
            rowsM[gid_pxVx (gi)][gid_pxVx (gi)] += A3;
            rowsM[gid_pzVz (gi)][gid_pzVz (gi)] += A3;
            rowsM[gid_pxVz (gi)][gid_pxVz (gi)] += A3;
            rowsM[gid_pzVx (gi)][gid_pzVx (gi)] += A3;

            // ψ - field 耦合。只有在 d > 0 时才产生非零系数
            const double hd_x = 0.5 * dt * dx_i;
            const double hd_z = 0.5 * dt * dz_i;
            if (hd_x != 0.0 || hd_z != 0.0) {
                for (int j = 0; j < 3; j++) {
                    const int gj = gnode[j];
                    const double gxj = dNdx[j] * A3;   // G_S,x_ij = ∫ N_i ∂N_j/∂x  (行 i 无关)
                    const double gzj = dNdz[j] * A3;

                    if (hd_x != 0.0) {
                        rowsA[gid_pxSxx(gi)][gid_xx(gj)] += hd_x * gxj;   // ∂x σxx
                        rowsA[gid_pxSxz(gi)][gid_xz(gj)] += hd_x * gxj;   // ∂x σxz
                        rowsA[gid_pxVx (gi)][gid_vx(gj)] += hd_x * gxj;   // ∂x vx
                        rowsA[gid_pxVz (gi)][gid_vz(gj)] += hd_x * gxj;   // ∂x vz
                    }
                    if (hd_z != 0.0) {
                        rowsA[gid_pzSxz(gi)][gid_xz(gj)] += hd_z * gzj;   // ∂z σxz
                        rowsA[gid_pzSzz(gi)][gid_zz(gj)] += hd_z * gzj;   // ∂z σzz
                        rowsA[gid_pzVz (gi)][gid_vz(gj)] += hd_z * gzj;   // ∂z vz
                        rowsA[gid_pzVx (gi)][gid_vx(gj)] += hd_z * gzj;   // ∂z vx
                    }
                }
            }
        }
    }
#endif  // USE_Q1_MESH

    // ----- 构造 CSR -----
    out.Ntot = Ntot;
    out.A.build_from_rows(Ntot, rowsA);
    out.M.build_from_rows(Ntot, rowsM);
    rowsA.clear(); rowsA.shrink_to_fit();
    rowsM.clear(); rowsM.shrink_to_fit();

    printf("[assemble] done.  nnz(A)=%d  nnz(M)=%d   (%.2f s)\n",
           out.A.nnz(), out.M.nnz(),
           (double)(clock() - tA0) / CLOCKS_PER_SEC);
}

// -----------------------------------------------------------------------------
// 组装本步 b
// -----------------------------------------------------------------------------
void assemble_rhs(const FemSystem &sys,
                  const femmesh::Mesh &mesh,
                  const std::vector<double> &q_n,
                  double dt,
                  double sou_t,
                  std::vector<double> &b_out)
{
    if ((int)b_out.size() != Ntot) b_out.assign(Ntot, 0.0);

    // (1) 主 5 场 RHS
    //
    //   Backward Euler : b = M·q^n
    //   Crank–Nicolson : b = (M - hc·K)·q^n = 2·M·q^n - A·q^n
    //                    （A = M + hc·K，hc = dt/2；与 FDM/ImplicitV2 一致）
    sys.M.multiply(q_n, b_out);

#ifndef USE_BACKWARD_EULER
    {
        static std::vector<double> Aq;
        if ((int)Aq.size() != Ntot) Aq.assign(Ntot, 0.0);
        sys.A.multiply(q_n, Aq);
        // 只对主 5 场块做 CN 修正；ψ 块仍走下方专用公式
        const int n_main = 5 * Nnodes;
        for (int i = 0; i < n_main; i++)
            b_out[i] = 2.0 * b_out[i] - Aq[i];
    }
#endif

    // (2) 声源：与 FDM/ImplicitV2 完全一致 —— 直接加 dt·sou_t 到 σxx / σzz
    //
    //   FDM:  b[id_xx(jsou,isou)] += dt * sou_t;   (identity mass → Δσ ≈ dt·f)
    //
    //   FEM Galerkin: M·σ̇ = f  →  Δσ = dt·f / M_ii，要让 Δσ ≈ dt·f 需 f = M_ii·sou_t
    //   即 b += dt·sou_t·M0_diag[k]（consistent mass 对角，≠ M_L）
    {
        // 默认：**集中源**，不做空间平滑（用户在 SESSION 5 明确要求恢复裸点源，
        // 因为 FDM 同参数没有此问题，说明散射的根源不在源，而在离散离散化）。
        //
        // 若需高频抑制诊断，编译期传：
        //   -DSOURCE_SMOOTH_RADIUS=3 -DSOURCE_SMOOTH_SIGMA=1.7   → 7×7 高斯
        // 强制关闭平滑（等价于默认）：
        //   -DNO_SOURCE_SMOOTH
#ifdef NO_SOURCE_SMOOTH
        const int    SR    = 0;
        const double SIGMA = 1.0;   // 未使用
#else
    #ifndef SOURCE_SMOOTH_RADIUS
        #define SOURCE_SMOOTH_RADIUS 0   // 默认：集中源（无平滑）
    #endif
    #ifndef SOURCE_SMOOTH_SIGMA
        #define SOURCE_SMOOTH_SIGMA 1.0
    #endif
        const int    SR    = SOURCE_SMOOTH_RADIUS;
        const double SIGMA = SOURCE_SMOOTH_SIGMA;
#endif
        // 先累计 Gaussian 未归一权重和真实有效节点，再归一化 → 域内 ∑w = 1
        const int MAX_SW = 15;                              // 上限
        double wsum = 0.0;
        int    nd_list[MAX_SW * MAX_SW];
        double w_list [MAX_SW * MAX_SW];
        int    NW = 0;

        const double two_sigma2 = (SR == 0) ? 1.0 : (2.0 * SIGMA * SIGMA);
        for (int dj = -SR; dj <= SR; dj++) {
            for (int di = -SR; di <= SR; di++) {
                const int ix = isou + di;
                const int iz = jsou + dj;
                if (ix < 0 || ix >= Nx || iz < 0 || iz >= Nz) continue;
                double w;
                if (SR == 0) {
                    w = 1.0;                                // 集中源
                } else {
                    w = std::exp(-(double)(di * di + dj * dj) / two_sigma2);
                }
                nd_list[NW] = ix + iz * Nx;
                w_list [NW] = w;
                wsum += w;
                NW++;
            }
        }
        // 归一化 + 施加
        const double inv_wsum = 1.0 / wsum;
        for (int k = 0; k < NW; k++) {
            const int    nd = nd_list[k];
            const double sc = dt * sou_t * (w_list[k] * inv_wsum) * mesh.M0_diag[nd];
            b_out[gid_xx(nd)] += sc;
            b_out[gid_zz(nd)] += sc;
        }
    }

    // (3) 8 个 ψ 行的额外修正：
    //     b[ψ_k] = E_α(x_k) · ( M_L·ψ^n  -  (dt/2)·d_α·(G·field^n)_k )
    //
    //     其中 (G·field^n)_k 通过一次 element loop 汇总。
    //     8 个 (field, direction) 组合都在同一循环内完成。
    //
    //     实现顺序：
    //       a) 累计 8 个 “G·field^n” 数组 grad_*[Nnodes]
    //       b) 对每个 ψ 块的每个节点做 b[ψ_k] := E · (b[ψ_k] - (dt/2)·d·grad_k)

    const double half_dt = 0.5 * dt;

    // 8 个 gradient 数组
    std::vector<double> g_xSxx(Nnodes, 0.0);   // ∂x σxx
    std::vector<double> g_zSxz(Nnodes, 0.0);   // ∂z σxz
    std::vector<double> g_xSxz(Nnodes, 0.0);   // ∂x σxz
    std::vector<double> g_zSzz(Nnodes, 0.0);   // ∂z σzz
    std::vector<double> g_xVx (Nnodes, 0.0);   // ∂x vx
    std::vector<double> g_zVz (Nnodes, 0.0);   // ∂z vz
    std::vector<double> g_xVz (Nnodes, 0.0);   // ∂x vz
    std::vector<double> g_zVx (Nnodes, 0.0);   // ∂z vx

#ifdef USE_Q1_MESH
    // ---- Q1 版本的 (G_α · field)_k 累加 ----
    //   g_α_field[i] = Σ_e Σ_j Dα[i][j]_e · field[j]   (i,j 走遍 element 4 节点)
    //   注意 Q1 中 ∂N/∂x 依赖 z 分量，所以 Dα[i][j] 依赖 i（与 P1 不同）。
    static const int mDx_r [4][4] = { {-2,2,1,-1},{-2,2,1,-1},{-1,1,2,-2},{-1,1,2,-2} };
    static const int mDz_r [4][4] = { {-2,-1,1,2},{-1,-2,2,1},{-1,-2,2,1},{-2,-1,1,2} };

    for (const auto &qe : mesh.quads) {
        const double kDx = qe.dz_val / 12.0;
        const double kDz = qe.dx_val / 12.0;
        const int gn[4] = { qe.nodes[0], qe.nodes[1], qe.nodes[2], qe.nodes[3] };
        for (int i = 0; i < 4; i++) {
            const int ni = gn[i];
            double axx=0, azz=0, axz_x=0, axz_z=0, avx_x=0, avx_z=0, avz_x=0, avz_z=0;
            for (int j = 0; j < 4; j++) {
                const double Dxij = kDx * (double)mDx_r[i][j];
                const double Dzij = kDz * (double)mDz_r[i][j];
                const int nj = gn[j];
                axx    += Dxij * q_n[off_xx + nj];
                azz    += Dzij * q_n[off_zz + nj];
                axz_x  += Dxij * q_n[off_xz + nj];   // ∂x σxz
                axz_z  += Dzij * q_n[off_xz + nj];   // ∂z σxz
                avx_x  += Dxij * q_n[off_vx + nj];   // ∂x vx
                avx_z  += Dzij * q_n[off_vx + nj];   // ∂z vx
                avz_x  += Dxij * q_n[off_vz + nj];   // ∂x vz
                avz_z  += Dzij * q_n[off_vz + nj];   // ∂z vz
            }
            g_xSxx[ni] += axx;
            g_zSzz[ni] += azz;
            g_xSxz[ni] += axz_x;
            g_zSxz[ni] += axz_z;
            g_xVx [ni] += avx_x;
            g_zVx [ni] += avx_z;
            g_xVz [ni] += avz_x;
            g_zVz [ni] += avz_z;
        }
    }
#else
    for (const auto &t : mesh.triangles) {
        const double A3 = t.A / 3.0;
        const int n0 = t.nodes[0], n1 = t.nodes[1], n2 = t.nodes[2];

        // 每个 field 在这个 element 上的常数梯度
        auto grad_of = [&](int off, double dNa[3]) -> double {
            return dNa[0] * q_n[off + n0]
                 + dNa[1] * q_n[off + n1]
                 + dNa[2] * q_n[off + n2];
        };

        double dNdx[3] = { t.dN_dx[0], t.dN_dx[1], t.dN_dx[2] };
        double dNdz[3] = { t.dN_dz[0], t.dN_dz[1], t.dN_dz[2] };

        const double dxSxx = grad_of(off_xx, dNdx);   // ∂x σxx
        const double dzSxz = grad_of(off_xz, dNdz);   // ∂z σxz
        const double dxSxz = grad_of(off_xz, dNdx);   // ∂x σxz
        const double dzSzz = grad_of(off_zz, dNdz);   // ∂z σzz
        const double dxVx  = grad_of(off_vx, dNdx);   // ∂x vx
        const double dzVz  = grad_of(off_vz, dNdz);   // ∂z vz
        const double dxVz  = grad_of(off_vz, dNdx);   // ∂x vz
        const double dzVx  = grad_of(off_vx, dNdz);   // ∂z vx

        // 分散 (|A|/3)·(∂α f) 到 3 个节点
        const int ns[3] = { n0, n1, n2 };
        for (int a = 0; a < 3; a++) {
            const int nk = ns[a];
            g_xSxx[nk] += A3 * dxSxx;
            g_zSxz[nk] += A3 * dzSxz;
            g_xSxz[nk] += A3 * dxSxz;
            g_zSzz[nk] += A3 * dzSzz;
            g_xVx [nk] += A3 * dxVx;
            g_zVz [nk] += A3 * dzVz;
            g_xVz [nk] += A3 * dxVz;
            g_zVx [nk] += A3 * dzVx;
        }
    }
#endif  // USE_Q1_MESH

    // 对每个节点应用 E·(...)
    const std::vector<double> &Ex = sys.Ex_at_ix;
    const std::vector<double> &Ez = sys.Ez_at_iz;
    const std::vector<double> &dx_ = sys.dampx_at_ix;
    const std::vector<double> &dz_ = sys.dampz_at_iz;

    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const int k    = ix + iz * Nx;
            const double Exk = Ex[ix];
            const double Ezk = Ez[iz];
            const double dxk = dx_[ix];
            const double dzk = dz_[iz];

            // b[ψ_k] := E · ( b[ψ_k] - (dt/2)·d·grad_k )
            //   b[ψ_k] 现在等于 M_L_k · ψ^n_k
            // 8 个块统一处理：
            b_out[gid_pxSxx(k)] = Exk * (b_out[gid_pxSxx(k)] - half_dt * dxk * g_xSxx[k]);
            b_out[gid_pzSxz(k)] = Ezk * (b_out[gid_pzSxz(k)] - half_dt * dzk * g_zSxz[k]);
            b_out[gid_pxSxz(k)] = Exk * (b_out[gid_pxSxz(k)] - half_dt * dxk * g_xSxz[k]);
            b_out[gid_pzSzz(k)] = Ezk * (b_out[gid_pzSzz(k)] - half_dt * dzk * g_zSzz[k]);
            b_out[gid_pxVx (k)] = Exk * (b_out[gid_pxVx (k)] - half_dt * dxk * g_xVx [k]);
            b_out[gid_pzVz (k)] = Ezk * (b_out[gid_pzVz (k)] - half_dt * dzk * g_zVz [k]);
            b_out[gid_pxVz (k)] = Exk * (b_out[gid_pxVz (k)] - half_dt * dxk * g_xVz [k]);
            b_out[gid_pzVx (k)] = Ezk * (b_out[gid_pzVx (k)] - half_dt * dzk * g_zVx [k]);
        }
    }
}

} // namespace femasm
