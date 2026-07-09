/* =============================================================================
 *  assembler.cpp   —  位移二阶格式装配（方案 B / SESSION 7）
 *  ---------------------------------------------------------------------------
 *  组装 K（标准刚度阵）、M（consistent mass）、C（吸收层质量比例阻尼），
 *  并合成 Newmark 有效矩阵  A = M + γ·dt·C + β·dt²·K。
 *
 *  为什么这样能根治菱形：K = ∫ Bᵀ D B 里梯度以 ∫∇N·∇N 出现，是连中心节点的
 *  正规二阶差分（不像 velocity-stress collocated 那样跳中心 → 2Δx 各向异性）。
 * ===========================================================================*/
#include "assembler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <map>

namespace femasm {

using femcfg::Nx;
using femcfg::Nz;
using femcfg::Nnodes;
using femcfg::Ntot;
using femcfg::off_ux;
using femcfg::off_uz;
using femcfg::isou;
using femcfg::jsou;
using femcfg::dx;
using femcfg::dz;
using femcfg::NABS;

using femla::RowMap;

// 局部 dof (顶点 a ∈ {0,1,2}, 分量 comp ∈ {0=ux,1=uz}) → 全局 dof
static inline int gdof(int gnode, int comp)
{
    return (comp == 0 ? off_ux : off_uz) + gnode;
}

// -----------------------------------------------------------------------------
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  ElastSystem &out)
{
    printf("[assemble] Ntot=%d  Ntri=%d  (P1 tri, displacement 2nd-order, Newmark)\n",
           Ntot, mesh.n_tris());
    fflush(stdout);
    clock_t tA0 = clock();

    const double beta  = out.beta;
    const double gamma = out.gamma;

    // ---- PML 方向性剖面（节点级）----
    out.dampx_ix.assign(Nx, 0.0); out.dampz_iz.assign(Nz, 0.0);
    out.Ex_ix.assign(Nx, 1.0);    out.Ez_iz.assign(Nz, 1.0);
    out.facx_ix.assign(Nx, dt);   out.facz_iz.assign(Nz, dt);
    double dmax = 0.0; int n_pml = 0;
    for (int ix = 0; ix < Nx; ix++) {
        const double d = femcfg::dampx_at(ix * dx);
        out.dampx_ix[ix] = d;
        if (d > 0.0) { out.Ex_ix[ix] = std::exp(-d * dt);
                       out.facx_ix[ix] = (1.0 - out.Ex_ix[ix]) / d; }
        if (d > dmax) dmax = d;
    }
    for (int iz = 0; iz < Nz; iz++) {
        const double d = femcfg::dampz_at(iz * dz);
        out.dampz_iz[iz] = d;
        if (d > 0.0) { out.Ez_iz[iz] = std::exp(-d * dt);
                       out.facz_iz[iz] = (1.0 - out.Ez_iz[iz]) / d; }
        if (d > dmax) dmax = d;
    }
    for (int iz = 0; iz < Nz; iz++)
        for (int ix = 0; ix < Nx; ix++)
            if (out.dampx_ix[ix] > 0.0 || out.dampz_iz[iz] > 0.0) n_pml++;
    out.pml_nodes = n_pml;
    out.damp_max  = dmax;

    std::vector<RowMap> rowsA(Ntot);
    std::vector<RowMap> rowsKeff(Ntot);   // 标准 K + 角区弹簧 P（进 A、组 b）
    std::vector<RowMap> rowsC(Ntot);      // PML 阻尼
    std::vector<RowMap> rowsKxx(Ntot);    // x 向刚度块（驱动 ψ_xx）
    std::vector<RowMap> rowsKzz(Ntot);    // z 向刚度块（驱动 ψ_zz）

    // ---- 遍历三角形装配 ----
    for (const auto &t : mesh.triangles) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double rho = mp.rho;
        const double lam = mp.lambda;
        const double mu  = mp.mu;
        const double Cc  = lam + 2.0 * mu;               // λ+2μ

        const double A_e = t.A;
        const double dNx[3] = { t.dN_dx[0], t.dN_dx[1], t.dN_dx[2] };
        const double dNz[3] = { t.dN_dz[0], t.dN_dz[1], t.dN_dz[2] };
        const int    gn [3] = { t.nodes[0], t.nodes[1], t.nodes[2] };

        // 单元级 PML 系数（节点平均）：
        //   d_sum = <d_x + d_z>（用于 C = ∫ρ(d_x+d_z)NN）
        //   d_pro = <d_x·d_z>  （用于 P = ∫ρ d_x d_z NN，仅角区非零）
        double d_sum = 0.0, d_pro = 0.0;
        if (NABS > 0) {
            for (int a = 0; a < 3; a++) {
                const double dxk = out.dampx_ix[gn[a] % Nx];
                const double dzk = out.dampz_iz[gn[a] / Nx];
                d_sum += (dxk + dzk);
                d_pro += (dxk * dzk);
            }
            d_sum /= 3.0; d_pro /= 3.0;
        }

        // B 矩阵 (3×6)
        double B[3][6] = {{0}};
        for (int a = 0; a < 3; a++) {
            B[0][2*a+0] = dNx[a];
            B[1][2*a+1] = dNz[a];
            B[2][2*a+0] = dNz[a];
            B[2][2*a+1] = dNx[a];
        }
        const double D[3][3] = { { Cc,  lam, 0.0 },
                                 { lam, Cc,  0.0 },
                                 { 0.0, 0.0, mu  } };
        double DB[3][6];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 6; c++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += D[r][k] * B[k][c];
                DB[r][c] = s;
            }
        double Ke[6][6];                                 // 完整单元刚度
        for (int p = 0; p < 6; p++)
            for (int q = 0; q < 6; q++) {
                double s = 0.0;
                for (int r = 0; r < 3; r++) s += B[r][p] * DB[r][q];
                Ke[p][q] = A_e * s;
            }

        // ---- 按【检验函数导数方向】拆分刚度：Ke = Ke_dx + Ke_dz ----
        // 内力 ∫(∇w):σ = ∫(∂x w)·σ_·x + ∫(∂z w)·σ_·z。PML 复坐标拉伸对
        // ∂x 项乘 1/s_x、∂z 项乘 1/s_z → 需把内力按检验函数的 ∂x / ∂z 分开。
        //   检验侧 B 行：row0(εxx)=∂x w_ux, row1(εzz)=∂z w_uz, row2(γxz)=∂z w_ux+∂x w_uz
        //   ∂x 部分：row0 全部 + row2 的 w_uz 分量；∂z 部分：row1 全部 + row2 的 w_ux 分量
        // 逐单元本地 dof p=2a+comp（comp0=ux,comp1=uz）。
        // Bx[r][p]/Bz[r][p]：把 B 的检验侧按导数方向劈开（Bx+Bz=B）。
        double Bx[3][6] = {{0}}, Bz[3][6] = {{0}};
        for (int a = 0; a < 3; a++) {
            // ux 检验(comp0)：εxx 用 ∂x，γxz 用 ∂z
            Bx[0][2*a+0] = dNx[a];      // ∂x w_ux → εxx 行
            Bz[2][2*a+0] = dNz[a];      // ∂z w_ux → γxz 行
            // uz 检验(comp1)：εzz 用 ∂z，γxz 用 ∂x
            Bz[1][2*a+1] = dNz[a];      // ∂z w_uz → εzz 行
            Bx[2][2*a+1] = dNx[a];      // ∂x w_uz → γxz 行
        }
        double Ke_dx[6][6], Ke_dz[6][6];   // 检验-∂x / 检验-∂z 分解（Ke_dx+Ke_dz=Ke）
        for (int p = 0; p < 6; p++)
            for (int q = 0; q < 6; q++) {
                double sx = 0.0, sz = 0.0;
                for (int r = 0; r < 3; r++) { sx += Bx[r][p]*DB[r][q]; sz += Bz[r][p]*DB[r][q]; }
                Ke_dx[p][q] = A_e * sx;
                Ke_dz[p][q] = A_e * sz;
            }

        // 混合质量矩阵（抗 P1 频散）：M = (1-α)·consistent + α·lumped
        //   consistent: diag=A_e/6, off=A_e/12   lumped: diag=A_e/3, off=0
        //   α=0 纯 consistent（默认历史行为）；α=1 纯 lumped。
        //   两者 P1 频散误差符号相反，适当 α 抵消主导误差、压制频散尾。
        //   仍对称 SPD（lumped 对角正、consistent SPD）→ A 仍对称 SPD，不影响 CIM。
        const double A12  = A_e / 12.0;
        const double mblend = femcfg::mass_blend;
        const double m_diag_c = A_e / 6.0, m_off_c = A_e / 12.0;   // consistent
        const double m_diag_l = A_e / 3.0;                         // lumped diag
        const double m_diag_blend = (1.0 - mblend) * m_diag_c + mblend * m_diag_l;
        const double m_off_blend  = (1.0 - mblend) * m_off_c;
        (void)A12;

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                const double m0 = (a == b) ? m_diag_blend : m_off_blend;
                const double me = rho * m0;              // 质量
                const double ce = d_sum * me;            // PML 阻尼 C
                const double pe = d_pro * me;            // PML 弹簧 P

                // 分量对角块（ux-ux, uz-uz）
                for (int comp = 0; comp < 2; comp++) {
                    const int gi = gdof(gn[a], comp);
                    const int gj = gdof(gn[b], comp);
                    const double kij = Ke[2*a+comp][2*b+comp];   // 完整刚度该分量项
                    const double keff = kij + pe;                // K_eff = K + P
                    rowsKeff[gi][gj] += keff;
                    rowsC   [gi][gj] += ce;
                    rowsA   [gi][gj] += me + gamma * dt * ce + beta * dt * dt * keff;

                    // 方向分解内力 K_dirx / K_dirz（驱动 ψ，只算 RHS，不进 A）
                    // 对角分量部分。K_dirx+K_dirz = K（完整刚度）→ PML 外精确退化。
                    if (NABS > 0) {
                        rowsKxx[gi][gj] += Ke_dx[2*a+comp][2*b+comp];
                        rowsKzz[gi][gj] += Ke_dz[2*a+comp][2*b+comp];
                    }
                }
                // 分量交叉项 (ux-uz, uz-ux)：K_eff 里 + 方向分解内力里都要有
                {
                    const int gi_x = gdof(gn[a], 0);
                    const int gj_z = gdof(gn[b], 1);
                    const double kxz = Ke[2*a+0][2*b+1];
                    rowsKeff[gi_x][gj_z] += kxz;
                    rowsA   [gi_x][gj_z] += beta * dt * dt * kxz;

                    const int gi_z = gdof(gn[a], 1);
                    const int gj_x = gdof(gn[b], 0);
                    const double kzx = Ke[2*a+1][2*b+0];
                    rowsKeff[gi_z][gj_x] += kzx;
                    rowsA   [gi_z][gj_x] += beta * dt * dt * kzx;

                    if (NABS > 0) {
                        rowsKxx[gi_x][gj_z] += Ke_dx[2*a+0][2*b+1];
                        rowsKzz[gi_x][gj_z] += Ke_dz[2*a+0][2*b+1];
                        rowsKxx[gi_z][gj_x] += Ke_dx[2*a+1][2*b+0];
                        rowsKzz[gi_z][gj_x] += Ke_dz[2*a+1][2*b+0];
                    }
                }
            }
        }
    }

    // ---- 声源等效节点力系数：F_a = M0·∇N_a ----
    //
    // 【SESSION 12】P1 声源空间平滑：单节点点源本质宽带（含 λ<2dx 分量），P1
    // 网格无法干净传播这些超短波长 → 拖出后期尾波频散。修法：把矩张量源在源
    // 周围节点上做 **3×3 高斯平滑**（展宽点源、削掉 λ<2dx 分量）。对每个"子源
    // 节点"仍施加标准 F_a=M0·∇N_a，按归一化高斯权重（∑w=1，保持总矩不变）叠加。
    //   -DSRC_SMOOTH_R=R    平滑半径（节点数，窗口 (2R+1)²）。默认 1（3×3）。
    //   -DSRC_SMOOTH_SIGMA=σ  高斯 σ（节点单位）。默认 1.0。
    //   R=0 退回原始集中点源。
    // 只改 RHS 的声源力系数，不动 A（A 仍对称 SPD，CIM 不受影响）。
#ifndef SRC_SMOOTH_R
    #define SRC_SMOOTH_R 1
#endif
#ifndef SRC_SMOOTH_SIGMA
    #define SRC_SMOOTH_SIGMA 1.0
#endif
    const int    SR   = SRC_SMOOTH_R;
    const double SSIG = SRC_SMOOTH_SIGMA;
    std::map<int, double> srcmap;

    // 收集子源节点及归一化高斯权重
    std::vector<int>    svnode;
    std::vector<double> sw;
    double wsum = 0.0;
    for (int dj = -SR; dj <= SR; dj++)
        for (int di = -SR; di <= SR; di++) {
            const int ix = isou + di, iz = jsou + dj;
            if (ix < 1 || ix >= Nx-1 || iz < 1 || iz >= Nz-1) continue;
            const double w = (SR == 0) ? 1.0
                           : std::exp(-(double)(di*di+dj*dj) / (2.0*SSIG*SSIG));
            svnode.push_back(ix + iz * Nx); sw.push_back(w); wsum += w;
        }
    for (double &w : sw) w /= wsum;

    // 对每个子源节点施加标准矩张量等效力 F_a = w·∇N_a
    for (size_t s = 0; s < svnode.size(); s++) {
        const int src_node = svnode[s];
        const double wgt = sw[s];
        for (const auto &t : mesh.triangles) {
            int loc = -1;
            for (int a = 0; a < 3; a++) if (t.nodes[a] == src_node) loc = a;
            if (loc < 0) continue;
            for (int a = 0; a < 3; a++) {
                srcmap[gdof(t.nodes[a], 0)] += wgt * t.dN_dx[a];
                srcmap[gdof(t.nodes[a], 1)] += wgt * t.dN_dz[a];
            }
        }
    }
    out.src_dof.clear();
    out.src_coef.clear();
    for (const auto &kv : srcmap) {
        out.src_dof .push_back(kv.first);
        out.src_coef.push_back(kv.second);
    }

    // ---- 构造 CSR ----
    out.Ntot = Ntot;
    out.A.build_from_rows(Ntot, rowsA);
    out.K_eff.build_from_rows(Ntot, rowsKeff);
    out.C.build_from_rows(Ntot, rowsC);
    out.K_xx.build_from_rows(Ntot, rowsKxx);
    out.K_zz.build_from_rows(Ntot, rowsKzz);
    rowsA.clear(); rowsA.shrink_to_fit();
    rowsKeff.clear(); rowsKeff.shrink_to_fit();
    rowsC.clear(); rowsC.shrink_to_fit();
    rowsKxx.clear(); rowsKxx.shrink_to_fit();
    rowsKzz.clear(); rowsKzz.shrink_to_fit();

    printf("[assemble] nnz(A)=%d nnz(Keff)=%d nnz(C)=%d nnz(Kxx)=%d nnz(Kzz)=%d src=%d  (%.2f s)\n",
           out.A.nnz(), out.K_eff.nnz(), out.C.nnz(), out.K_xx.nnz(), out.K_zz.nnz(),
           (int)out.src_dof.size(), (double)(clock() - tA0) / CLOCKS_PER_SEC);
    if (NABS > 0)
        printf("[pml] ADE-PML unsplit  NABS=%d  d0_x=%g d0_z=%g  pml_nodes=%d  max(d)=%g\n",
               NABS, femcfg::compute_d0_x(), femcfg::compute_d0_z(), n_pml, dmax);
    else
        printf("[pml] NABS=0 -> no PML (degenerates to pure displacement, A unchanged)\n");
}

// -----------------------------------------------------------------------------
// PML 记忆变量显式更新 + F_ψ = ψ_x + ψ_z   （SESSION 13：正确的 unsplit CFS）
//
//   复坐标拉伸：内力里 ∂x 项 ×(1/s_x)、∂z 项 ×(1/s_z)，s_α=1+d_α/(iω)。
//   1/s_α = 1 - d_α/(iω+d_α)。设方向分解内力 g_x=K_dirx·u, g_z=K_dirz·u
//   （g_x+g_z=K·u，已在 K_eff 内）。拉伸后内力 = K·u - (ψ_x + ψ_z)，其中
//        ψ̇_x + d_x ψ_x = d_x g_x      （x 方向记忆，只用 d_x 驱动 g_x）
//        ψ̇_z + d_z ψ_z = d_z g_z      （z 方向记忆，只用 d_z 驱动 g_z）
//   —— 关键：ψ_x、ψ_z 各自独立、**不做 (dz-dx) 差**，消除旧版抵消 bug。
//   指数积分（逐 DOF）：ψ^{n+1}=E·ψ^n + (1-E)·g   （E=exp(-d·dt)，稳定精确）
//        因为稳态 ψ→g，(1-E) 使离散稳态也 →g，正确。
//   d=0 处 E=1 → ψ^{n+1}=ψ^n，若初值 0 则恒 0，F_ψ=0，精确退化回无 PML。
//   F_ψ = ψ_x + ψ_z，在 main 里作为内力修正从 RHS 减去（b -= F_ψ）。不进 A。
// -----------------------------------------------------------------------------
void pml_update_force(const ElastSystem &sys,
                      const std::vector<double> &u_pred,
                      std::vector<double> &psi_x,
                      std::vector<double> &psi_z,
                      std::vector<double> &Fpsi_out)
{
    if ((int)Fpsi_out.size() != Ntot) Fpsi_out.assign(Ntot, 0.0);
    if (femcfg::NABS <= 0) { std::fill(Fpsi_out.begin(), Fpsi_out.end(), 0.0); return; }
    if ((int)psi_x.size() != Ntot) psi_x.assign(Ntot, 0.0);
    if ((int)psi_z.size() != Ntot) psi_z.assign(Ntot, 0.0);

    static std::vector<double> gx, gz;
    sys.K_xx.multiply(u_pred, gx);   // g_x = K_dirx · u（x 方向内力）
    sys.K_zz.multiply(u_pred, gz);   // g_z = K_dirz · u（z 方向内力）

    const int MX = Nx;
    for (int i = 0; i < Ntot; i++) {
        const int node = i % Nnodes;
        const int ix = node % MX;
        const int iz = node / MX;
        const double Ex = sys.Ex_ix[ix];
        const double Ez = sys.Ez_iz[iz];
        // ψ_x^{n+1} = E_x·ψ_x^n + (1-E_x)·g_x   （d=0 → E=1，ψ 不变；恒 0）
        psi_x[i] = Ex * psi_x[i] + (1.0 - Ex) * gx[i];
        psi_z[i] = Ez * psi_z[i] + (1.0 - Ez) * gz[i];
        Fpsi_out[i] = psi_x[i] + psi_z[i];
    }
}

// -----------------------------------------------------------------------------
void assemble_source_force(const ElastSystem &sys,
                           double sou_t,
                           std::vector<double> &F_out)
{
    if ((int)F_out.size() != Ntot) F_out.assign(Ntot, 0.0);
    else std::fill(F_out.begin(), F_out.end(), 0.0);

    const double M0 = femcfg::src_scale * sou_t;
    for (size_t k = 0; k < sys.src_dof.size(); k++)
        F_out[sys.src_dof[k]] += M0 * sys.src_coef[k];
}

// -----------------------------------------------------------------------------
void recover_nodal_stress(const femmesh::Mesh &mesh,
                          const femmat::Material &mat,
                          const std::vector<double> &u,
                          std::vector<double> &sxx,
                          std::vector<double> &szz,
                          std::vector<double> &sxz)
{
    sxx.assign(Nnodes, 0.0);
    szz.assign(Nnodes, 0.0);
    sxz.assign(Nnodes, 0.0);
    std::vector<double> wsum(Nnodes, 0.0);

    for (const auto &t : mesh.triangles) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double lam = mp.lambda, mu = mp.mu, Cc = lam + 2.0 * mu;
        const double A_e = t.A;
        const int gn[3] = { t.nodes[0], t.nodes[1], t.nodes[2] };

        double dux_dx = 0, dux_dz = 0, duz_dx = 0, duz_dz = 0;
        for (int a = 0; a < 3; a++) {
            const double ux = u[off_ux + gn[a]];
            const double uz = u[off_uz + gn[a]];
            dux_dx += t.dN_dx[a] * ux;
            dux_dz += t.dN_dz[a] * ux;
            duz_dx += t.dN_dx[a] * uz;
            duz_dz += t.dN_dz[a] * uz;
        }
        const double exx = dux_dx, ezz = duz_dz, gxz = dux_dz + duz_dx;
        const double s_xx = Cc * exx + lam * ezz;
        const double s_zz = lam * exx + Cc * ezz;
        const double s_xz = mu * gxz;

        for (int a = 0; a < 3; a++) {
            const int nd = gn[a];
            sxx[nd] += A_e * s_xx;
            szz[nd] += A_e * s_zz;
            sxz[nd] += A_e * s_xz;
            wsum[nd] += A_e;
        }
    }
    for (int k = 0; k < Nnodes; k++) {
        if (wsum[k] > 0.0) {
            const double inv = 1.0 / wsum[k];
            sxx[k] *= inv; szz[k] *= inv; sxz[k] *= inv;
        }
    }
}

} // namespace femasm
