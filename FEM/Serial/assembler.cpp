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

        const double A12 = A_e / 12.0;                   // consistent mass 模板基数

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                const double m0 = A12 * (a == b ? 2.0 : 1.0);
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

                    // K_xx / K_zz 分量对角部分（驱动 ψ，不进 A）
                    if (NABS > 0) {
                        const double coef_xx = (comp == 0) ? Cc : mu;   // ux用λ+2μ, uz用μ
                        const double coef_zz = (comp == 0) ? mu : Cc;
                        rowsKxx[gi][gj] += A_e * coef_xx * dNx[a] * dNx[b];
                        rowsKzz[gi][gj] += A_e * coef_zz * dNz[a] * dNz[b];
                    }
                }
                // 分量交叉项 (ux-uz, uz-ux)：只在 K_eff 里（M、C、P 无交叉）
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
                }
            }
        }
    }

    // ---- 声源等效节点力系数：F_a = M0·∇N_a，对源邻接三角求和 ----
    const int src_node = isou + jsou * Nx;
    std::map<int, double> srcmap;
    for (const auto &t : mesh.triangles) {
        int loc = -1;
        for (int a = 0; a < 3; a++) if (t.nodes[a] == src_node) loc = a;
        if (loc < 0) continue;
        for (int a = 0; a < 3; a++) {
            srcmap[gdof(t.nodes[a], 0)] += t.dN_dx[a];
            srcmap[gdof(t.nodes[a], 1)] += t.dN_dz[a];
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
// PML 记忆变量显式更新 + F_ψ = ψ_xx + ψ_zz
//   ODE:  ψ̇_xx + d_x·ψ_xx = (d_z - d_x)·(K_xx·u)
//         ψ̇_zz + d_z·ψ_zz = (d_x - d_z)·(K_zz·u)
//   指数积分（逐 DOF，d 取该 DOF 所在节点的方向阻尼）：
//         ψ^{n+1} = E·ψ^n + fac·rhs,  E=exp(-d·dt), fac=(1-E)/d (d→0 取 dt)
//   驱动量 K_xx·u_pred / K_zz·u_pred 用 Newmark 预测位移（已知，保持全显式）。
// -----------------------------------------------------------------------------
void pml_update_force(const ElastSystem &sys,
                      const std::vector<double> &u_pred,
                      std::vector<double> &psi_xx,
                      std::vector<double> &psi_zz,
                      std::vector<double> &Fpsi_out)
{
    if ((int)Fpsi_out.size() != Ntot) Fpsi_out.assign(Ntot, 0.0);
    if (femcfg::NABS <= 0) { std::fill(Fpsi_out.begin(), Fpsi_out.end(), 0.0); return; }
    if ((int)psi_xx.size() != Ntot) psi_xx.assign(Ntot, 0.0);
    if ((int)psi_zz.size() != Ntot) psi_zz.assign(Ntot, 0.0);

    static std::vector<double> gxx, gzz;
    sys.K_xx.multiply(u_pred, gxx);
    sys.K_zz.multiply(u_pred, gzz);

    // DOF i 属于块 (comp,node)，node=(i%Nnodes)，其 (ix,iz) 决定 d_x,d_z
    for (int i = 0; i < Ntot; i++) {
        const int node = i % Nnodes;
        const int ix = node % Nx;
        const int iz = node / Nx;
        const double dx_i = sys.dampx_ix[ix];
        const double dz_i = sys.dampz_iz[iz];

        psi_xx[i] = sys.Ex_ix[ix] * psi_xx[i]
                  + sys.facx_ix[ix] * (dz_i - dx_i) * gxx[i];
        psi_zz[i] = sys.Ez_iz[iz] * psi_zz[i]
                  + sys.facz_iz[iz] * (dx_i - dz_i) * gzz[i];
        Fpsi_out[i] = psi_xx[i] + psi_zz[i];
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
