/* =============================================================================
 *  assembler.h   —  位移二阶格式装配层（方案 B / SESSION 7）
 *  ---------------------------------------------------------------------------
 *  从「位移二阶弹性波弱形式 + P1 三角单元」得到
 *
 *      M ü + C u̇ + K u = f
 *
 *  以及 Newmark 平均加速度的有效系统
 *
 *      A a^{n+1} = b,     A = M + γ·dt·C + β·dt²·K   （时间无关，装配一次）
 *
 *  未知向量 2 * Nnodes 个自由度：块布局 [ux; uz]（见 config.h）。
 *
 *  单元量（P1 三角，B 与 D 在单元内常数）：
 *      ε = [εxx, εzz, γxz]ᵀ = B u_e,   B(3×6)，每节点列 [dNx,0; 0,dNz; dNz,dNx]
 *      D = [[λ+2μ, λ,    0 ],
 *           [λ,    λ+2μ, 0 ],
 *           [0,    0,    μ ]]                （平面应变各向同性）
 *      Ke = |A_e| · Bᵀ D B                   （6×6 单元刚度）
 *      Me = ρ_e · ∫ N_i N_j                  （consistent，块对角 ux/uz 各一份）
 *      Ce = d_e · Me                         （质量比例阻尼，d_e = 单元平均阻尼）
 *
 *  声源：各向同性矩张量（爆炸源），等效节点力 F_a = M0(t)·∇N_a，
 *      对源节点邻接的所有三角求和（自平衡、圆对称）。
 * ===========================================================================*/
#ifndef FEM_ASSEMBLER_H
#define FEM_ASSEMBLER_H

#include <vector>

#include "config.h"
#include "material.h"
#include "mesh.h"
#include "sparse_matrix.h"

namespace femasm {

struct ElastSystem {
    // ---- 进 CIM 的核心：对称正定、常系数、装配一次 ----
    //   A = M + γ·dt·C + β·dt²·K_eff
    //   M、C、K_eff 全对称 → A 对称；M SPD + 其余 PSD → A SPD。
    //   CIM/QUBO 后端把 A x = b 映射成 min ½xᵀAx − bᵀx（A SPD 时凸）即可。
    femla::SparseMatrixCSR A;

    // ---- 每步 RHS 需要的矩阵（不进 CIM，只做 SpMV 组 b） ----
    femla::SparseMatrixCSR K_eff;  // 标准刚度 K + PML 角区弹簧 P（对称）
    femla::SparseMatrixCSR C;      // PML 阻尼 ∫ρ(d_x+d_z)NN（对称，仅 PML 区非零）
    femla::SparseMatrixCSR K_xx;   // x 方向刚度块（驱动 ψ_xx，不进 A）
    femla::SparseMatrixCSR K_zz;   // z 方向刚度块（驱动 ψ_zz，不进 A）
    int Ntot = 0;

    double beta = femcfg::nm_beta;
    double gamma = femcfg::nm_gamma;
    double dt = femcfg::dt;

    // ---- PML 节点级剖面（方向性） ----
    std::vector<double> dampx_ix;  // d_x(ix·dx)，size Nx
    std::vector<double> dampz_iz;  // d_z(iz·dz)，size Nz
    std::vector<double> Ex_ix;     // exp(-d_x·dt)，size Nx
    std::vector<double> Ez_iz;     // exp(-d_z·dt)，size Nz
    std::vector<double> facx_ix;   // (1-Ex)/d_x（d→0 取 dt），size Nx
    std::vector<double> facz_iz;   // (1-Ez)/d_z（d→0 取 dt），size Nz
    int    pml_nodes = 0;
    double damp_max = 0.0;

    // 声源等效节点力系数（时间无关；每步乘 M0(t)）
    std::vector<int>    src_dof;
    std::vector<double> src_coef;
};

// 装配 A / K_eff / C / K_xx / K_zz + PML 剖面 + 声源力系数。
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  ElastSystem &out);

// PML 记忆变量显式更新 + 生成 F_ψ（加到内力，即 RHS 减去它）。
//   输入 u_pred（Newmark 预测位移，已知），in/out ψ_xx、ψ_zz（各 Ntot）。
//   算法：g = K_xx·u_pred / K_zz·u_pred（SpMV），逐节点指数积分更新 ψ，
//         F_ψ = ψ_xx + ψ_zz。全部显式，不触碰 A。
void pml_update_force(const ElastSystem &sys,
                      const std::vector<double> &u_pred,
                      std::vector<double> &psi_xx,
                      std::vector<double> &psi_zz,
                      std::vector<double> &Fpsi_out);

// 组装本步声源等效力向量 F（大小 Ntot）：F = M0(t) · Σ ∇N_a。
void assemble_source_force(const ElastSystem &sys,
                           double sou_t,
                           std::vector<double> &F_out);

// 从位移场恢复节点应力（面积加权投影），用于与 FDM 快照对照。
//   sxx, szz, sxz 各大小 Nnodes。
void recover_nodal_stress(const femmesh::Mesh &mesh,
                          const femmat::Material &mat,
                          const std::vector<double> &u,
                          std::vector<double> &sxx,
                          std::vector<double> &szz,
                          std::vector<double> &sxz);

} // namespace femasm

#endif // FEM_ASSEMBLER_H
