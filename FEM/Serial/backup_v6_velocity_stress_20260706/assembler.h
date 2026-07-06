/* =============================================================================
 *  assembler.h
 *  ---------------------------------------------------------------------------
 *  从「PDE 弱形式 + P1 三角单元 + ADE-PML」得到全局
 *
 *      A x = b
 *
 *  的组装层（README §3.2 的“离散层”）。
 *
 *  未知向量共 13 * Nnodes 个自由度（详见 config.h 的块偏移表）：
 *      5 主场：vx, vz, σxx, σzz, σxz
 *      8 ψ  ：ψ_xSxx, ψ_zSxz, ψ_xSxz, ψ_zSzz,
 *              ψ_xVx , ψ_zVz , ψ_xVz , ψ_zVx
 *
 *  时间格式：
 *      - 主 5 方程：Crank–Nicolson（默认，与 FDM/ImplicitV2 一致）
 *                   A = M + (dt/2)K，b = (M - dt/2·K)q^n + dt·f
 *      - 可选：-DUSE_BACKWARD_EULER 退回 A = M + dt·K，b = M·q^n + dt·f
 *      - 8 个 ψ  ：指数 + 梯形 (与 FDM/ImplicitV2 一致)。
 *          ψ_α^{n+1} = E_α · ψ_α^n + (-d_α·dt/2)·[E_α·g^n + g^{n+1}]
 *          → M_L·ψ_α^{n+1} + (d·dt/2)·G·field^{n+1}
 *              = E·[M_L·ψ_α^n - (d·dt/2)·G·field^n]
 *      其中 M_L 是 lumped mass 对角矩阵，保证 d=0 处 ψ ≡ ψ^0。
 *
 *  关键块结构：
 *    Row vx (主 5 场耦合已有，加 PML 项)：
 *        A[vx , vx ] += Mρ
 *        A[vx , xx ] += dt·Gx           (Gx_ij = ∫ ∂N_i/∂x N_j)
 *        A[vx , xz ] += dt·Gz
 *        A[vx , ψ_xSxx] += -dt·M0       (ψ 引入 stretched-coord 修正)
 *        A[vx , ψ_zSxz] += -dt·M0
 *
 *    Row σxx：
 *        A[xx , xx ] += M0
 *        A[xx , vx ] += -dt·DxC         (DxC_ij = ∫ N_i C ∂N_j/∂x)
 *        A[xx , vz ] += -dt·DzL
 *        A[xx , ψ_xVx] += -dt·M_C       (M_C_ij = C_e * M0^e_ij)
 *        A[xx , ψ_zVz] += -dt·M_L
 *
 *    Row ψ_xSxx (以 lumped mass 为对角)：
 *        A[ψ , ψ ] += M_L(node k)
 *        A[ψ , xx ] += (dt/2)·d_x(x_k)·G_S       (G_S = ∫ N_i ∂N_j/∂α)
 *
 *  声源施加：在 (isou, jsou) 节点处，对 σxx / σzz 的 RHS 加 dt·sou_t。
 * ===========================================================================*/
#ifndef FEM_ASSEMBLER_H
#define FEM_ASSEMBLER_H

#include <vector>

#include "config.h"
#include "material.h"
#include "mesh.h"
#include "sparse_matrix.h"

namespace femasm {

struct FemSystem {
    femla::SparseMatrixCSR A;    // (M + dt K + PML 耦合)，与时间无关，只装配一次
    femla::SparseMatrixCSR M;    // 用于每步 RHS 的 mass 算子（含 ψ 行的 M_L 对角）
    int Ntot = 0;

    // ---- PML 相关（每节点一维缓存） ----
    std::vector<double> dampx_at_ix;   // d_x(ix·dx),  size Nx
    std::vector<double> dampz_at_iz;   // d_z(iz·dz),  size Nz
    std::vector<double> Ex_at_ix;      // exp(-d_x·dt), size Nx
    std::vector<double> Ez_at_iz;      // exp(-d_z·dt), size Nz
    int    pml_active_x = 0;
    int    pml_active_z = 0;
};

// 组装 A 与 M。全部按块 (13 * Nnodes) 布局，DOF 排序同 config.h。
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  FemSystem &out);

// 组装本步右端 b。分两部分：
//   1) 主 5 场：CN 时 b = 2·M·q^n - A·q^n；BE 时 b = M·q^n
//                加上 σxx, σzz 上的 dt·sou_t 声源（与 FDM 相同，不乘 M_L）
//   2) 8 个 ψ 场：b = E·(M_L·ψ^n  -  (dt/2)·d·G·field^n)
//        其中 (G·field)_k 通过一次 element loop 累加计算（"nodal gradient"）
void assemble_rhs(const FemSystem &sys,
                  const femmesh::Mesh &mesh,
                  const std::vector<double> &q_n,
                  double dt,
                  double sou_t,
                  std::vector<double> &b_out);

} // namespace femasm

#endif // FEM_ASSEMBLER_H
