/* =============================================================================
 *  assembler.h
 *  ---------------------------------------------------------------------------
 *  从「PDE 弱形式 + P1 三角单元」得到全局
 *
 *      A x = b
 *
 *  的组装层。这是 README §3.2 中所说的“离散层 (Discretization Layer)”。
 *
 *  重要设计（README §9.2）：
 *      - build_system(...) 只填 A 和 M；
 *      - assemble_rhs(...)  只根据 q^n 与 sou_t 填 b；
 *      - 求解器抽象保持接口一致，不入本层。
 *
 *  块结构（5 场，都用连续 P1）:
 *      q = [ vx ; vz ; sxx ; szz ; sxz ]
 *      M = block-diag(Mρ, Mρ, M0, M0, M0)
 *      K = [  0     0     Gx    0    Gz  ]         Gx = ∫ ∂N_i/∂x  N_j dA
 *          [  0     0     0    Gz    Gx  ]         Gz = ∫ ∂N_i/∂z  N_j dA
 *          [-DxC  -DzL    0     0     0  ]         DxC_ij = ∫ N_i C ∂N_j/∂x dA
 *          [-DxL  -DzC    0     0     0  ]         DzL_ij = ∫ N_i λ ∂N_j/∂z dA  等
 *          [-DzM  -DxM    0     0     0  ]
 *      A = M + dt * K            （Backward Euler）
 *      b = M q^n + dt * f^{n+1}
 *
 *  声源施加 (README §7.3 boundary 项省略 => 自然自由面):
 *      在 (isou, jsou) 节点行处，对 σxx / σzz 的 RHS 加 dt * sou_t（同 FDM）。
 * ===========================================================================*/
#ifndef FEM_ASSEMBLER_H
#define FEM_ASSEMBLER_H

#include <vector>

#include "config.h"
#include "material.h"
#include "mesh.h"
#include "sparse_matrix.h"

namespace femasm {

// A = M + dt K     和  M 都返回（M 用于每步 RHS）
struct FemSystem {
    femla::SparseMatrixCSR A;     // (M + dt K)，与时间无关，只装配一次
    femla::SparseMatrixCSR M;     // 只是 mass 部分，用于每步右端 M * q^n
    int Ntot = 0;
};

// 组装 A 与 M。全部按块 (5 * Nnodes) 布局，DOF 排序同 config.h。
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  FemSystem &out);

// 组装本步右端 b = M * q^n + dt * f^{n+1}
//   q^n 传入按“5 * Nnodes”的整个块向量
//   sou_t : 本时间步 sou_t = -(1-2a^2)e^{-a^2}，与 FDM 相同
void assemble_rhs(const femla::SparseMatrixCSR &M,
                  const std::vector<double>    &q_n,
                  double dt,
                  double sou_t,
                  std::vector<double>          &b_out);

} // namespace femasm

#endif // FEM_ASSEMBLER_H
