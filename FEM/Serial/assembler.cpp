/* =============================================================================
 *  assembler.cpp
 *  ---------------------------------------------------------------------------
 *  P1 三角单元的闭式积分公式（|A| 为该单元面积）：
 *
 *      ∫_e N_i N_j dA          = (|A|/12) * (1 + δ_ij)     (consistent mass)
 *      ∫_e N_j dA              = |A|/3
 *      ∂N_i/∂x                 = 常数 (dN_dx[i])
 *      ∫_e (∂N_i/∂x) N_j dA    = dN_dx[i] * |A|/3
 *      ∫_e N_i (∂N_j/∂x) dA    = dN_dx[j] * |A|/3          (= 上式转置)
 *
 *  说明：
 *      - Mρ^e_ij  = ρ_e * (|A|/12) * (1 + δ_ij)             ← 密度加权 mass
 *      - M0^e_ij  = (|A|/12) * (1 + δ_ij)                    ← unit mass
 *      - Gx^e_ij  = dN_dx[i] * |A|/3      (行 i 与列 j 无关)
 *      - Gz^e_ij  = dN_dz[i] * |A|/3
 *      - DxC^e_ij = C_e * dN_dx[j] * |A|/3    (与行 i 无关)
 *      - DzL^e_ij = L_e * dN_dz[j] * |A|/3    …等等
 *
 *  这些量都是 3x3 局部矩阵；下面通过一个 (row, col, val) accumulator
 *  三重循环把它们装配到 std::map<int,double>[Ntot] 中。
 * ===========================================================================*/
#include "assembler.h"

#include <cstdio>
#include <cstring>
#include <ctime>

namespace femasm {

using femcfg::Nnodes;
using femcfg::Ntot;
using femcfg::off_vx;
using femcfg::off_vz;
using femcfg::off_xx;
using femcfg::off_zz;
using femcfg::off_xz;
using femcfg::isou;
using femcfg::jsou;
using femcfg::dof_xx;
using femcfg::dof_zz;

using femla::RowMap;
using femla::SparseMatrixCSR;
using femmesh::Triangle;

static inline int gid_vx(int node) { return off_vx + node; }
static inline int gid_vz(int node) { return off_vz + node; }
static inline int gid_xx(int node) { return off_xx + node; }
static inline int gid_zz(int node) { return off_zz + node; }
static inline int gid_xz(int node) { return off_xz + node; }

// -----------------------------------------------------------------------------
// 把两个 CSR 一起装配（A 与 M）。装配顺序：先做 mass 块，再做 K 耦合块。
// -----------------------------------------------------------------------------
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  FemSystem &out)
{
    printf("[assemble] Ntot=%d  Ntri=%d ... ", Ntot, mesh.n_tris());
    fflush(stdout);
    clock_t tA0 = clock();

    std::vector<RowMap> rowsA(Ntot);
    std::vector<RowMap> rowsM(Ntot);

    // ----- 遍历每个三角形 -----
    const auto &tris = mesh.triangles;
    for (const auto &t : tris) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double rho   = mp.rho;
        const double lamda = mp.lambda;         // L = λ
        const double mu    = mp.mu;             // M = μ
        const double C     = lamda + 2.0 * mu;  // C = λ + 2μ

        const double Aabs = t.A;               // |A|
        const double A_over_3 = Aabs / 3.0;
        const double A_over_12 = Aabs / 12.0;
        const double dNdx[3] = { t.dN_dx[0], t.dN_dx[1], t.dN_dx[2] };
        const double dNdz[3] = { t.dN_dz[0], t.dN_dz[1], t.dN_dz[2] };
        const int n0 = t.nodes[0], n1 = t.nodes[1], n2 = t.nodes[2];
        const int gnode[3] = { n0, n1, n2 };

        // ===================================================================
        // (1) Mass 部分：填入 A 的对角块 (5x)，同时填入独立 M 矩阵
        // ===================================================================
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                const double mij_unit = A_over_12 * (i == j ? 2.0 : 1.0);
                const double mij_rho  = rho * mij_unit;

                const int gi = gnode[i], gj = gnode[j];

                rowsA[gid_vx(gi)][gid_vx(gj)] += mij_rho;   // Mρ 块
                rowsA[gid_vz(gi)][gid_vz(gj)] += mij_rho;
                rowsA[gid_xx(gi)][gid_xx(gj)] += mij_unit;  // M0 块
                rowsA[gid_zz(gi)][gid_zz(gj)] += mij_unit;
                rowsA[gid_xz(gi)][gid_xz(gj)] += mij_unit;

                rowsM[gid_vx(gi)][gid_vx(gj)] += mij_rho;
                rowsM[gid_vz(gi)][gid_vz(gj)] += mij_rho;
                rowsM[gid_xx(gi)][gid_xx(gj)] += mij_unit;
                rowsM[gid_zz(gi)][gid_zz(gj)] += mij_unit;
                rowsM[gid_xz(gi)][gid_xz(gj)] += mij_unit;
            }
        }

        // ===================================================================
        // (2) 耦合部分 dt * K：
        //
        //   vx 行 (row = gid_vx(vi)):
        //       +dt * Gx^e_ij * σxx_j    => A[vx_i, xx_j] += dt * dN_dx[i] * A/3
        //       +dt * Gz^e_ij * σxz_j    => A[vx_i, xz_j] += dt * dN_dz[i] * A/3
        //
        //   vz 行:
        //       +dt * Gz^e_ij * σzz_j    => A[vz_i, zz_j] += dt * dN_dz[i] * A/3
        //       +dt * Gx^e_ij * σxz_j    => A[vz_i, xz_j] += dt * dN_dx[i] * A/3
        //
        //   σxx 行 (row = gid_xx(vi)):
        //       -dt * DxC^e_ij * vx_j    => A[xx_i, vx_j] += -dt * C * dN_dx[j] * A/3
        //       -dt * DzL^e_ij * vz_j    => A[xx_i, vz_j] += -dt * L * dN_dz[j] * A/3
        //
        //   σzz 行:
        //       -dt * DxL^e_ij * vx_j
        //       -dt * DzC^e_ij * vz_j
        //
        //   σxz 行:
        //       -dt * DzM^e_ij * vx_j    => A[xz_i, vx_j] += -dt * M * dN_dz[j] * A/3
        //       -dt * DxM^e_ij * vz_j    => A[xz_i, vz_j] += -dt * M * dN_dx[j] * A/3
        // ===================================================================
        const double dt_A3 = dt * A_over_3;

        for (int i = 0; i < 3; i++) {
            const int gi = gnode[i];
            const double dt_dx_i = dt_A3 * dNdx[i];   // 用于 vx/vz 行（行系数）
            const double dt_dz_i = dt_A3 * dNdz[i];

            for (int j = 0; j < 3; j++) {
                const int gj = gnode[j];
                const double dt_dx_j = dt_A3 * dNdx[j];  // 用于 σ 行（列系数）
                const double dt_dz_j = dt_A3 * dNdz[j];

                // vx 行
                rowsA[gid_vx(gi)][gid_xx(gj)] += +dt_dx_i;
                rowsA[gid_vx(gi)][gid_xz(gj)] += +dt_dz_i;

                // vz 行
                rowsA[gid_vz(gi)][gid_zz(gj)] += +dt_dz_i;
                rowsA[gid_vz(gi)][gid_xz(gj)] += +dt_dx_i;

                // σxx 行
                rowsA[gid_xx(gi)][gid_vx(gj)] += -C     * dt_dx_j;
                rowsA[gid_xx(gi)][gid_vz(gj)] += -lamda * dt_dz_j;

                // σzz 行
                rowsA[gid_zz(gi)][gid_vx(gj)] += -lamda * dt_dx_j;
                rowsA[gid_zz(gi)][gid_vz(gj)] += -C     * dt_dz_j;

                // σxz 行
                rowsA[gid_xz(gi)][gid_vx(gj)] += -mu * dt_dz_j;
                rowsA[gid_xz(gi)][gid_vz(gj)] += -mu * dt_dx_j;
            }
        }
    }

    // ----- 构造 CSR -----
    out.Ntot = Ntot;
    out.A.build_from_rows(Ntot, rowsA);
    out.M.build_from_rows(Ntot, rowsM);

    rowsA.clear(); rowsA.shrink_to_fit();
    rowsM.clear(); rowsM.shrink_to_fit();

    printf("done. nnz(A)=%d, nnz(M)=%d, %.2f s\n",
           out.A.nnz(), out.M.nnz(),
           (double)(clock() - tA0) / CLOCKS_PER_SEC);
}

// -----------------------------------------------------------------------------
// 组装本步 b = M q^n + dt f
//    - M q^n 是一次稀疏矩阵向量乘
//    - 声源加到 σxx / σzz 的 RHS
// -----------------------------------------------------------------------------
void assemble_rhs(const SparseMatrixCSR &M,
                  const std::vector<double> &q_n,
                  double dt,
                  double sou_t,
                  std::vector<double> &b_out)
{
    if ((int)b_out.size() != Ntot) b_out.assign(Ntot, 0.0);
    M.multiply(q_n, b_out);

    // 声源施加：与 FDM 相同 (isou, jsou) 位置的 σxx / σzz
    const int src_node = isou + jsou * femcfg::Nx;
    b_out[gid_xx(src_node)] += dt * sou_t;
    b_out[gid_zz(src_node)] += dt * sou_t;
}

} // namespace femasm
