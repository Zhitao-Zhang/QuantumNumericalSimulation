/* =============================================================================
 *  sparse_matrix.h
 *  ---------------------------------------------------------------------------
 *  CSR 稀疏矩阵 + 装配临时结构 (RowMap)。
 *
 *  README §9 明确要求：
 *      LinearSystem = { A, b, x }
 *      assemble_linear_system(...) 只填 A, b
 *      solve_linear_system(...)  只求 x
 *  因此本模块**不含**任何 PDE / FEM 语义，只提供矩阵数据结构和运算原语。
 * ===========================================================================*/
#ifndef FEM_SPARSE_MATRIX_H
#define FEM_SPARSE_MATRIX_H

#include <map>
#include <vector>

namespace femla {

// 装配阶段的行式稀疏容器
using RowMap = std::map<int, double>;

// 只读 CSR
struct SparseMatrixCSR {
    int                 n     = 0;
    std::vector<int>    rowptr;
    std::vector<int>    colind;
    std::vector<double> values;
    std::vector<double> diag;    // 对角元（Jacobi 预条件用）

    // 从 std::vector<RowMap> 一次性构造 CSR
    void build_from_rows(int n_, std::vector<RowMap> &rows);

    // y = A * x
    void multiply(const std::vector<double> &x, std::vector<double> &y) const;

    int   nnz()  const { return (int)values.size(); }
    void  clear();
};

// 统一「线性系统」容器（README §9）
struct LinearSystem {
    SparseMatrixCSR     A;
    std::vector<double> b;
    std::vector<double> x;    // 也承担初值 x0 与解 x 的双重角色
};

} // namespace femla

#endif // FEM_SPARSE_MATRIX_H
