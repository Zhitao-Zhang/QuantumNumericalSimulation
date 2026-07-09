/* =============================================================================
 *  sparse_matrix.cpp
 * ===========================================================================*/
#include "sparse_matrix.h"

namespace femla {

void SparseMatrixCSR::build_from_rows(int n_, std::vector<RowMap> &rows)
{
    n = n_;
    rowptr.assign(n + 1, 0);
    for (int i = 0; i < n; i++)
        rowptr[i + 1] = rowptr[i] + (int)rows[i].size();
    const int nnz_ = rowptr[n];
    colind.resize(nnz_);
    values.resize(nnz_);
    diag.assign(n, 0.0);
    int k = 0;
    for (int i = 0; i < n; i++) {
        for (auto &kv : rows[i]) {
            colind[k] = kv.first;
            values[k] = kv.second;
            if (kv.first == i) diag[i] = kv.second;
            k++;
        }
    }
    // 防 diag=0 → Jacobi 会崩，保底填 1.0
    for (int i = 0; i < n; i++)
        if (diag[i] == 0.0) diag[i] = 1.0;
}

void SparseMatrixCSR::multiply(const std::vector<double> &x,
                               std::vector<double> &y) const
{
    if ((int)y.size() != n) y.assign(n, 0.0);
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int k = rowptr[i]; k < rowptr[i + 1]; k++)
            s += values[k] * x[colind[k]];
        y[i] = s;
    }
}

void SparseMatrixCSR::clear()
{
    n = 0;
    rowptr.clear(); colind.clear(); values.clear(); diag.clear();
}

} // namespace femla
