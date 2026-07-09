/* =============================================================================
 *  material.h
 *  ---------------------------------------------------------------------------
 *  逐 cell 的物性 (ρ, λ, μ)。
 *
 *  与 FDM 一致（Quantum/FDM/ImplicitV2/main.cpp 里的 build_medium）：
 *      cell (icx, icz) 若 iz∈[100, 120) 则取固体 2（高速夹层），否则固体 1。
 *  两个三角形共享同一 cell 的物性。
 *
 *  README §7.2 的建议是“P1 空间上 5 场都用连续 P1”，这与逐 cell 的常物性
 *  完全兼容：每个单元内 (λ+2μ), λ, μ, ρ 都取常数，形函数积分闭式给出。
 *
 *  预留接口:
 *      Material::lambda(icx, icz), mu(icx, icz), rho(icx, icz)
 *      即使后续换成任意 2D 模型（VTI、井孔、断层等），只要重写这三个
 *      函数或直接给数组，装配层无需改动。
 * ===========================================================================*/
#ifndef FEM_MATERIAL_H
#define FEM_MATERIAL_H

#include <vector>

#include "config.h"

namespace femmat {

struct CellProps {
    double rho;
    double lambda;
    double mu;
    double lambda_plus_2mu() const { return lambda + 2.0 * mu; }
};

class Material {
public:
    Material();

    // 按 cell_id 查
    const CellProps& at_cell(int cell_id) const { return cells_[cell_id]; }

    // 按 (icx, icz) 查
    const CellProps& at(int icx, int icz) const {
        return cells_[icx + icz * femcfg::NcellX];
    }

    // 供调试输出：统计每种介质的 cell 数
    int  n_layer_cells() const { return n_layer_cells_; }
    int  n_background() const  { return femcfg::Ncells - n_layer_cells_; }

    void print_summary() const;

private:
    std::vector<CellProps> cells_;
    int n_layer_cells_ = 0;
};

} // namespace femmat

#endif // FEM_MATERIAL_H
