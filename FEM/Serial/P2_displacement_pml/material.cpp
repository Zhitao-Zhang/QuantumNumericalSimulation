/* =============================================================================
 *  material.cpp
 * ===========================================================================*/
#include "material.h"

#include <cstdio>

namespace femmat {

using femcfg::NcellX;
using femcfg::NcellZ;
using femcfg::Ncells;
using femcfg::Vp1;  using femcfg::Vs1;  using femcfg::rho1;
using femcfg::Vp2;  using femcfg::Vs2;  using femcfg::rho2;
using femcfg::borehole_iz_lo;
using femcfg::borehole_iz_hi;
using femcfg::NBH_val;

static CellProps make_props(double rho, double vp, double vs)
{
    CellProps p;
    p.rho    = rho;
    p.mu     = rho * vs * vs;
    p.lambda = rho * vp * vp - 2.0 * p.mu;
    return p;
}

Material::Material()
{
    cells_.assign(Ncells, CellProps{});
    n_layer_cells_ = 0;

    const CellProps bg    = make_props(rho1, Vp1, Vs1);
    const CellProps layer = make_props(rho2, Vp2, Vs2);

    for (int icz = 0; icz < NcellZ; icz++) {
        for (int icx = 0; icx < NcellX; icx++) {
            const bool in_layer = (NBH_val > 0)
                                && (icz >= borehole_iz_lo)
                                && (icz <  borehole_iz_hi);
            cells_[icx + icz * NcellX] = in_layer ? layer : bg;
            if (in_layer) n_layer_cells_++;
        }
    }
}

void Material::print_summary() const
{
    printf("[material] Background : Vp=%g Vs=%g rho=%g  (mu=%.3e, lambda=%.3e)\n",
           Vp1, Vs1, rho1,
           rho1 * Vs1 * Vs1,
           rho1 * Vp1 * Vp1 - 2.0 * rho1 * Vs1 * Vs1);
    if (NBH_val > 0) {
        printf("[material] High-V layer: Vp=%g Vs=%g rho=%g  (mu=%.3e, lambda=%.3e)\n",
               Vp2, Vs2, rho2,
               rho2 * Vs2 * Vs2,
               rho2 * Vp2 * Vp2 - 2.0 * rho2 * Vs2 * Vs2);
        printf("[material] Layer occupies iz∈[%d, %d)  (%d / %d cells, %.2f%%)\n",
               borehole_iz_lo, borehole_iz_hi,
               n_layer_cells_, Ncells,
               100.0 * n_layer_cells_ / (double)Ncells);
    } else {
        printf("[material] NBH=0  -> homogeneous medium\n");
    }
}

} // namespace femmat
