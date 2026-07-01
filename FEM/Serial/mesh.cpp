/* =============================================================================
 *  mesh.cpp
 * ===========================================================================*/
#include "mesh.h"

#include <cstdio>

namespace femmesh {

using femcfg::Nx;
using femcfg::Nz;
using femcfg::dx;
using femcfg::dz;
using femcfg::NcellX;
using femcfg::NcellZ;
using femcfg::Ncells;
using femcfg::Ntri;
using femcfg::Nnodes;

Mesh::Mesh()
{
    build_nodes();
    build_triangles();

    printf("[mesh] Nx=%d Nz=%d  Nnodes=%d  Ncells=%d  Ntri=%d  (dx=%g, dz=%g)\n",
           Nx, Nz, Nnodes, Ncells, Ntri, dx, dz);
    printf("[mesh] Domain: x∈[0, %g] m, z∈[0, %g] m\n",
           (Nx - 1) * dx, (Nz - 1) * dz);
}

void Mesh::cell_ij(int cell_id, int &icx, int &icz)
{
    icz = cell_id / NcellX;
    icx = cell_id - icz * NcellX;
}

void Mesh::build_nodes()
{
    x_node.assign(Nnodes, 0.0);
    z_node.assign(Nnodes, 0.0);
    for (int iz = 0; iz < Nz; iz++) {
        for (int ix = 0; ix < Nx; ix++) {
            const int nid = ix + iz * Nx;
            x_node[nid] = ix * dx;
            z_node[nid] = iz * dz;
        }
    }
}

void Mesh::build_triangles()
{
    triangles.clear();
    triangles.reserve(Ntri);
    cell_id_of_tri.reserve(Ntri);

    // 每个 cell (icx, icz) 切两个 CCW 三角
    //   记 SW=(icx,icz), SE=(icx+1,icz), NW=(icx,icz+1), NE=(icx+1,icz+1)
    //   T_lower: SW, SE, NE   （CCW，因为 z 轴“向上”）
    //   T_upper: SW, NE, NW   （CCW）
    //
    //   注: 这里我们采用“z 值增大方向 = 数学 y 增大方向”约定去判 CCW。
    //   带入 A_signed = 0.5 dx dz > 0，符合 CCW。
    for (int icz = 0; icz < NcellZ; icz++) {
        for (int icx = 0; icx < NcellX; icx++) {
            const int cid = icx + icz * NcellX;
            const int sw = (icx    ) + (icz    ) * Nx;
            const int se = (icx + 1) + (icz    ) * Nx;
            const int nw = (icx    ) + (icz + 1) * Nx;
            const int ne = (icx + 1) + (icz + 1) * Nx;

            // T_lower
            {
                Triangle t;
                t.nodes[0] = sw; t.nodes[1] = se; t.nodes[2] = ne;
                t.cell_id  = cid;
                const double px[3] = { x_node[sw], x_node[se], x_node[ne] };
                const double pz[3] = { z_node[sw], z_node[se], z_node[ne] };
                compute_p1_geometry(px, pz, t.A_signed, t.A, t.dN_dx, t.dN_dz);
                triangles.push_back(t);
            }
            // T_upper
            {
                Triangle t;
                t.nodes[0] = sw; t.nodes[1] = ne; t.nodes[2] = nw;
                t.cell_id  = cid;
                const double px[3] = { x_node[sw], x_node[ne], x_node[nw] };
                const double pz[3] = { z_node[sw], z_node[ne], z_node[nw] };
                compute_p1_geometry(px, pz, t.A_signed, t.A, t.dN_dx, t.dN_dz);
                triangles.push_back(t);
            }
        }
    }
}

} // namespace femmesh
