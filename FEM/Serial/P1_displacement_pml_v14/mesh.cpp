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
#ifdef USE_Q1_MESH
    build_quads();
    printf("[mesh] element type = Q1 quad  (bilinear, 1 element / cell)\n");
#else
    build_triangles();
    printf("[mesh] element type = P1 tri (Union Jack)   (2 elements / cell)\n");
#endif
    build_lumped_mass();
    build_M0_diagonal();

#ifdef USE_Q1_MESH
    printf("[mesh] Nx=%d Nz=%d  Nnodes=%d  Ncells=%d  Nelem(quad)=%d  (dx=%g, dz=%g)\n",
           Nx, Nz, Nnodes, Ncells, (int)quads.size(), dx, dz);
#else
    printf("[mesh] Nx=%d Nz=%d  Nnodes=%d  Ncells=%d  Ntri=%d  (dx=%g, dz=%g)\n",
           Nx, Nz, Nnodes, Ncells, Ntri, dx, dz);
#endif
    printf("[mesh] Domain: x in [0, %g] m, z in [0, %g] m\n",
           (Nx - 1) * dx, (Nz - 1) * dz);

    double ml_min = M_L[0], ml_max = M_L[0], ml_sum = 0.0;
    for (double v : M_L) {
        if (v < ml_min) ml_min = v;
        if (v > ml_max) ml_max = v;
        ml_sum += v;
    }
    printf("[mesh] Lumped mass: min=%g  max=%g  sum=%g  (expected sum = (Nx-1)(Nz-1)*dx*dz = %g)\n",
           ml_min, ml_max, ml_sum, (double)(Nx-1)*(Nz-1)*dx*dz);
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

void Mesh::build_lumped_mass()
{
    // M_L[k] = ∫ N_k dΩ = Σ_{elements ∋ k} ∫_e N_k dxdz
    //   P1 tri : ∫_e N_k dxdz = |A_e|/3           (each of 3 nodes gets |A_e|/3)
    //   Q1 quad: ∫_e N_k dxdz = |A_e|/4 = dx·dz/4 (each of 4 nodes gets |A_e|/4)
    M_L.assign(Nnodes, 0.0);
    if (!triangles.empty()) {
        for (const auto &t : triangles) {
            const double contrib = t.A / 3.0;
            M_L[t.nodes[0]] += contrib;
            M_L[t.nodes[1]] += contrib;
            M_L[t.nodes[2]] += contrib;
        }
    }
    if (!quads.empty()) {
        for (const auto &q : quads) {
            const double contrib = 0.25 * q.dx_val * q.dz_val;
            M_L[q.nodes[0]] += contrib;
            M_L[q.nodes[1]] += contrib;
            M_L[q.nodes[2]] += contrib;
            M_L[q.nodes[3]] += contrib;
        }
    }
}

void Mesh::build_triangles()
{
    triangles.clear();
    triangles.reserve(Ntri);
    cell_id_of_tri.reserve(Ntri);

    // 每个 cell (icx, icz) 切成两个 CCW 三角。
    //   记 SW=(icx,icz), SE=(icx+1,icz), NW=(icx,icz+1), NE=(icx+1,icz+1)
    //
    //   Union Jack 网格：奇偶 cell 交替对角，避免 all-SW-NE 网格的方向偏置
    //     - (icx+icz) 偶：SW-NE 对角
    //         T_a: SW → SE → NE   (CCW)
    //         T_b: SW → NE → NW   (CCW)
    //     - (icx+icz) 奇：NW-SE 对角
    //         T_a: SW → SE → NW   (CCW)
    //         T_b: SE → NE → NW   (CCW)
    //
    //   编译期用 `-DNO_UNION_JACK` 可退回旧的“全 SW-NE 对角”行为。
    //
    //   注: 使用“z 值增大方向”作 CCW 判据；两种拆分对每个三角 A_signed = 0.5 dx dz > 0。

    auto push_tri = [&](int cid, int n0, int n1, int n2)
    {
        Triangle t;
        t.nodes[0] = n0; t.nodes[1] = n1; t.nodes[2] = n2;
        t.cell_id  = cid;
        const double px[3] = { x_node[n0], x_node[n1], x_node[n2] };
        const double pz[3] = { z_node[n0], z_node[n1], z_node[n2] };
        compute_p1_geometry(px, pz, t.A_signed, t.A, t.dN_dx, t.dN_dz);
        triangles.push_back(t);
    };

    for (int icz = 0; icz < NcellZ; icz++) {
        for (int icx = 0; icx < NcellX; icx++) {
            const int cid = icx + icz * NcellX;
            const int sw = (icx    ) + (icz    ) * Nx;
            const int se = (icx + 1) + (icz    ) * Nx;
            const int nw = (icx    ) + (icz + 1) * Nx;
            const int ne = (icx + 1) + (icz + 1) * Nx;

#ifdef NO_UNION_JACK
            // 老版本：全 cell 都 SW-NE 对角
            push_tri(cid, sw, se, ne);
            push_tri(cid, sw, ne, nw);
#else
            const bool parity_even = ((icx + icz) & 1) == 0;
            if (parity_even) {
                // SW-NE 对角
                push_tri(cid, sw, se, ne);
                push_tri(cid, sw, ne, nw);
            } else {
                // NW-SE 对角
                push_tri(cid, sw, se, nw);
                push_tri(cid, se, ne, nw);
            }
#endif
        }
    }
}

void Mesh::build_M0_diagonal()
{
    // Consistent mass 对角：M0_diag[k] = Σ_{e∋k} (|A_e|/6)  for P1 tri
    //                        = Σ_{e∋k} (|A_e|/4)  for Q1 quad
    M0_diag.assign(Nnodes, 0.0);
    if (!triangles.empty()) {
        for (const auto &t : triangles) {
            const double d = t.A / 6.0;
            M0_diag[t.nodes[0]] += d;
            M0_diag[t.nodes[1]] += d;
            M0_diag[t.nodes[2]] += d;
        }
    }
    if (!quads.empty()) {
        for (const auto &q : quads) {
            const double d = 0.25 * q.dx_val * q.dz_val;
            M0_diag[q.nodes[0]] += d;
            M0_diag[q.nodes[1]] += d;
            M0_diag[q.nodes[2]] += d;
            M0_diag[q.nodes[3]] += d;
        }
    }
}

// -----------------------------------------------------------------------------
// Q1 四边形网格生成（USE_Q1_MESH）
//   每个矩形 cell → 1 个 4 节点双线性 quad 单元
//   节点 local 顺序：0=SW, 1=SE, 2=NE, 3=NW
// -----------------------------------------------------------------------------
void Mesh::build_quads()
{
    quads.clear();
    quads.reserve(Ncells);

    for (int icz = 0; icz < NcellZ; icz++) {
        for (int icx = 0; icx < NcellX; icx++) {
            Quad q;
            q.cell_id = icx + icz * NcellX;
            q.nodes[0] = (icx    ) + (icz    ) * Nx;  // SW
            q.nodes[1] = (icx + 1) + (icz    ) * Nx;  // SE
            q.nodes[2] = (icx + 1) + (icz + 1) * Nx;  // NE
            q.nodes[3] = (icx    ) + (icz + 1) * Nx;  // NW
            q.dx_val = dx;
            q.dz_val = dz;
            quads.push_back(q);
        }
    }
}

} // namespace femmesh
