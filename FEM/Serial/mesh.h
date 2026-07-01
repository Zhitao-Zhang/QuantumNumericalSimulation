/* =============================================================================
 *  mesh.h
 *  ---------------------------------------------------------------------------
 *  结构化三角 P1 网格与单元几何。
 *
 *  网格布局（Nx × Nz 个节点）：
 *      节点 (ix, iz) 位置 (x=ix*dx, z=iz*dz)，node_id = ix + iz*Nx
 *      每个矩形 cell (icx, icz)（SW 角在 (icx,icz)）切成 2 个三角形：
 *          T_lower : (icx,   icz)  , (icx+1, icz)  , (icx+1, icz+1)     CCW
 *          T_upper : (icx,   icz)  , (icx+1, icz+1), (icx,   icz+1)     CCW
 *      两个三角形共用 cell 材料。
 *
 *  P1 形函数几何量（每个三角形常数）：
 *      A_signed = 0.5 * ((x2-x1)(z3-z1) - (x3-x1)(z2-z1))
 *      b_i = z_j - z_k, c_i = x_k - x_j  (i,j,k 循环)
 *      ∂N_i/∂x = b_i / (2 A_signed)
 *      ∂N_i/∂z = c_i / (2 A_signed)
 *      |A|      = |A_signed|
 *
 *  README §7.1: 采用三角形非结构化描述，即使结构化网格也统一保留
 *      - 节点坐标表
 *      - 三角形连接表
 *      - 基于参考三角形的数值积分（本文件里给出闭式公式）
 * ===========================================================================*/
#ifndef FEM_MESH_H
#define FEM_MESH_H

#include <cstddef>
#include <vector>

#include "config.h"

namespace femmesh {

// 单个三角形单元
struct Triangle {
    int  nodes[3];   // 三个顶点的 global node id (CCW 顺序)
    int  cell_id;    // 所属 cell 编号 (用于查材料)
    double A;        // |面积|
    double A_signed; // 有符号面积（这里一律 > 0）
    double dN_dx[3]; // 形函数梯度 x 分量（每个三角形常数）
    double dN_dz[3];
};

class Mesh {
public:
    Mesh();

    // 节点坐标（大小 = Nnodes，先 ix 后 iz 遍历）
    std::vector<double> x_node;
    std::vector<double> z_node;

    // 三角形连接表（大小 = Ntri）
    std::vector<Triangle> triangles;

    // 单元 -> cell 属性索引（大小 = Ncells），仅结构上暴露
    std::vector<int> cell_id_of_tri;

    // 每个 cell 的 (icx, icz) —— 供 material 层读取
    static void cell_ij(int cell_id, int &icx, int &icz);

    int n_nodes() const { return (int)x_node.size(); }
    int n_tris () const { return (int)triangles.size(); }

private:
    void build_nodes();
    void build_triangles();
};

// 计算 P1 三角形的形函数梯度与面积
// - 假定 3 个顶点已经是 CCW 顺序（本 mesh 生成器保证）
// - 用有符号面积作正定分母；若外部传入 CW 三角形，梯度公式仍成立，但 |A| 用于积分。
inline void compute_p1_geometry(const double *px, const double *pz,
                                double &A_signed, double &A_abs,
                                double dN_dx[3], double dN_dz[3])
{
    const double x1 = px[0], x2 = px[1], x3 = px[2];
    const double z1 = pz[0], z2 = pz[1], z3 = pz[2];

    const double det = (x2 - x1) * (z3 - z1) - (x3 - x1) * (z2 - z1);
    A_signed = 0.5 * det;
    A_abs    = A_signed >= 0.0 ? A_signed : -A_signed;

    // b_i = z_j - z_k, c_i = x_k - x_j （循环 (1,2,3)->(j,k) = (2,3), (3,1), (1,2)）
    const double b1 = z2 - z3, b2 = z3 - z1, b3 = z1 - z2;
    const double c1 = x3 - x2, c2 = x1 - x3, c3 = x2 - x1;

    const double inv2A = 1.0 / (2.0 * A_signed);
    dN_dx[0] = b1 * inv2A;
    dN_dx[1] = b2 * inv2A;
    dN_dx[2] = b3 * inv2A;
    dN_dz[0] = c1 * inv2A;
    dN_dz[1] = c2 * inv2A;
    dN_dz[2] = c3 * inv2A;
}

} // namespace femmesh

#endif // FEM_MESH_H
