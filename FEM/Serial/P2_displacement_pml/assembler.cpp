/* =============================================================================
 *  assembler.cpp   —  位移二阶格式装配（方案 B / SESSION 7）
 *  ---------------------------------------------------------------------------
 *  组装 K（标准刚度阵）、M（consistent mass）、C（吸收层质量比例阻尼），
 *  并合成 Newmark 有效矩阵  A = M + γ·dt·C + β·dt²·K。
 *
 *  为什么这样能根治菱形：K = ∫ Bᵀ D B 里梯度以 ∫∇N·∇N 出现，是连中心节点的
 *  正规二阶差分（不像 velocity-stress collocated 那样跳中心 → 2Δx 各向异性）。
 * ===========================================================================*/
#include "assembler.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <map>
#include <utility>
#include <vector>

namespace femasm {

using femcfg::Nx;
using femcfg::Nz;
using femcfg::Nnodes;
using femcfg::Ntot;
using femcfg::off_ux;
using femcfg::off_uz;
using femcfg::isou;
using femcfg::jsou;
using femcfg::dx;
using femcfg::dz;
using femcfg::NABS;

using femla::RowMap;

// 局部 dof (顶点 a ∈ {0,1,2}, 分量 comp ∈ {0=ux,1=uz}) → 全局 dof
static inline int gdof(int gnode, int comp)
{
    return (comp == 0 ? off_ux : off_uz) + gnode;
}

// -----------------------------------------------------------------------------
void build_system(const femmesh::Mesh &mesh,
                  const femmat::Material &mat,
                  double dt,
                  ElastSystem &out)
{
#ifdef USE_P1_MESH
    printf("[assemble] Ntot=%d  Ntri(P1)=%d  (P1 tri, displacement 2nd-order, Newmark)\n",
           Ntot, mesh.n_tris());
#else
    printf("[assemble] Ntot=%d  Ntri6(P2)=%d  (P2 tri, displacement 2nd-order, Newmark)\n",
           Ntot, mesh.n_tris6());
#endif
    fflush(stdout);
    clock_t tA0 = clock();

    const double beta  = out.beta;
    const double gamma = out.gamma;

    // ---- PML 方向性剖面（按细网格 IX∈[0,MX), IZ∈[0,MZ) 采样）----
    //   物理坐标 = IX·hgrid；P1 时 MX=Nx、hgrid=dx，退化回顶点级。
    using femcfg::MX; using femcfg::MZ; using femcfg::hgrid; using femcfg::hgridz;
    out.dampx_ix.assign(MX, 0.0); out.dampz_iz.assign(MZ, 0.0);
    out.Ex_ix.assign(MX, 1.0);    out.Ez_iz.assign(MZ, 1.0);
    out.facx_ix.assign(MX, dt);   out.facz_iz.assign(MZ, dt);
    double dmax = 0.0; int n_pml = 0;
    for (int IX = 0; IX < MX; IX++) {
        const double d = femcfg::dampx_at(IX * hgrid);
        out.dampx_ix[IX] = d;
        if (d > 0.0) { out.Ex_ix[IX] = std::exp(-d * dt);
                       out.facx_ix[IX] = (1.0 - out.Ex_ix[IX]) / d; }
        if (d > dmax) dmax = d;
    }
    for (int IZ = 0; IZ < MZ; IZ++) {
        const double d = femcfg::dampz_at(IZ * hgridz);
        out.dampz_iz[IZ] = d;
        if (d > 0.0) { out.Ez_iz[IZ] = std::exp(-d * dt);
                       out.facz_iz[IZ] = (1.0 - out.Ez_iz[IZ]) / d; }
        if (d > dmax) dmax = d;
    }
    for (int IZ = 0; IZ < MZ; IZ++)
        for (int IX = 0; IX < MX; IX++)
            if (out.dampx_ix[IX] > 0.0 || out.dampz_iz[IZ] > 0.0) n_pml++;
    out.pml_nodes = n_pml;
    out.damp_max  = dmax;

    std::vector<RowMap> rowsA(Ntot);
    std::vector<RowMap> rowsKeff(Ntot);   // 标准 K + 角区弹簧 P（进 A、组 b）
    std::vector<RowMap> rowsC(Ntot);      // PML 阻尼
    std::vector<RowMap> rowsKxx(Ntot);    // x 向刚度块（驱动 ψ_xx）
    std::vector<RowMap> rowsKzz(Ntot);    // z 向刚度块（驱动 ψ_zz）

#ifdef USE_P1_MESH
    // =========================== P1 装配（3 节点，闭式） =========================
    for (const auto &t : mesh.triangles) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double rho = mp.rho;
        const double lam = mp.lambda;
        const double mu  = mp.mu;
        const double Cc  = lam + 2.0 * mu;               // λ+2μ

        const double A_e = t.A;
        const double dNx[3] = { t.dN_dx[0], t.dN_dx[1], t.dN_dx[2] };
        const double dNz[3] = { t.dN_dz[0], t.dN_dz[1], t.dN_dz[2] };
        const int    gn [3] = { t.nodes[0], t.nodes[1], t.nodes[2] };

        // 单元级 PML 系数（节点平均）：
        //   d_sum = <d_x + d_z>（用于 C = ∫ρ(d_x+d_z)NN）
        //   d_pro = <d_x·d_z>  （用于 P = ∫ρ d_x d_z NN，仅角区非零）
        double d_sum = 0.0, d_pro = 0.0;
        if (NABS > 0) {
            for (int a = 0; a < 3; a++) {
                const double dxk = out.dampx_ix[gn[a] % MX];
                const double dzk = out.dampz_iz[gn[a] / MX];
                d_sum += (dxk + dzk);
                d_pro += (dxk * dzk);
            }
            d_sum /= 3.0; d_pro /= 3.0;
        }

        // B 矩阵 (3×6)
        double B[3][6] = {{0}};
        for (int a = 0; a < 3; a++) {
            B[0][2*a+0] = dNx[a];
            B[1][2*a+1] = dNz[a];
            B[2][2*a+0] = dNz[a];
            B[2][2*a+1] = dNx[a];
        }
        const double D[3][3] = { { Cc,  lam, 0.0 },
                                 { lam, Cc,  0.0 },
                                 { 0.0, 0.0, mu  } };
        double DB[3][6];
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 6; c++) {
                double s = 0.0;
                for (int k = 0; k < 3; k++) s += D[r][k] * B[k][c];
                DB[r][c] = s;
            }
        double Ke[6][6];                                 // 完整单元刚度
        for (int p = 0; p < 6; p++)
            for (int q = 0; q < 6; q++) {
                double s = 0.0;
                for (int r = 0; r < 3; r++) s += B[r][p] * DB[r][q];
                Ke[p][q] = A_e * s;
            }

        // 混合质量矩阵（抗 P1 频散）：M = (1-α)·consistent + α·lumped
        //   consistent: diag=A_e/6, off=A_e/12   lumped: diag=A_e/3, off=0
        //   α=0 纯 consistent（默认历史行为）；α=1 纯 lumped。
        //   两者 P1 频散误差符号相反，适当 α 抵消主导误差、压制频散尾。
        //   仍对称 SPD（lumped 对角正、consistent SPD）→ A 仍对称 SPD，不影响 CIM。
        const double A12  = A_e / 12.0;
        const double mblend = femcfg::mass_blend;
        const double m_diag_c = A_e / 6.0, m_off_c = A_e / 12.0;   // consistent
        const double m_diag_l = A_e / 3.0;                         // lumped diag
        const double m_diag_blend = (1.0 - mblend) * m_diag_c + mblend * m_diag_l;
        const double m_off_blend  = (1.0 - mblend) * m_off_c;
        (void)A12;

        for (int a = 0; a < 3; a++) {
            for (int b = 0; b < 3; b++) {
                const double m0 = (a == b) ? m_diag_blend : m_off_blend;
                const double me = rho * m0;              // 质量
                const double ce = d_sum * me;            // PML 阻尼 C
                const double pe = d_pro * me;            // PML 弹簧 P

                // 分量对角块（ux-ux, uz-uz）
                for (int comp = 0; comp < 2; comp++) {
                    const int gi = gdof(gn[a], comp);
                    const int gj = gdof(gn[b], comp);
                    const double kij = Ke[2*a+comp][2*b+comp];   // 完整刚度该分量项
                    const double keff = kij + pe;                // K_eff = K + P
                    rowsKeff[gi][gj] += keff;
                    rowsC   [gi][gj] += ce;
                    rowsA   [gi][gj] += me + gamma * dt * ce + beta * dt * dt * keff;

                    // K_xx / K_zz 分量对角部分（驱动 ψ，不进 A）
                    if (NABS > 0) {
                        const double coef_xx = (comp == 0) ? Cc : mu;   // ux用λ+2μ, uz用μ
                        const double coef_zz = (comp == 0) ? mu : Cc;
                        rowsKxx[gi][gj] += A_e * coef_xx * dNx[a] * dNx[b];
                        rowsKzz[gi][gj] += A_e * coef_zz * dNz[a] * dNz[b];
                    }
                }
                // 分量交叉项 (ux-uz, uz-ux)：只在 K_eff 里（M、C、P 无交叉）
                {
                    const int gi_x = gdof(gn[a], 0);
                    const int gj_z = gdof(gn[b], 1);
                    const double kxz = Ke[2*a+0][2*b+1];
                    rowsKeff[gi_x][gj_z] += kxz;
                    rowsA   [gi_x][gj_z] += beta * dt * dt * kxz;

                    const int gi_z = gdof(gn[a], 1);
                    const int gj_x = gdof(gn[b], 0);
                    const double kzx = Ke[2*a+1][2*b+0];
                    rowsKeff[gi_z][gj_x] += kzx;
                    rowsA   [gi_z][gj_x] += beta * dt * dt * kzx;
                }
            }
        }
    }
#else
    // =========================== P2 装配（6 节点，高斯求积） =====================
    //   参考三角面积坐标 (L0,L1,L2)，L0=1-ξ-η, L1=ξ, L2=η。
    //   P2 形函数：顶点 N_i=L_i(2L_i-1)，边中点 N_{3}=4L0L1, N_{4}=4L1L2, N_{5}=4L2L0。
    //   ∂N/∂x = Σ (∂N/∂L_k)·(∂L_k/∂x)，∂L_k/∂x = dL_dx[k]（常数）。
    //   积分用 Strang 6 点 4 阶规则（对 mass 的四次、K 的二次均精确）。
    // 6 点规则（面积坐标 b0=1-b1-b2）：两组各 3 点
    static const double gpa = 0.108103018168070, gpb = 0.445948490915965;  // w=0.223381589678011
    static const double gpc = 0.816847572980459, gpd = 0.091576213509771;  // w=0.109951743655322
    struct GP { double L0,L1,L2,w; };
    const GP gps[6] = {
        {1-gpb-gpb, gpb, gpb, 0.223381589678011},
        {1-gpa-gpb, gpa, gpb, 0.223381589678011},
        {1-gpb-gpa, gpb, gpa, 0.223381589678011},
        {1-gpd-gpd, gpd, gpd, 0.109951743655322},
        {1-gpc-gpd, gpc, gpd, 0.109951743655322},
        {1-gpd-gpc, gpd, gpc, 0.109951743655322},
    };
    const double mblend = femcfg::mass_blend;   // P2 也支持混合质量（默认 0.5）

    for (const auto &t : mesh.tris6) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double rho = mp.rho, lam = mp.lambda, mu = mp.mu, Cc = lam + 2.0 * mu;
        const double A_e = t.A;
        const double dLx[3] = { t.dL_dx[0], t.dL_dx[1], t.dL_dx[2] };
        const double dLz[3] = { t.dL_dz[0], t.dL_dz[1], t.dL_dz[2] };
        const int    gn[6]  = { t.nodes[0],t.nodes[1],t.nodes[2],t.nodes[3],t.nodes[4],t.nodes[5] };

        // 6 个节点的 PML 阻尼（细网格索引）
        double dxn[6], dzn[6];
        if (NABS > 0)
            for (int a = 0; a < 6; a++) { dxn[a] = out.dampx_ix[gn[a] % MX];
                                          dzn[a] = out.dampz_iz[gn[a] / MX]; }
        else for (int a = 0; a < 6; a++) { dxn[a] = 0; dzn[a] = 0; }

        // 单元 12×12 局部矩阵：Me(质量,标量→每分量), Ke(刚度), Ce, Pe, Kxx, Kzz
        double Mc[6][6] = {{0}};   // consistent mass (标量)
        double Ke[12][12] = {{0}};
        double Ce[6][6] = {{0}}, Pe[6][6] = {{0}};       // 标量（每分量对角）
        double Kxxs[6][6] = {{0}}, Kzzs[6][6] = {{0}};   // 标量刚度块（驱动 ψ）

        for (int g = 0; g < 6; g++) {
            const double L0 = gps[g].L0, L1 = gps[g].L1, L2 = gps[g].L2;
            const double wA = gps[g].w * A_e;            // 权 × 面积
            const double Lv[3] = { L0, L1, L2 };
            // 形函数 N[6]
            double N[6] = { L0*(2*L0-1), L1*(2*L1-1), L2*(2*L2-1),
                            4*L0*L1, 4*L1*L2, 4*L2*L0 };
            // ∂N/∂L_k （3 个面积坐标方向）
            //   顶点 i: dNi/dLi = 4Li-1，其余 0
            //   边中点: N3=4L0L1 → dN3/dL0=4L1, dN3/dL1=4L0；类推
            double dNdL[6][3] = {{0}};
            dNdL[0][0] = 4*L0-1; dNdL[1][1] = 4*L1-1; dNdL[2][2] = 4*L2-1;
            dNdL[3][0] = 4*L1; dNdL[3][1] = 4*L0;
            dNdL[4][1] = 4*L2; dNdL[4][2] = 4*L1;
            dNdL[5][0] = 4*L2; dNdL[5][2] = 4*L0;
            // ∂N/∂x, ∂N/∂z = Σ_k dNdL[·][k]·dL{x,z}[k]
            double Nx_[6], Nz_[6];
            for (int a = 0; a < 6; a++) {
                Nx_[a] = dNdL[a][0]*dLx[0] + dNdL[a][1]*dLx[1] + dNdL[a][2]*dLx[2];
                Nz_[a] = dNdL[a][0]*dLz[0] + dNdL[a][1]*dLz[1] + dNdL[a][2]*dLz[2];
            }
            // 高斯点处 PML 阻尼（用 N 插值节点值）
            double dxg = 0, dzg = 0;
            if (NABS > 0)
                for (int a = 0; a < 6; a++) { dxg += N[a]*dxn[a]; dzg += N[a]*dzn[a]; }
            (void)Lv;

            for (int a = 0; a < 6; a++) {
                for (int b = 0; b < 6; b++) {
                    const double NaNb = N[a]*N[b]*wA;
                    Mc  [a][b] += NaNb;
                    Ce  [a][b] += (dxg + dzg) * NaNb;
                    Pe  [a][b] += (dxg * dzg) * NaNb;
                    Kxxs[a][b] += Nx_[a]*Nx_[b]*wA;   // 待乘 coef(comp)
                    Kzzs[a][b] += Nz_[a]*Nz_[b]*wA;
                    // 完整刚度 Bᵀ D B（12×12）：本地 dof = 2a+comp
                    // εxx=Σ Nx·ux, εzz=Σ Nz·uz, γxz=Σ(Nz·ux+Nx·uz)
                    // Ke[2a+cx][2b+cy] = ∫ Ba^T D Bb
                    const double bax = Nx_[a], baz = Nz_[a];
                    const double bbx = Nx_[b], bbz = Nz_[b];
                    // Ba = [[bax,0],[0,baz],[baz,bax]] (3x2)
                    // D·Bb 列0=(Cc*bbx, lam*bbx? ...) 直接展开：
                    // ux-ux: Cc*bax*bbx + mu*baz*bbz
                    Ke[2*a+0][2*b+0] += (Cc*bax*bbx + mu*baz*bbz) * wA;
                    // ux-uz: lam*bax*bbz + mu*baz*bbx
                    Ke[2*a+0][2*b+1] += (lam*bax*bbz + mu*baz*bbx) * wA;
                    // uz-ux: lam*baz*bbx + mu*bax*bbz
                    Ke[2*a+1][2*b+0] += (lam*baz*bbx + mu*bax*bbz) * wA;
                    // uz-uz: Cc*baz*bbz + mu*bax*bbx
                    Ke[2*a+1][2*b+1] += (Cc*baz*bbz + mu*bax*bbx) * wA;
                }
            }
        }

        // 混合质量：lumped 用行和（对 P2 consistent 行和为正，安全）
        double rowsum[6];
        for (int a = 0; a < 6; a++) { double s=0; for(int b=0;b<6;b++) s+=Mc[a][b]; rowsum[a]=s; }

        for (int a = 0; a < 6; a++) {
            for (int b = 0; b < 6; b++) {
                // 混合 mass（标量）
                double mab = (1.0-mblend)*Mc[a][b] + (a==b ? mblend*rowsum[a] : 0.0);
                const double me = rho * mab;
                const double ce = Ce[a][b];          // 阻尼（consistent 型）
                const double pe = Pe[a][b];          // 弹簧
                for (int comp = 0; comp < 2; comp++) {
                    const int gi = gdof(gn[a], comp);
                    const int gj = gdof(gn[b], comp);
                    const double kij  = Ke[2*a+comp][2*b+comp];
                    const double keff = kij + pe;
                    rowsKeff[gi][gj] += keff;
                    rowsC   [gi][gj] += ce;
                    rowsA   [gi][gj] += me + gamma*dt*ce + beta*dt*dt*keff;
                    if (NABS > 0) {
                        const double coef_xx = (comp==0) ? Cc : mu;
                        const double coef_zz = (comp==0) ? mu : Cc;
                        rowsKxx[gi][gj] += coef_xx * Kxxs[a][b];
                        rowsKzz[gi][gj] += coef_zz * Kzzs[a][b];
                    }
                }
                // 分量交叉（只进 K_eff）
                const int gi_x = gdof(gn[a],0), gj_z = gdof(gn[b],1);
                const double kxz = Ke[2*a+0][2*b+1];
                rowsKeff[gi_x][gj_z] += kxz;
                rowsA   [gi_x][gj_z] += beta*dt*dt*kxz;
                const int gi_z = gdof(gn[a],1), gj_x = gdof(gn[b],0);
                const double kzx = Ke[2*a+1][2*b+0];
                rowsKeff[gi_z][gj_x] += kzx;
                rowsA   [gi_z][gj_x] += beta*dt*dt*kzx;
            }
        }
    }
#endif  // USE_P1_MESH

    // ---- 声源等效节点力系数：F_a = M0·∇N_a ----
    //
    // 【SESSION 11 修复】点源(单顶点)在 P2 上会激发源节点棋盘伪影(顶点/边中点
    // 形函数梯度差异大 → 最高频网格模式)。修法：把矩张量源在源周围顶点上做
    // **高斯空间平滑**(展宽点源、削掉 λ<2dx 分量)。对每个"子源顶点"仍施加
    // 标准 F_a=M0·∇N_a，按高斯权重叠加、权重归一化(∑w=1，保持总矩不变)。
    //   -DSRC_SMOOTH_R=R  平滑半径(顶点数，窗口 (2R+1)²)。默认 P2=3(7×7)，P1=0。
    //   -DSRC_SMOOTH_SIGMA=σ  高斯 σ(顶点单位)。默认 P2=1.5，P1=1.0。
    //   R=0 退回原始集中点源。
#ifdef USE_P1_MESH
    #ifndef SRC_SMOOTH_R
        #define SRC_SMOOTH_R 0
    #endif
    #ifndef SRC_SMOOTH_SIGMA
        #define SRC_SMOOTH_SIGMA 1.0
    #endif
#else
    #ifndef SRC_SMOOTH_R
        #define SRC_SMOOTH_R 3
    #endif
    #ifndef SRC_SMOOTH_SIGMA
        #define SRC_SMOOTH_SIGMA 1.5
    #endif
#endif
    const int    SR    = SRC_SMOOTH_R;
    const double SSIG  = SRC_SMOOTH_SIGMA;
    std::map<int, double> srcmap;

    // 收集子源顶点及归一化高斯权重
    std::vector<std::pair<int,int>> sv;   // (ix,iz)
    std::vector<double> sw;
    double wsum = 0.0;
    for (int dj = -SR; dj <= SR; dj++)
        for (int di = -SR; di <= SR; di++) {
            const int ix = isou + di, iz = jsou + dj;
            if (ix < 1 || ix >= Nx-1 || iz < 1 || iz >= Nz-1) continue;
            const double w = (SR == 0) ? 1.0
                           : std::exp(-(double)(di*di+dj*dj) / (2.0*SSIG*SSIG));
            sv.push_back({ix, iz}); sw.push_back(w); wsum += w;
        }
    for (double &w : sw) w /= wsum;

    // 对每个子源顶点施加标准矩张量等效力 F_a = w·∇N_a
    for (size_t s = 0; s < sv.size(); s++) {
        const int svnode = femcfg::vert_node_id(sv[s].first, sv[s].second);
        const double wgt = sw[s];
#ifdef USE_P1_MESH
        for (const auto &t : mesh.triangles) {
            int loc = -1;
            for (int a = 0; a < 3; a++) if (t.nodes[a] == svnode) loc = a;
            if (loc < 0) continue;
            for (int a = 0; a < 3; a++) {
                srcmap[gdof(t.nodes[a], 0)] += wgt * t.dN_dx[a];
                srcmap[gdof(t.nodes[a], 1)] += wgt * t.dN_dz[a];
            }
        }
#else
        for (const auto &t : mesh.tris6) {
            int loc = -1;
            for (int a = 0; a < 3; a++) if (t.nodes[a] == svnode) loc = a;
            if (loc < 0) continue;
            double L[3] = {0,0,0}; L[loc] = 1.0;
            double dNdL[6][3] = {{0}};
            dNdL[0][0] = 4*L[0]-1; dNdL[1][1] = 4*L[1]-1; dNdL[2][2] = 4*L[2]-1;
            dNdL[3][0] = 4*L[1]; dNdL[3][1] = 4*L[0];
            dNdL[4][1] = 4*L[2]; dNdL[4][2] = 4*L[1];
            dNdL[5][0] = 4*L[2]; dNdL[5][2] = 4*L[0];
            for (int a = 0; a < 6; a++) {
                const double nx = dNdL[a][0]*t.dL_dx[0] + dNdL[a][1]*t.dL_dx[1] + dNdL[a][2]*t.dL_dx[2];
                const double nz = dNdL[a][0]*t.dL_dz[0] + dNdL[a][1]*t.dL_dz[1] + dNdL[a][2]*t.dL_dz[2];
                srcmap[gdof(t.nodes[a], 0)] += wgt * nx;
                srcmap[gdof(t.nodes[a], 1)] += wgt * nz;
            }
        }
#endif
    }
    out.src_dof.clear();
    out.src_coef.clear();
    for (const auto &kv : srcmap) {
        out.src_dof .push_back(kv.first);
        out.src_coef.push_back(kv.second);
    }

    // ---- 构造 CSR ----
    out.Ntot = Ntot;
    out.A.build_from_rows(Ntot, rowsA);
    out.K_eff.build_from_rows(Ntot, rowsKeff);
    out.C.build_from_rows(Ntot, rowsC);
    out.K_xx.build_from_rows(Ntot, rowsKxx);
    out.K_zz.build_from_rows(Ntot, rowsKzz);
    rowsA.clear(); rowsA.shrink_to_fit();
    rowsKeff.clear(); rowsKeff.shrink_to_fit();
    rowsC.clear(); rowsC.shrink_to_fit();
    rowsKxx.clear(); rowsKxx.shrink_to_fit();
    rowsKzz.clear(); rowsKzz.shrink_to_fit();

    printf("[assemble] nnz(A)=%d nnz(Keff)=%d nnz(C)=%d nnz(Kxx)=%d nnz(Kzz)=%d src=%d  (%.2f s)\n",
           out.A.nnz(), out.K_eff.nnz(), out.C.nnz(), out.K_xx.nnz(), out.K_zz.nnz(),
           (int)out.src_dof.size(), (double)(clock() - tA0) / CLOCKS_PER_SEC);
    if (NABS > 0)
        printf("[pml] ADE-PML unsplit  NABS=%d  d0_x=%g d0_z=%g  pml_nodes=%d  max(d)=%g\n",
               NABS, femcfg::compute_d0_x(), femcfg::compute_d0_z(), n_pml, dmax);
    else
        printf("[pml] NABS=0 -> no PML (degenerates to pure displacement, A unchanged)\n");
}

// -----------------------------------------------------------------------------
// PML 记忆变量显式更新 + F_ψ = ψ_xx + ψ_zz
//   ODE:  ψ̇_xx + d_x·ψ_xx = (d_z - d_x)·(K_xx·u)
//         ψ̇_zz + d_z·ψ_zz = (d_x - d_z)·(K_zz·u)
//   指数积分（逐 DOF，d 取该 DOF 所在节点的方向阻尼）：
//         ψ^{n+1} = E·ψ^n + fac·rhs,  E=exp(-d·dt), fac=(1-E)/d (d→0 取 dt)
//   驱动量 K_xx·u_pred / K_zz·u_pred 用 Newmark 预测位移（已知，保持全显式）。
// -----------------------------------------------------------------------------
void pml_update_force(const ElastSystem &sys,
                      const std::vector<double> &u_pred,
                      std::vector<double> &psi_xx,
                      std::vector<double> &psi_zz,
                      std::vector<double> &Fpsi_out)
{
    if ((int)Fpsi_out.size() != Ntot) Fpsi_out.assign(Ntot, 0.0);
    if (femcfg::NABS <= 0) { std::fill(Fpsi_out.begin(), Fpsi_out.end(), 0.0); return; }
    if ((int)psi_xx.size() != Ntot) psi_xx.assign(Ntot, 0.0);
    if ((int)psi_zz.size() != Ntot) psi_zz.assign(Ntot, 0.0);

    static std::vector<double> gxx, gzz;
    sys.K_xx.multiply(u_pred, gxx);
    sys.K_zz.multiply(u_pred, gzz);

    // DOF i 属于块 (comp,node)，node=(i%Nnodes)，其细网格 (IX,IZ) 决定 d_x,d_z
    const int MX = femcfg::MX;
    for (int i = 0; i < Ntot; i++) {
        const int node = i % Nnodes;
        const int IX = node % MX;
        const int IZ = node / MX;
        const double dx_i = sys.dampx_ix[IX];
        const double dz_i = sys.dampz_iz[IZ];

        psi_xx[i] = sys.Ex_ix[IX] * psi_xx[i]
                  + sys.facx_ix[IX] * (dz_i - dx_i) * gxx[i];
        psi_zz[i] = sys.Ez_iz[IZ] * psi_zz[i]
                  + sys.facz_iz[IZ] * (dx_i - dz_i) * gzz[i];
        Fpsi_out[i] = psi_xx[i] + psi_zz[i];
    }
}

// -----------------------------------------------------------------------------
void assemble_source_force(const ElastSystem &sys,
                           double sou_t,
                           std::vector<double> &F_out)
{
    if ((int)F_out.size() != Ntot) F_out.assign(Ntot, 0.0);
    else std::fill(F_out.begin(), F_out.end(), 0.0);

    const double M0 = femcfg::src_scale * sou_t;
    for (size_t k = 0; k < sys.src_dof.size(); k++)
        F_out[sys.src_dof[k]] += M0 * sys.src_coef[k];
}

// -----------------------------------------------------------------------------
void recover_nodal_stress(const femmesh::Mesh &mesh,
                          const femmat::Material &mat,
                          const std::vector<double> &u,
                          std::vector<double> &sxx,
                          std::vector<double> &szz,
                          std::vector<double> &sxz)
{
    sxx.assign(Nnodes, 0.0);
    szz.assign(Nnodes, 0.0);
    sxz.assign(Nnodes, 0.0);
    std::vector<double> wsum(Nnodes, 0.0);

#ifdef USE_P1_MESH
    for (const auto &t : mesh.triangles) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double lam = mp.lambda, mu = mp.mu, Cc = lam + 2.0 * mu;
        const double A_e = t.A;
        const int gn[3] = { t.nodes[0], t.nodes[1], t.nodes[2] };

        double dux_dx = 0, dux_dz = 0, duz_dx = 0, duz_dz = 0;
        for (int a = 0; a < 3; a++) {
            const double ux = u[off_ux + gn[a]];
            const double uz = u[off_uz + gn[a]];
            dux_dx += t.dN_dx[a] * ux;
            dux_dz += t.dN_dz[a] * ux;
            duz_dx += t.dN_dx[a] * uz;
            duz_dz += t.dN_dz[a] * uz;
        }
        const double exx = dux_dx, ezz = duz_dz, gxz = dux_dz + duz_dx;
        const double s_xx = Cc * exx + lam * ezz;
        const double s_zz = lam * exx + Cc * ezz;
        const double s_xz = mu * gxz;

        for (int a = 0; a < 3; a++) {
            const int nd = gn[a];
            sxx[nd] += A_e * s_xx;
            szz[nd] += A_e * s_zz;
            sxz[nd] += A_e * s_xz;
            wsum[nd] += A_e;
        }
    }
#else
    // P2：应力在单元内线性变化。对每个 Tri6，在 3 个顶点的面积坐标处评估
    // ∂N_a/∂x,z（6 节点），得顶点应力，面积加权投影到顶点节点。
    for (const auto &t : mesh.tris6) {
        const femmat::CellProps &mp = mat.at_cell(t.cell_id);
        const double lam = mp.lambda, mu = mp.mu, Cc = lam + 2.0 * mu;
        const double A_e = t.A;
        const int gn[6] = { t.nodes[0],t.nodes[1],t.nodes[2],t.nodes[3],t.nodes[4],t.nodes[5] };
        const double dLx[3] = { t.dL_dx[0], t.dL_dx[1], t.dL_dx[2] };
        const double dLz[3] = { t.dL_dz[0], t.dL_dz[1], t.dL_dz[2] };

        for (int v = 0; v < 3; v++) {
            double L[3] = {0,0,0}; L[v] = 1.0;   // 顶点 v 的面积坐标
            double dNdL[6][3] = {{0}};
            dNdL[0][0]=4*L[0]-1; dNdL[1][1]=4*L[1]-1; dNdL[2][2]=4*L[2]-1;
            dNdL[3][0]=4*L[1]; dNdL[3][1]=4*L[0];
            dNdL[4][1]=4*L[2]; dNdL[4][2]=4*L[1];
            dNdL[5][0]=4*L[2]; dNdL[5][2]=4*L[0];
            double dux_dx=0,dux_dz=0,duz_dx=0,duz_dz=0;
            for (int a = 0; a < 6; a++) {
                const double nx = dNdL[a][0]*dLx[0]+dNdL[a][1]*dLx[1]+dNdL[a][2]*dLx[2];
                const double nz = dNdL[a][0]*dLz[0]+dNdL[a][1]*dLz[1]+dNdL[a][2]*dLz[2];
                const double ux = u[off_ux + gn[a]], uz = u[off_uz + gn[a]];
                dux_dx += nx*ux; dux_dz += nz*ux; duz_dx += nx*uz; duz_dz += nz*uz;
            }
            const double exx=dux_dx, ezz=duz_dz, gxz=dux_dz+duz_dx;
            const int nd = gn[v];
            sxx[nd] += A_e * (Cc*exx + lam*ezz);
            szz[nd] += A_e * (lam*exx + Cc*ezz);
            sxz[nd] += A_e * (mu*gxz);
            wsum[nd] += A_e;
        }
    }
#endif
    for (int k = 0; k < Nnodes; k++) {
        if (wsum[k] > 0.0) {
            const double inv = 1.0 / wsum[k];
            sxx[k] *= inv; szz[k] *= inv; sxz[k] *= inv;
        }
    }
}

} // namespace femasm
