/* =============================================================================
 *  main.cpp   —  2D 位移二阶弹性波 FEM（P1 三角，Newmark 隐式，串行）
 *  ---------------------------------------------------------------------------
 *  方案 B（SESSION 7）：丢弃 velocity-stress collocated 格式（波前菱形，已用
 *  diag_*.py 确诊为 collocation 各向异性），改成标准位移二阶格式：
 *
 *      M ü + C u̇ + K u = f
 *
 *  Newmark 平均加速度（β=1/4, γ=1/2，无条件稳定）每步解：
 *
 *      A a^{n+1} = f^{n+1} - C v_pred - K u_pred,   A = M + γ·dt·C + β·dt²·K
 *      u^{n+1} = u_pred + β·dt²·a^{n+1}
 *      v^{n+1} = v_pred + γ·dt·a^{n+1}
 *      u_pred  = u^n + dt·v^n + dt²(½-β)·a^n
 *      v_pred  = v^n + dt(1-γ)·a^n
 *
 *  求解层仍走 ILinearSolver / BiCGSTAB（CIM 落位口不变）。
 *  参数与 Quantum/FDM/ImplicitV2 对齐；输出 vx/vz/σxx/σzz/σxz 同名同格式对照。
 * ===========================================================================*/
#include <cstdio>
#include <ctime>
#include <vector>
#include <memory>
#include <algorithm>

#include "config.h"
#include "mesh.h"
#include "material.h"
#include "sparse_matrix.h"
#include "linear_solver.h"
#include "assembler.h"
#include "field_io.h"

using namespace femcfg;

int main(void)
{
    printf("==== 2D Displacement (2nd-order) Elastic FEM  |  P1 tri  |  Newmark  |  serial ====\n");
    printf("Grid : Nx=%d Nz=%d  dx=%g dz=%g  dt=%g  NT=%d  NABS=%d\n",
           Nx, Nz, dx, dz, dt, Ntime, NABS);
    printf("DOF  : 2 fields/node (ux,uz)  off_ux=%d off_uz=%d  Ntot=%d\n",
           off_ux, off_uz, Ntot);
    printf("Newmark: beta=%g gamma=%g   (avg-accel, unconditionally stable)\n",
           nm_beta, nm_gamma);
    printf("Source (isou,jsou)=(%d,%d)  explosive moment  Recv (irec,jrec)=(%d,%d)\n",
           isou, jsou, irec, jrec);

    // 1. 网格 + 介质
    femmesh::Mesh mesh;
    femmat::Material mat;
    mat.print_summary();

    // 2. 装配 K/M/C + Newmark 有效矩阵 A（时间无关）
    femasm::ElastSystem sys;
    sys.dt = dt; sys.beta = nm_beta; sys.gamma = nm_gamma;
    femasm::build_system(mesh, mat, dt, sys);

    // 3. 求解器
    std::unique_ptr<femla::ILinearSolver> solver =
        femla::make_solver(femla::SolverKind::BiCGSTAB);
    printf("[solver] %s  (tol=%.1e, max_iter=%d)\n",
           solver->name(), bicg_tol, bicg_max_it);

    // 4. 状态向量
    std::vector<double> u(Ntot, 0.0), v(Ntot, 0.0), a(Ntot, 0.0);
    std::vector<double> u_pred(Ntot, 0.0), v_pred(Ntot, 0.0);
    std::vector<double> F(Ntot, 0.0), Ku(Ntot, 0.0), Cv(Ntot, 0.0);
    std::vector<double> b(Ntot, 0.0), x(Ntot, 0.0);
    // PML 记忆变量（显式更新，不进 A）
    std::vector<double> psi_xx(Ntot, 0.0), psi_zz(Ntot, 0.0), Fpsi(Ntot, 0.0);

    // 应力恢复缓冲 + 快照用的速度分量
    std::vector<double> sxx, szz, sxz;
    std::vector<double> vx_snap(Nnodes), vz_snap(Nnodes);

    // 5. 输出容器
    std::vector<float> sou_record(Ntime, 0.0f);
    std::vector<float> trace_txx (Ntime, 0.0f);
    std::vector<float> trace_tzz (Ntime, 0.0f);
    std::vector<float> trace_p   (Ntime, 0.0f);
    std::vector<float> trace_zaxis_txx((size_t)Nz * (size_t)Ntime, 0.0f);

    // 接收/源节点：顶点 (ix,iz) 映射到细网格全局 id（P1 时等于 ix+iz*Nx）
    const int rec_node = node_id(irec, jrec);
    printf("[trace] rec_node=%d  src_node=%d\n", rec_node, node_id(isou, jsou));

    const double beta = nm_beta, gamma = nm_gamma;
    const double bdt2 = beta * dt * dt;
    const double gdt  = gamma * dt;

    // 6. 时间循环（Newmark 隐式）
    clock_t t0 = clock();
    for (int it = 0; it < Ntime; it++) {
        const double sou_t = ricker(it);
        sou_record[it] = (float)sou_t;

        // 6a. 声源等效力 F^{n+1}
        femasm::assemble_source_force(sys, sou_t, F);

        // 6b. 预测量
        const double c_upred_a = dt * dt * (0.5 - beta);
        const double c_vpred_a = dt * (1.0 - gamma);
        for (int i = 0; i < Ntot; i++) {
            u_pred[i] = u[i] + dt * v[i] + c_upred_a * a[i];
            v_pred[i] = v[i] + c_vpred_a * a[i];
        }

        // 6b'. PML 记忆变量显式更新，得到 F_ψ（NABS=0 时恒 0）
        femasm::pml_update_force(sys, u_pred, psi_xx, psi_zz, Fpsi);

        // 6c. RHS = F - C·v_pred - K_eff·u_pred - F_ψ
        //     （A a = b 里 A、K_eff、C 都是常系数对称阵；F_ψ 是已知显式修正）
        sys.K_eff.multiply(u_pred, Ku);
        sys.C.multiply(v_pred, Cv);
        for (int i = 0; i < Ntot; i++)
            b[i] = F[i] - Cv[i] - Ku[i] - Fpsi[i];

        // 6d. 解 A·a^{n+1} = b（warm-start：a^{n+1} ← a^n）
        std::copy(a.begin(), a.end(), x.begin());
        femla::SolveRequest req;
        req.A = &sys.A; req.b = &b; req.x = &x;
        req.use_initial_guess = true;
        req.max_iter = bicg_max_it; req.tol = bicg_tol;
        femla::SolveReport rep = solver->solve(req);

        // 6e. 校正 u, v, a
        for (int i = 0; i < Ntot; i++) {
            const double an = x[i];
            u[i] = u_pred[i] + bdt2 * an;
            v[i] = v_pred[i] + gdt  * an;
            a[i] = an;
        }

        // 6f. 恢复应力（用于 trace 与快照，与 FDM σ 对照）
        femasm::recover_nodal_stress(mesh, mat, u, sxx, szz, sxz);
        trace_txx[it] = (float)sxx[rec_node];
        trace_tzz[it] = (float)szz[rec_node];
        trace_p  [it] = 0.5f * (trace_txx[it] + trace_tzz[it]);
        for (int iz = 0; iz < Nz; iz++)
            trace_zaxis_txx[(size_t)iz * (size_t)Ntime + (size_t)it]
                = (float)sxx[node_id(isou, iz)];   // 源列，顶点→细网格映射

        if (it % 20 == 0)
            printf("it=%4d t=%.4f [%s] iters=%4d res=%.2e  sou=%.2e txx(rec)=%.3e\n",
                   it, it * dt, rep.solver_name, rep.iterations, rep.final_residual,
                   sou_t, trace_txx[it]);

        if (it % snap_stride == 0) {
            for (int k = 0; k < Nnodes; k++) {
                vx_snap[k] = v[off_ux + k];
                vz_snap[k] = v[off_uz + k];
            }
            femio::write_snapshot_fields(it, Nx, Nz, vx_snap, vz_snap, sxx, szz, sxz);
        }
    }
    printf("Time loop done. Elapsed %.2f s\n",
           (double)(clock() - t0) / CLOCKS_PER_SEC);

    // 7. Trace 落盘（与 FDM 同名）
    femio::wfile2d("source_time.dat",     sou_record.data(), 1, Ntime);
    femio::wfile2d("trace_txx.dat",       trace_txx.data(),  1, Ntime);
    femio::wfile2d("trace_tzz.dat",       trace_tzz.data(),  1, Ntime);
    femio::wfile2d("trace_pressure.dat",  trace_p.data(),    1, Ntime);
    femio::wfile2d("trace_zaxis_txx.dat", trace_zaxis_txx.data(), Nz, Ntime);
    printf("Done. snap_<field>_itXXXX.dat ; trace_*.dat ; source_time.dat\n");
    return 0;
}
