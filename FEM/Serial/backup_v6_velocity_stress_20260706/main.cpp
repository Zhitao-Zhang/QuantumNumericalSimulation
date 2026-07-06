/* =============================================================================
 *  main.cpp  (2D velocity–stress elastic wave FEM  |  Serial CPU  |  P1 tri)
 *  ---------------------------------------------------------------------------
 *  第二版：加入 ADE-PML（8 个 ψ 辅助场作为额外主未知量）
 *
 *  层次架构（保持文档 CIM_FEM_Ax_b_guide_velocity_stress_v2.md 的分层）：
 *
 *      [PDE 层]        —— 二维速度–应力一阶弹性波 + stretched-coord PML
 *      [离散层]        —— P1 三角单元
 *                        主 5 方程用 Crank–Nicolson（默认，对齐 FDM/ImplicitV2）
 *                        -DUSE_BACKWARD_EULER 可退回 Backward Euler
 *                        8 个 ψ 辅助方程用 exp+trap 时间积分
 *                        13 * Nnodes DOF 组成的 A x = b
 *      [求解层]        —— ILinearSolver 抽象接口
 *                        当前用 Jacobi-BiCGSTAB；未来可替换为 CIM / QUBO
 *      [子域预留]      —— LinearSystem { A, b, x } 已显式暴露
 *
 *  参数（config.h）与 Quantum/FDM/ImplicitV2/main.cpp 完全对齐：
 *      - 网格   200 × 200，dx=dz=10 m，dt=1e-4 s，NT=4500
 *      - 声源   Ricker 30 Hz，位于节点 (isou, jsou)=(100, 60)
 *      - 接收器 节点 (irec, jrec)=(100, 110) 记录 σxx / σzz
 *      - 介质   固体背景 (Vp=5000, Vs=3300, ρ=2450)
 *              iz∈[100,120) 一条高速夹层 (Vp=6000, Vs=3400, ρ=2770)
 *      - PML    NP=30 cells，RC=1e-6，V0=Vp1
 * ===========================================================================*/
#include <cstdio>
#include <cstdlib>
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
    printf("==== 2D Velocity-Stress FEM + ADE-PML (P1 tri, Crank-Nicolson, serial) ====\n");
    printf("Grid : Nx=%d, Nz=%d,  dx=%g, dz=%g,  dt=%g,  NT=%d,  NP=%d\n",
           Nx, Nz, dx, dz, dt, Ntime, NP);
    printf("Blocks (13):  vx=%d vz=%d xx=%d zz=%d xz=%d\n",
           off_vx, off_vz, off_xx, off_zz, off_xz);
    printf("             pxSxx=%d pzSxz=%d pxSxz=%d pzSzz=%d\n",
           off_pxSxx, off_pzSxz, off_pxSxz, off_pzSzz);
    printf("             pxVx=%d  pzVz=%d  pxVz=%d  pzVx=%d   Ntot=%d\n",
           off_pxVx, off_pzVz, off_pxVz, off_pzVx, Ntot);
    printf("Source (isou, jsou) = (%d, %d)     Recv (irec, jrec) = (%d, %d)\n",
           isou, jsou, irec, jrec);

    // -------------------------------------------------------------------------
    // 1. 网格 + 介质
    // -------------------------------------------------------------------------
    femmesh::Mesh mesh;
    femmat::Material mat;
    mat.print_summary();

    // -------------------------------------------------------------------------
    // 2. 装配 A 与 M（含 PML 阻尼曲线 + 全部 ψ 耦合，时间无关）
    // -------------------------------------------------------------------------
    femasm::FemSystem sys;
    femasm::build_system(mesh, mat, dt, sys);

    // -------------------------------------------------------------------------
    // 3. 求解器（README §13：只依赖抽象接口）
    // -------------------------------------------------------------------------
    std::unique_ptr<femla::ILinearSolver> solver =
        femla::make_solver(femla::SolverKind::BiCGSTAB);
    printf("[solver] using %s   (tol=%.1e, max_iter=%d)\n",
           solver->name(), bicg_tol, bicg_max_it);

    // -------------------------------------------------------------------------
    // 4. 状态向量（13*Nnodes）、右端、warm-start 初值
    // -------------------------------------------------------------------------
    std::vector<double> q(Ntot, 0.0);
    std::vector<double> x(Ntot, 0.0);
    std::vector<double> b(Ntot, 0.0);

    // -------------------------------------------------------------------------
    // 5. 输出容器（与 FDM 保持相同的数组和文件名）
    // -------------------------------------------------------------------------
    std::vector<float> sou_record(Ntime, 0.0f);
    std::vector<float> trace_txx (Ntime, 0.0f);
    std::vector<float> trace_tzz (Ntime, 0.0f);
    std::vector<float> trace_p   (Ntime, 0.0f);
    std::vector<float> trace_zaxis_txx((size_t)Nz * (size_t)Ntime, 0.0f);

    const int rec_node = irec + jrec * Nx;
    const int src_node = isou + jsou * Nx;
    printf("[trace] rec_node=%d (dof_xx=%d)  src_node=%d (dof_xx=%d)\n",
           rec_node, off_xx + rec_node, src_node, off_xx + src_node);

    // -------------------------------------------------------------------------
    // 6. 时间循环 (Backward Euler + ADE-PML)
    //
    //       A * q^{n+1} = b(q^n, sou_t)
    //
    //   b 的形状：
    //       主 5 场   ： M q^n + [声源项到 σxx / σzz]
    //       8 个 ψ 场 ： E · ( M_L·ψ^n  -  (dt/2)·d·(G·field^n) )
    // -------------------------------------------------------------------------
    clock_t t0 = clock();
    for (int it = 0; it < Ntime; it++) {
        const double sou_t = ricker(it);
        sou_record[it] = (float)sou_t;

        // 6a. 组装右端
        femasm::assemble_rhs(sys, mesh, q, dt, sou_t, b);

        // 6b. warm-start：x0 ← q^n
        std::copy(q.begin(), q.end(), x.begin());

        // 6c. 求解 A x = b
        femla::SolveRequest req;
        req.A                 = &sys.A;
        req.b                 = &b;
        req.x                 = &x;
        req.use_initial_guess = true;
        req.max_iter          = bicg_max_it;
        req.tol               = bicg_tol;
        femla::SolveReport rep = solver->solve(req);

        // 6d. q ← x  (下一时步的旧解)
        std::swap(q, x);

        // 6e. 记录 trace（只关心主 5 场）
        trace_txx[it] = (float)q[off_xx + rec_node];
        trace_tzz[it] = (float)q[off_zz + rec_node];
        trace_p  [it] = 0.5f * (trace_txx[it] + trace_tzz[it]);
        for (int iz = 0; iz < Nz; iz++)
            trace_zaxis_txx[(size_t)iz * (size_t)Ntime + (size_t)it]
                = (float)q[off_xx + (isou + iz * Nx)];

        if (it % 20 == 0) {
            printf("it=%4d  t=%.4f s  [%s] iters=%4d  res=%.3e  sou=%.3e  txx(rec)=%.3e\n",
                   it, it * dt, rep.solver_name, rep.iterations, rep.final_residual,
                   sou_t, trace_txx[it]);
        }

        if (it % snap_stride == 0) {
            femio::write_snapshot(it, q, Nx, Nz,
                                  off_vx, off_vz, off_xx, off_zz, off_xz);
        }
    }
    printf("Time loop done. Elapsed %.2f s\n",
           (double)(clock() - t0) / CLOCKS_PER_SEC);

    // -------------------------------------------------------------------------
    // 7. Trace 落盘（与 FDM 同名）
    // -------------------------------------------------------------------------
    femio::wfile2d("source_time.dat",    sou_record.data(), 1, Ntime);
    femio::wfile2d("trace_txx.dat",      trace_txx.data(),  1, Ntime);
    femio::wfile2d("trace_tzz.dat",      trace_tzz.data(),  1, Ntime);
    femio::wfile2d("trace_pressure.dat", trace_p.data(),    1, Ntime);
    femio::wfile2d("trace_zaxis_txx.dat", trace_zaxis_txx.data(), Nz, Ntime);
    printf("Saved trace_zaxis_txx.dat  shape = Nz x NT = %d x %d  (col x_src = %d)\n",
           Nz, Ntime, isou);

    printf("Done. snap_<field>_itXXXX.dat ; trace_*.dat ; trace_zaxis_txx.dat ; source_time.dat\n");
    return 0;
}
