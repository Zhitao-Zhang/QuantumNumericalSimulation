/* =============================================================================
 *  solver_bicgstab.cpp
 *  ---------------------------------------------------------------------------
 *  Jacobi 预条件 BiCGSTAB。算法与 Quantum/FDM/ImplicitV2/main.cpp 中的
 *  实现一致，方便 FEM 和 FDM 使用同一族迭代器做端到端对照。
 * ===========================================================================*/
#include "linear_solver.h"

#include <cmath>
#include <algorithm>

namespace femla {

SolveReport BiCGSTABSolver::solve(const SolveRequest &req)
{
    SolveReport rep;
    rep.solver_name = name();

    const SparseMatrixCSR &A = *req.A;
    const std::vector<double> &b = *req.b;
    std::vector<double> &x       = *req.x;

    const int n = A.n;
    if (!req.use_initial_guess) std::fill(x.begin(), x.end(), 0.0);

    std::vector<double> r(n), r0(n), p(n, 0.0), v(n, 0.0),
                        s(n), t(n), ph(n), sh(n);

    A.multiply(x, r);
    for (int i = 0; i < n; i++) r[i] = b[i] - r[i];
    r0 = r;

    double bnorm = 0.0;
    for (int i = 0; i < n; i++) bnorm += b[i] * b[i];
    bnorm = std::sqrt(bnorm);
    if (bnorm == 0.0) bnorm = 1.0;

    double rho_prev = 1.0, alpha = 1.0, omega = 1.0;

    for (int iter = 0; iter < req.max_iter; iter++) {
        double rho = 0.0;
        for (int i = 0; i < n; i++) rho += r0[i] * r[i];
        if (std::fabs(rho) < 1e-300) {
            r0 = r;
            rho = 0.0;
            for (int i = 0; i < n; i++) rho += r0[i] * r[i];
            rho_prev = 1.0; alpha = 1.0; omega = 1.0;
            std::fill(p.begin(), p.end(), 0.0);
            std::fill(v.begin(), v.end(), 0.0);
        }
        double beta = (rho / rho_prev) * (alpha / omega);
        for (int i = 0; i < n; i++)
            p[i] = r[i] + beta * (p[i] - omega * v[i]);

        for (int i = 0; i < n; i++) ph[i] = p[i] / A.diag[i];
        A.multiply(ph, v);

        double r0v = 0.0;
        for (int i = 0; i < n; i++) r0v += r0[i] * v[i];
        if (std::fabs(r0v) < 1e-300) {
            rep.iterations = iter + 1;
            rep.final_residual = 1.0;
            rep.success = false;
            return rep;
        }
        alpha = rho / r0v;

        for (int i = 0; i < n; i++) s[i] = r[i] - alpha * v[i];
        double snorm = 0.0;
        for (int i = 0; i < n; i++) snorm += s[i] * s[i];
        snorm = std::sqrt(snorm);
        if (snorm / bnorm < req.tol) {
            for (int i = 0; i < n; i++) x[i] += alpha * ph[i];
            rep.iterations     = iter + 1;
            rep.final_residual = snorm / bnorm;
            rep.success        = true;
            return rep;
        }

        for (int i = 0; i < n; i++) sh[i] = s[i] / A.diag[i];
        A.multiply(sh, t);
        double ts = 0.0, tt = 0.0;
        for (int i = 0; i < n; i++) { ts += t[i] * s[i]; tt += t[i] * t[i]; }
        if (tt < 1e-300) {
            rep.iterations = iter + 1;
            rep.final_residual = 1.0;
            rep.success = false;
            return rep;
        }
        omega = ts / tt;

        for (int i = 0; i < n; i++) x[i] += alpha * ph[i] + omega * sh[i];
        for (int i = 0; i < n; i++) r[i] = s[i] - omega * t[i];

        double rnorm = 0.0;
        for (int i = 0; i < n; i++) rnorm += r[i] * r[i];
        rnorm = std::sqrt(rnorm);
        if (rnorm / bnorm < req.tol) {
            rep.iterations     = iter + 1;
            rep.final_residual = rnorm / bnorm;
            rep.success        = true;
            return rep;
        }
        if (std::fabs(omega) < 1e-300) {
            rep.iterations     = iter + 1;
            rep.final_residual = rnorm / bnorm;
            rep.success        = false;
            return rep;
        }
        rho_prev = rho;
    }

    std::vector<double> tmp(n);
    A.multiply(x, tmp);
    double rn = 0.0;
    for (int i = 0; i < n; i++) { double d = tmp[i] - b[i]; rn += d * d; }
    rep.iterations     = req.max_iter;
    rep.final_residual = std::sqrt(rn) / bnorm;
    rep.success        = false;
    return rep;
}

SolveReport DummyCIMSolver::solve(const SolveRequest &req)
{
    // TODO: 未来在此接入 CIM/QUBO 后端。
    // 现阶段只是把控制流打到 BiCGSTAB，保证程序端到端可运行，
    // 但**求解调用点已经是 ILinearSolver::solve**，替换后无需改动上层。
    SolveReport rep = fallback_.solve(req);
    rep.solver_name = name();
    return rep;
}

std::unique_ptr<ILinearSolver> make_solver(SolverKind kind)
{
    switch (kind) {
        case SolverKind::BiCGSTAB : return std::unique_ptr<ILinearSolver>(new BiCGSTABSolver());
        case SolverKind::DummyCIM : return std::unique_ptr<ILinearSolver>(new DummyCIMSolver());
    }
    return std::unique_ptr<ILinearSolver>(new BiCGSTABSolver());
}

} // namespace femla
