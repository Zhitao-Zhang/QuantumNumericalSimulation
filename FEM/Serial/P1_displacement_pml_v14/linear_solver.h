/* =============================================================================
 *  linear_solver.h
 *  ---------------------------------------------------------------------------
 *  README §13 要求的抽象求解器接口。
 *
 *  设计要点：
 *      - 时间推进 / PDE 装配层 **不** 直接调用 BiCGSTAB。
 *      - 只依赖 ILinearSolver::solve(...) 接口，任何后端（BiCGSTAB / GMRES /
 *        CIM / QUBO / 子域块迭代）都只需实现同一接口。
 *      - SolveRequest 里持有 A / b / x 的 const 指针 + is_initial_guess 标志，
 *        允许 warm-start；SolveReport 返回收敛信息。
 *
 *  README §11 预留：Subdomain 数据结构与块解求器（当前先不实现，
 *      但接口签名保持一致，未来 BlockJacobiSolver::solve(...) 可直接落位）。
 * ===========================================================================*/
#ifndef FEM_LINEAR_SOLVER_H
#define FEM_LINEAR_SOLVER_H

#include <memory>
#include <vector>

#include "sparse_matrix.h"

namespace femla {

struct SolveRequest {
    const SparseMatrixCSR* A            = nullptr;
    const std::vector<double>* b        = nullptr;
    std::vector<double>*       x        = nullptr;
    bool use_initial_guess              = false;

    int    max_iter                     = 5000;
    double tol                          = 1.0e-9;
};

struct SolveReport {
    bool   success       = false;
    int    iterations    = 0;
    double final_residual= 0.0;
    const char* solver_name = "none";
};

class ILinearSolver {
public:
    virtual ~ILinearSolver() = default;
    virtual SolveReport solve(const SolveRequest &req) = 0;
    virtual const char* name() const = 0;
};

// -----------------------------------------------------------------------------
// 具体后端 1：Jacobi 预条件 BiCGSTAB（与 FDM 中的实现一致）
// -----------------------------------------------------------------------------
class BiCGSTABSolver : public ILinearSolver {
public:
    SolveReport solve(const SolveRequest &req) override;
    const char* name() const override { return "BiCGSTAB(Jacobi)"; }
};

// -----------------------------------------------------------------------------
// 具体后端 2：DummyCIMSolver（当前占位，后续 CIM/QUBO 落位于此）
//   现阶段直接调用 BiCGSTAB 完成求解，以便端到端跑通；接口和调用点已到位。
// -----------------------------------------------------------------------------
class DummyCIMSolver : public ILinearSolver {
public:
    SolveReport solve(const SolveRequest &req) override;
    const char* name() const override { return "DummyCIM(fallback=BiCGSTAB)"; }

private:
    BiCGSTABSolver fallback_;
};

// 工厂：按名字或枚举取一个具体求解器实现
enum class SolverKind {
    BiCGSTAB,
    DummyCIM,
};
std::unique_ptr<ILinearSolver> make_solver(SolverKind kind);

} // namespace femla

#endif // FEM_LINEAR_SOLVER_H
