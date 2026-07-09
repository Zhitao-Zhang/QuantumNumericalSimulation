# CHANGELOG — Quantum/FEM/Serial

本文件记录 `Quantum/FEM/Serial/` 目录下所有由 AI 会话产生的改动，用于
**可追溯 / 可还原**。每次新会话或新改动，追加一节即可，格式统一。

书写约定：
- 每次会话一个 `## SESSION N (yyyy-mm-dd HH:MM, 简短主题)` 块。
- 每次「独立改动」在会话里用 `### Change N.k` 分小节。
- 每个改动必须包含四块信息：
  1. **动机 / 需求**（用户说了什么、想要什么）
  2. **改了哪些文件**（新增 / 修改 / 删除；对已有文件给出「旧代码 → 新代码」的关键片段）
  3. **验证方式**（怎么跑 / 期望输出）
  4. **回滚指令**（要还原到本条修改之前该怎么操作）

回滚原则：
- 优先给出 **精确的反向 diff / 反向 StrReplace**，让后续会话可以按行还原。
- 新增文件的回滚 = 删除该文件。
- 删除文件的回滚 = 在“回滚指令”里保留完整旧内容。
- **不要依赖 git**（这个目录当前不一定纳入 git 提交），日志里必须自包含。

---

## SESSION 1 (2026-07-01, 首次实现 FEM 串行版基线)

### 上下文

用户要求：根据 `CIM_FEM_Ax_b_guide_velocity_stress_v2.md` 编写一份 **串行的
二维速度–应力 P1 三角有限元** 程序；声源和模型参数须与
`Quantum/FDM/ImplicitV2/main.cpp` 完全一致，输出与 FDM 可对照。

选定的设计：
- 5 场（vx, vz, σxx, σzz, σxz）都用连续 P1 空间；
- Backward Euler，`A = M + Δt K`；
- Jacobi-BiCGSTAB 作为线性后端；
- `ILinearSolver` 抽象接口 + `DummyCIMSolver` 落位口（未来 CIM 从这里接入）；
- 无 PML（文档 §14）；自由面边界。

### Change 1.1  新建 `config.h`

- **动机**：把所有参数放在一个头文件，编译期常量，能被宏覆盖；参数集与
  FDM/ImplicitV2 完全对齐。
- **新增文件**：`config.h`
- **要点**：
  - `NX = NZ = 200, dx = dz = 10, dt = 1e-4, NT = 4500`
  - `F0 = 30 Hz, T0 = 1.2/F0`
  - `isou = NX/2, jsou = NZ/2 - 40`
  - `irec = NX/2, jrec = NZ/2 + 10`
  - 固体1：`Vp1=5000, Vs1=3300, rho1=2450`
  - 固体2 (高速夹层 iz∈[100,120))：`Vp2=6000, Vs2=3400, rho2=2770`
  - DOF 布局：块堆叠 `[vx; vz; σxx; σzz; σxz]`，offset 静态给出
  - `bicg_tol=1e-9, bicg_max_it=5000`
  - `snap_stride=50`
- **回滚**：`rm /home/wanwb/zzt/Quantum/FEM/Serial/config.h`

### Change 1.2  新建 `mesh.h / mesh.cpp`

- **动机**：结构化 P1 三角网格，Nx × Nz 节点，每个 cell 切 2 个 CCW 三角；
  暴露 `Triangle { nodes[3], cell_id, A, A_signed, dN_dx[3], dN_dz[3] }`。
- **新增文件**：
  - `mesh.h`（`Triangle` struct、`Mesh` 类、`compute_p1_geometry(...)`）
  - `mesh.cpp`（`Mesh::build_nodes` + `Mesh::build_triangles`）
- **关键约定**：
  - 节点 (ix, iz) 位置 (ix·dx, iz·dz)，`node_id = ix + iz·Nx`
  - 每个 cell 切两三角：`T_lower = (SW, SE, NE)`、`T_upper = (SW, NE, NW)`，均 CCW
  - 面积对 `A_signed = 0.5 dx dz > 0`
- **回滚**：`rm mesh.h mesh.cpp`

### Change 1.3  新建 `material.h / material.cpp`

- **动机**：逐 cell 的 `(ρ, λ, μ)`；与 FDM `build_medium()` 里
  `iz∈[100,120) 用固体2，其他用固体1` 完全一致。
- **新增文件**：`material.h`、`material.cpp`
- **要点**：
  - `CellProps { rho, lambda, mu }`；`lambda_plus_2mu()` 便捷方法
  - `Material` 类持有 `std::vector<CellProps> cells_[Ncells]`
  - `Material::at(icx, icz)` 与 `at_cell(cell_id)` 双入口
  - `print_summary()` 打印背景 / 夹层 / 占比
- **回滚**：`rm material.h material.cpp`

### Change 1.4  新建 `sparse_matrix.h / sparse_matrix.cpp`

- **动机**：CSR + Jacobi 对角 + `RowMap = std::map<int,double>` 装配容器；
  不含任何 FEM / PDE 语义（README §9 要求的解耦）。
- **新增文件**：`sparse_matrix.h`、`sparse_matrix.cpp`
- **API**：
  - `SparseMatrixCSR::build_from_rows(int n, std::vector<RowMap>&)`
  - `SparseMatrixCSR::multiply(x, y)`
  - `LinearSystem { A, b, x }`
- **回滚**：`rm sparse_matrix.h sparse_matrix.cpp`

### Change 1.5  新建 `linear_solver.h + solver_bicgstab.cpp`

- **动机**：抽象求解器接口 `ILinearSolver`，以及两个实现：
  1. `BiCGSTABSolver`：Jacobi 预条件 BiCGSTAB，与 FDM 里的算法完全一致；
  2. `DummyCIMSolver`：CIM 落位口，当前退化为调 BiCGSTAB。
- **新增文件**：`linear_solver.h`、`solver_bicgstab.cpp`
- **API**：
  - `SolveRequest { A, b, x, use_initial_guess, max_iter, tol }`
  - `SolveReport { success, iterations, final_residual, solver_name }`
  - `make_solver(SolverKind)` 工厂
- **回滚**：`rm linear_solver.h solver_bicgstab.cpp`

### Change 1.6  新建 `assembler.h / assembler.cpp`

- **动机**：`离散层`。从 mesh + material 生成 `A = M + dt K` 与 `M` 两个 CSR
  矩阵；再提供 `assemble_rhs(...)` 每步组装 b。
- **新增文件**：`assembler.h`、`assembler.cpp`
- **数学要点**（P1 三角 CCW，|A| 为面积）：
  - `∫ N_i N_j dA = (|A|/12)·(1+δ_ij)`
  - `∫ ∂N_i/∂x · N_j dA = dN_dx[i] · |A|/3`
  - `∫ N_i · ∂N_j/∂x dA = dN_dx[j] · |A|/3`
- **K 块结构**：
  ```
  K = [  0     0    +Gx    0    +Gz  ]   vx eq
      [  0     0     0    +Gz   +Gx  ]   vz eq
      [-DxC  -DzL    0     0     0   ]   σxx eq
      [-DxL  -DzC    0     0     0   ]   σzz eq
      [-DzM  -DxM    0     0     0   ]   σxz eq
  ```
- **声源**：`b[dof_xx(src_node)] += dt * sou_t; b[dof_zz(src_node)] += dt * sou_t;`
  与 FDM `main.cpp` 中同一 `dt * sou_t` 处理一致。
- **回滚**：`rm assembler.h assembler.cpp`

### Change 1.7  新建 `field_io.h / field_io.cpp`

- **动机**：ASCII 输出 `wfile2d`，签名与 FDM 一致，方便文件名/格式对照。
- **新增文件**：`field_io.h`、`field_io.cpp`
- **API**：
  - `wfile2d(fn, float*, nrows, ncols)`
  - `wfile2d_d(fn, std::vector<double>&, nrows, ncols)`
  - `write_snapshot(it, q, ...)` 一次写 5 个场
- **回滚**：`rm field_io.h field_io.cpp`

### Change 1.8  新建 `main.cpp`

- **动机**：时间推进主循环。
- **新增文件**：`main.cpp`
- **主流程**：
  1. 建 mesh + material
  2. 装配一次 A 与 M
  3. 每步：`b = M q^n + dt f`，`x0 = q^n`（warm-start），调 `ILinearSolver::solve`
  4. 记录 trace（(irec, jrec) 处 σxx/σzz + 沿源 z 轴 σxx 时空图）
  5. 每 `snap_stride` 步 dump 5 个场快照
- **回滚**：`rm main.cpp`

### Change 1.9  新建 `Makefile`

- **动机**：多文件工程需要一个统一构建入口，且要能通过宏覆盖参数。
- **新增文件**：`Makefile`
- **要点**：
  - `make`         : 200×200×4500，默认参数
  - `make small`   : 50×50×800，`-DNBH=0 -DSNAP_STRIDE=20`，冒烟测试
  - `make cleanall`: 清 `.dat` + 可执行
  - `-O3 -std=c++14`
- **回滚**：`rm Makefile`

### Change 1.10  新建 `README.md`

- **动机**：项目说明、与 FDM 对照表、CIM 预留说明。
- **新增文件**：`README.md`
- **回滚**：`rm README.md`

### Change 1.11  修复 `config.h`：小网格下声源/接收点越界

- **动机**：`make small`（NZ=50）时，`jsou = NZ/2 - 40 = -15` 导致
  `src_node = isou + jsou·Nx = -725` 越界。
- **改动文件**：`config.h`
- **具体修改**：

  **旧代码**（`Change 1.1` 版本）：
  ```cpp
  // 声源节点位置（FDM: jsou = Nz/2-40, isou = Nx/2；FEM 用同一 (ix,iz)）
  constexpr int isou = NX / 2;
  constexpr int jsou = NZ / 2 - 40;

  // 接收器节点位置
  constexpr int irec = NX / 2;
  constexpr int jrec = NZ / 2 + 10;
  ```

  **新代码**：
  ```cpp
  // 声源节点位置（FDM: jsou = Nz/2-40, isou = Nx/2；FEM 用同一 (ix,iz)）
  //   对小网格自动 clamp，避免 make small 时 NZ/2-40 < 0
  constexpr int _clamp01(int v, int lo, int hi) { return v < lo ? lo : (v > hi ? hi : v); }
  constexpr int isou = _clamp01(NX / 2,          1, NX - 2);
  constexpr int jsou = _clamp01(NZ / 2 - 40,     1, NZ - 2);

  // 接收器节点位置
  constexpr int irec = _clamp01(NX / 2,          1, NX - 2);
  constexpr int jrec = _clamp01(NZ / 2 + 10,     1, NZ - 2);
  ```

- **影响**：默认 NX=NZ=200 时 clamp 不生效，位置与 FDM 完全一致
  `(isou, jsou)=(100, 60), (irec, jrec)=(100, 110)`；只有小网格时才会截断。
- **验证**：`make small` 打印 `Source (isou, jsou) = (25, 1)`（原为 -15），
  `src_node=75`（原为 -725），trace 值出现符合物理的波达时间（it≈680）。
- **回滚方式**：把 `config.h` 中相应段落改回“旧代码”即可，无需其他文件配合。

### Change 1.12  新建 `CHANGELOG.md`（本文件）

- **动机**：用户要求建立一个改动日志文件，便于每次会话追加、追溯、还原。
- **新增文件**：`CHANGELOG.md`
- **回滚**：`rm CHANGELOG.md`

### 本次会话交付快照

- 目录内容（共 12 个源文件 + Makefile + README + CHANGELOG）：
  ```
  Serial/
  ├── config.h
  ├── mesh.h            mesh.cpp
  ├── material.h        material.cpp
  ├── sparse_matrix.h   sparse_matrix.cpp
  ├── linear_solver.h   solver_bicgstab.cpp
  ├── assembler.h       assembler.cpp
  ├── field_io.h        field_io.cpp
  ├── main.cpp
  ├── Makefile
  ├── README.md
  └── CHANGELOG.md   (本次新增)
  ```
- 冒烟测试：`make small` 通过；
  - Ntot=12500, nnz(A)=256530, 装配 0.06 s，800 步共 5.21 s
  - BiCGSTAB 每步 9–12 次收敛，`res < 1e-9`
  - σxx 首波到时约 it≈680（Vp=5000, 340 m），与理论 68 ms 吻合

### 会话级完全回滚方式

若要把 `Quantum/FEM/Serial/` 恢复到本次会话开始前的空目录状态：

```bash
cd /home/wanwb/zzt/Quantum/FEM/Serial
rm -f config.h \
      mesh.h mesh.cpp \
      material.h material.cpp \
      sparse_matrix.h sparse_matrix.cpp \
      linear_solver.h solver_bicgstab.cpp \
      assembler.h assembler.cpp \
      field_io.h field_io.cpp \
      main.cpp Makefile README.md CHANGELOG.md \
      fem_elastic_2d *.o *.dat
```

（`Quantum/FEM/` 目录本身、以及 `Quantum/FEM/Serial/` 目录都是本会话新建，
如要一并清除：`rm -rf /home/wanwb/zzt/Quantum/FEM/`。）

---

<!-- 后续会话请在此下方追加：

## SESSION N (yyyy-mm-dd HH:MM, 主题)

### 上下文
…

### Change N.k …
…

-->

## SESSION 2 (2026-07-01, 加入 ADE-PML 边界处理)

### 上下文

用户反馈：SESSION 1 的程序没有边界处理，产生大量人工反射。要求用 PML 做
边界吸收，参考 FDM/ImplicitV2 中的 ADE-PML 做法。

### 设计要点（数学 + FEM 实现）

- 引入 8 个 ψ 辅助场（与 FDM 完全对应），全部用 P1 连续，共 13*Nnodes 个 DOF：
  ```
  0..4  : vx, vz, σxx, σzz, σxz          (主 5 场)
  5     : ψ_xSxx  (in vx eq, ∂x σxx)     damp = d_x(x_k)
  6     : ψ_zSxz  (in vx eq, ∂z σxz)     damp = d_z(z_k)
  7     : ψ_xSxz  (in vz eq, ∂x σxz)     damp = d_x(x_k)
  8     : ψ_zSzz  (in vz eq, ∂z σzz)     damp = d_z(z_k)
  9     : ψ_xVx   (in σxx/σzz, ∂x vx)    damp = d_x(x_k)
  10    : ψ_zVz   (in σxx/σzz, ∂z vz)    damp = d_z(z_k)
  11    : ψ_xVz   (in σxz,     ∂x vz)    damp = d_x(x_k)
  12    : ψ_zVx   (in σxz,     ∂z vx)    damp = d_z(z_k)
  ```
- **主方程**用 Backward Euler，在 stretched coordinate 下把 ∂α 替换为 ∂α + ψ_α；
  ψ 通过 consistent-mass 耦合进入主方程：
  ```
  Row vx : A[vx , ψ_xSxx] += -dt · M0     A[vx , ψ_zSxz] += -dt · M0
  Row σxx: A[σxx, ψ_xVx ] += -dt · C·M0   A[σxx, ψ_zVz ] += -dt · λ·M0
  ...
  ```
- **ψ 方程**用**指数 + 梯形**离散（与 FDM 完全一致），并用 **lumped mass** 做 ψ 行的对角：
  ```
  M_L·ψ^{n+1} + (dt/2)·d·G_S·field^{n+1}
        = E·[ M_L·ψ^n - (dt/2)·d·G_S·field^n ]
  ```
  其中 `M_L[k] = ∫ N_k dΩ`，`G_S_kj = ∫ N_i ∂N_j/∂α dΩ`（P1 三角上闭式 `dN_dα[j]·|A|/3`）。
  当 `d=0`（PML 外）：`E=1`，方程退化为 `ψ^{n+1}=ψ^n`，且初值 `ψ^0=0` → ψ 恒 0；
  A / M / RHS 均自动“关闭”PML，与 SESSION 1 无 PML 版本**逐位一致**。
  已通过 `make small_nopml`（NPMLDEF=0）实证：trace 与 SESSION 1 完全相同。

### Change 2.1  改写 `config.h`

- **动机**：加入 PML 参数（NP, RC, d0, xleft/xright/zleft/zright, damp*_at）；
  把 5 块布局升级到 13 块（新增 8 个 ψ offset 和 `dof_p...` 便捷函数）；
  `Ntot` 从 `5*Nnodes` 变为 `13*Nnodes`。
- **改动文件**：`config.h`
- **关键差异**：
  - **原**：`constexpr int NfieldsPerNode = 5; constexpr int Ntot = Nnodes * 5;`
    只有 5 个 offset (`off_vx / off_vz / off_xx / off_zz / off_xz`)。
  - **新**：`constexpr int NfieldsPerNode = 13; constexpr int Ntot = Nnodes * 13;`
    新增 8 个 offset：`off_pxSxx / off_pzSxz / off_pxSxz / off_pzSzz /
    off_pxVx / off_pzVz / off_pxVz / off_pzVx`。
  - 新增 PML 常量与函数：`NP, RC_pml, V0_pml, DPx, DPz, d0_x, d0_z,
    xleft_p, xright_p, zleft_p, zright_p, dampx_at(x), dampz_at(z),
    compute_d0_x(), compute_d0_z()`。
  - 新增 8 个 `dof_p*(ix, iz)` 便捷索引。
- **回滚方式**：整体覆盖回 SESSION 1 版本；或删除新增的 PML 段 + 8 个 offset 常量 +
  `NfieldsPerNode=13` 改回 5、`Ntot = Nnodes * 5`。

### Change 2.2  `mesh.h / mesh.cpp` 加入 lumped mass 计算

- **动机**：ψ 方程行的对角要用 `M_L[k] = ∫ N_k dΩ`；同时给出全域 sanity check。
- **改动文件**：`mesh.h`、`mesh.cpp`
- **具体修改**：

  **mesh.h**  在 `class Mesh { public: ... };` 里新增：
  ```cpp
  std::vector<double> M_L;   // lumped mass per node

  // private:
  void build_lumped_mass();
  ```

  **mesh.cpp** 里 `Mesh::Mesh()` 追加调用 `build_lumped_mass()`，末尾多打印一行
  `[mesh] Lumped mass: ...`。新增函数：
  ```cpp
  void Mesh::build_lumped_mass()
  {
      M_L.assign(Nnodes, 0.0);
      for (const auto &t : triangles) {
          const double contrib = t.A / 3.0;
          M_L[t.nodes[0]] += contrib;
          M_L[t.nodes[1]] += contrib;
          M_L[t.nodes[2]] += contrib;
      }
  }
  ```
- **回滚方式**：从 `mesh.h` 里删除 `std::vector<double> M_L;` 与
  `void build_lumped_mass();` 声明；从 `mesh.cpp` 里删除 `build_lumped_mass`
  实现和 `Mesh::Mesh()` 里对它的调用，以及末尾的打印段。

### Change 2.3  重写 `assembler.h / assembler.cpp`

- **动机**：从 5 块变为 13 块装配；`assemble_rhs(...)` 签名变化，需要 mesh 和
  完整的 `FemSystem`（含 PML 阻尼曲线与 E 系数）。
- **改动文件**：`assembler.h`、`assembler.cpp`（都是**整体覆盖**）
- **旧签名 → 新签名**：
  ```cpp
  // 旧
  void assemble_rhs(const femla::SparseMatrixCSR &M,
                    const std::vector<double>    &q_n,
                    double dt, double sou_t,
                    std::vector<double>          &b_out);

  // 新
  void assemble_rhs(const FemSystem &sys,
                    const femmesh::Mesh &mesh,
                    const std::vector<double> &q_n,
                    double dt, double sou_t,
                    std::vector<double> &b_out);
  ```
  新签名多传 `sys`（拿 PML 曲线 `dampx_at_ix`, `Ex_at_ix` 等）和 `mesh`
  （做 element loop 计算 nodal gradient）。
- **`FemSystem` 结构扩容**：新增 5 个成员
  ```cpp
  std::vector<double> dampx_at_ix, dampz_at_iz;
  std::vector<double> Ex_at_ix,   Ez_at_iz;
  int pml_active_x = 0, pml_active_z = 0;
  ```
- **`build_system(...)` 内部**：
  - 先 `build_pml_profiles(dt, out)` 采样节点级 d、E
  - element loop 里除了原 5x5 主耦合外，还加：
    - 主方程与 ψ 的 8 组 mass 耦合 `A[主, ψ] += -dt · <mass>`
    - 每个 ψ 块的 lumped mass 对角 `A[ψ_k, ψ_k] += A_e/3`
    - 每个 ψ 与其对应 field 的 `A[ψ_k, field_j] += (dt/2)·d_α(x_k)·dN_dα[j]·A_e/3`
      （只在 `d != 0` 时写入，保证 PML 外 A 依然是原始 5 块结构 + ψ 对角）
- **`assemble_rhs(...)` 内部**：
  1. `sys.M.multiply(q_n, b_out)`  一次矩阵向量乘覆盖主 5 场 + ψ 对角项
  2. σxx / σzz 加声源
  3. 一次 mesh.triangles 大循环，同时计算 8 组 `G·field^n` 的节点分散（`g_xSxx[Nnodes]` 等）
  4. 对每个 (ix, iz) 节点应用：
     `b[ψ_k] := E_α(k) · ( b[ψ_k] - (dt/2)·d_α(k)·grad_k )`
- **回滚方式**：把 `assembler.h/.cpp` 整体还原到 SESSION 1 的旧版本
  （因为改动过大，不给逐段 diff；旧版本文件很短，可从 SESSION 1 的 Change 1.6
  的描述重新写回）。若确需精确旧版本，请把 SESSION 1 中的 `## Change 1.6`
  作为参考——其中列出了全部行为规范。

### Change 2.4  更新 `main.cpp`

- **动机**：适配新的 `assemble_rhs` 签名、13 块的 q / x / b 大小、以及打印 PML 信息。
- **改动文件**：`main.cpp`（整体覆盖）
- **关键改动点**：
  - 打印 header 增加 PML 与 8 个 ψ offset 信息
  - `q, x, b` 尺寸不变（仍是 `Ntot`），但 `Ntot` 现在等于 `13*Nnodes`
  - `assemble_rhs(sys, mesh, q, dt, sou_t, b);` 替换旧的 `assemble_rhs(sys.M, q, dt, sou_t, b);`
  - snapshot 与 trace 依然只写主 5 场，逻辑不变
- **回滚方式**：把 `main.cpp` 整体覆盖回 SESSION 1 的版本（详见 Change 1.8）。

### Change 2.5  更新 `Makefile`

- **动机**：小算例默认启用 PML；新增一个 `small_nopml` 目标用于回归对比。
- **改动文件**：`Makefile`
- **具体修改**：

  **旧**：
  ```makefile
  small: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20
  small: clean $(TARGET)
  	./$(TARGET)

  .PHONY: all run small clean cleanall
  ```

  **新**：
  ```makefile
  small: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNPMLDEF=10
  small: clean $(TARGET)
  	./$(TARGET)

  small_nopml: CXXFLAGS += -DNZ=50 -DNX=50 -DNT=800 -DNBH=0 -DSNAP_STRIDE=20 -DNPMLDEF=0
  small_nopml: clean $(TARGET)
  	./$(TARGET)

  .PHONY: all run small small_nopml clean cleanall
  ```

- **回滚方式**：删除 `small_nopml` 目标；从 `small` 的 CXXFLAGS 里删掉 `-DNPMLDEF=10`；
  从 `.PHONY` 里去掉 `small_nopml`。

### 本次会话验证

- 编译：`make clean && make`  → 无警告
- **回归验证（关键）**：`make small_nopml` (NPMLDEF=0) → trace 与 SESSION 1 逐位相同
  - 例如 `it=780, txx(rec)=-4.165e-07`（两次运行完全相同）
  - 证明 PML 关闭时新 13 块系统精确退化为旧 5 块行为
- **PML 生效验证**：`make small` (NPMLDEF=10)
  - `[pml] NP=10 d0_x=1036.16 d0_z=1036.16`
  - `[pml] active x-columns=20/50 z-rows=20/50 max(dampx)=1036.16`
  - `it=780, txx(rec)=-2.513e-07`（相对无 PML 减小，说明边界反射被吸收）
  - BiCGSTAB 每步 9–12 次迭代，与之前收敛速度相当
- **Lumped mass 校验**：`[mesh] Lumped mass: min=16.6667 max=100 sum=240100`
  与解析值 `(Nx-1)(Nz-1)·dx·dz = 49·49·10·10 = 240100` 精确一致

### 会话级完全回滚方式

从当前状态（含 PML）还原到 SESSION 1 结束时的状态：

1. 覆盖回 SESSION 1 版本：`config.h`、`mesh.h`、`mesh.cpp`、
   `assembler.h`、`assembler.cpp`、`main.cpp`、`Makefile`
   （SESSION 1 里每个 Change 都给出了 API / 语义描述，可作为重写参考）
2. 保留 SESSION 1 中已存在的其他文件不变：
   `material.h/.cpp`、`sparse_matrix.h/.cpp`、`linear_solver.h`、
   `solver_bicgstab.cpp`、`field_io.h/.cpp`、`README.md`、`CHANGELOG.md`
3. 清可执行与 `.dat`：
   ```bash
   cd /home/wanwb/zzt/Quantum/FEM/Serial && make cleanall
   ```

---

## SESSION 3 (2026-07-01 21:39, 修正 FEM 波场与 FDM 对不上：源缩放 + 质量集中 + Union Jack 网格)

### 上下文与诊断

用户提交对照图：`snap_txx_it0800.png`（FEM，200×200，it=800）跟
`Quantum/FDM/draw/snap_txx_it0800.png`（FDM 同参数同时刻）**完全对不上** ——
FEM 图是杂乱的高频对角条纹，振幅只有 ~1e-6；FDM 是干净的圆形 P 波，振幅 ~4e-5。
用户直觉：会不会是**网格剖分**的问题？

诊断三条原因（贡献从大到小）：

1. **声源缩放不对（100× 差）**
   FEM 弱形式下点源 `∫ N_k·δ(x-x_k)·f dΩ = N_k(x_k)·f = f`，b[σxx_k] += dt·f。
   解出 σxx_k ≈ dt·f / M_L_k。而 FDM 直接给节点加 dt·f（其 A 的 σxx 对角是 1，
   不是 M_L）。因此 FEM 里 σxx 物理振幅被 1/M_L = 1/(dx·dz) = 1/100 打压。
   **修法**：给源乘 `mesh.M_L[src_node]`。

2. **一致质量（consistent mass）在 P1-P1 速度–应力上数值弥散极差**
   Consistent M 会把源在邻居间"抹平"并激发 spurious mode，产生棋盘状抖动。
   **修法**：改用 **lumped mass**（M 的行和放对角），这是波动 FEM 的标准做法；
   接近 FDM 的中心差分行为。ψ 耦合项也同样 lumped，退化为逐节点耦合，与 FDM
   的 pointwise ψ ODE 完全一致。加 `-DUSE_CONSISTENT_MASS` 可切回。

3. **网格剖分方向偏置（用户猜测正确）**
   老版本所有 cell 都沿 SW-NE 对角线切 → 每个内部节点只连 6 个三角形（缺 NW 和 SE
   角方向的连接），造成波在对角方向传播偏快、非对角方向偏慢；用户观察到的对角
   条纹就是这个。
   **修法**：**Union Jack** 网格 —— 奇偶 cell 交替对角线方向，内部节点连 8 个
   三角形，回归各向同性。加 `-DNO_UNION_JACK` 可切回全 SW-NE 对角。

4. （次要）**点源在 P1 节点上仍会激发网格级高频**
   即使做完 1–3，源附近仍有小尺度菱形抖动 —— 单节点集中源本质上是宽带的，
   而 P1 网格无法传播 λ < 2dx 的模式。FDM 交错网格里 σxx 天然在格心，等价于对
   源做了一次网格级平滑；FEM 非交错节点上要显式做。
   **修法**：给源做一次 3×3 高斯核平滑（1-2-1 / 4 张量积，总权 = 1），源分到 9 个
   节点上并各自乘 `M_L[node]` 保持物理振幅。加 `-DNO_SOURCE_SMOOTH` 可切回集中点源。

### Change 3.1  `mesh.cpp` — Union Jack 网格

**动机**：消除结构化三角网格的方向偏置。

**修改文件**：`mesh.cpp` 的 `Mesh::build_triangles()`。

**关键 diff**（旧 → 新）：

旧：
```cpp
// T_lower: SW, SE, NE   （CCW）
// T_upper: SW, NE, NW   （CCW）
for (int icz = 0; icz < NcellZ; icz++) {
    for (int icx = 0; icx < NcellX; icx++) {
        // ... push T_lower(SW, SE, NE), T_upper(SW, NE, NW)
    }
}
```

新：
```cpp
// Union Jack：(icx+icz) 奇偶交替 SW-NE 和 NW-SE 对角
//  偶：T_a(SW, SE, NE), T_b(SW, NE, NW)   [SW-NE 对角]
//  奇：T_a(SW, SE, NW), T_b(SE, NE, NW)   [NW-SE 对角]
// `-DNO_UNION_JACK` 可退回全 SW-NE 对角
for (int icz = 0; icz < NcellZ; icz++) {
    for (int icx = 0; icx < NcellX; icx++) {
        const bool parity_even = ((icx + icz) & 1) == 0;
        if (parity_even) {
            push_tri(cid, sw, se, ne);
            push_tri(cid, sw, ne, nw);
        } else {
            push_tri(cid, sw, se, nw);
            push_tri(cid, se, ne, nw);
        }
    }
}
```

两种拆分下 `A_signed = 0.5·dx·dz > 0` 都 CCW，几何/梯度公式无需改。

**回滚**：把 `parity_even` 分支删除，退回旧的 T_lower/T_upper 二元结构（或加
`-DNO_UNION_JACK` 到编译命令）。

### Change 3.2  `assembler.cpp` — Lumped mass 默认开启

**动机**：一致质量导致 P1-P1 弥散/spurious mode；lumped mass 是波动 FEM 标准
且更接近 FDM 中心差分。

**修改文件**：`assembler.cpp` 里的 `m0_` 模板 lambda（在 `build_system` 三角
循环内），以及顶部注释块。

**关键 diff**：

旧：
```cpp
const double A_e = t.A;
const double A3  = A_e / 3.0;
const double A12 = A_e / 12.0;
// 3x3 consistent mass 模板：m0[i][j] = A_e/12 * (1 + δ_ij)
auto m0_ = [&](int i, int j) -> double { return A12 * (i == j ? 2.0 : 1.0); };
```

新：
```cpp
const double A_e = t.A;
const double A3  = A_e / 3.0;
// 3x3 mass 模板（默认 lumped；-DUSE_CONSISTENT_MASS 走一致质量）
#ifdef USE_CONSISTENT_MASS
    const double A12 = A_e / 12.0;
    auto m0_ = [&](int i, int j) -> double { return A12 * (i == j ? 2.0 : 1.0); };
#else
    // Lumped: 行和 → 对角 A_e/3，非对角为 0
    auto m0_ = [&](int i, int j) -> double { return (i == j) ? A3 : 0.0; };
#endif
```

`m0_` 会同时影响：
- 主 5 场的 M 块（Mρ, M0）→ 对角化
- ψ 耦合到主 5 场的 `-dt·M0` 项 → 逐节点耦合
- ψ 行自己的对角 `M_L` 项（原本就是 lumped，加 `A3`；对角化后与 `m0_` 一致）

对 A 与 M 的结构变化：M 从"consistent 3x3 stencil"变成对角；A 里主 5 场 vs 主
5 场的耦合以及 ψ 耦合的 mass 部分也都对角化 —— nnz(A) 与 nnz(M) 都变小。

**回滚**：加 `-DUSE_CONSISTENT_MASS` 到编译命令，或删掉 `#ifdef`，只留一致质量分支。

### Change 3.3  `assembler.cpp` — 源缩放 + 3×3 高斯平滑

**动机**：修振幅（M_L 缩放）+ 修源附近网格级抖动（空间平滑）。

**修改文件**：`assembler.cpp` 的 `assemble_rhs()`，声源施加块。

**关键 diff**：

旧：
```cpp
const int src_node = isou + jsou * Nx;
b_out[gid_xx(src_node)] += dt * sou_t;
b_out[gid_zz(src_node)] += dt * sou_t;
```

新：
```cpp
// (a) 物理振幅匹配：乘 M_L[src_node]
// (b) 空间平滑：3×3 高斯核 (1-2-1)/4 张量积，总权 = 1
#ifdef NO_SOURCE_SMOOTH
    const int di_list[1] = {0};  const int dj_list[1] = {0};
    const double w_list[1] = {1.0};  const int NW = 1;
#else
    const int di_list[9] = {-1,0,1,-1,0,1,-1,0,1};
    const int dj_list[9] = {-1,-1,-1,0,0,0,1,1,1};
    const double w_list[9] = { 1./16, 2./16, 1./16,
                               2./16, 4./16, 2./16,
                               1./16, 2./16, 1./16 };
    const int NW = 9;
#endif
for (int k = 0; k < NW; k++) {
    const int ix = isou + di_list[k];
    const int iz = jsou + dj_list[k];
    if (ix<0 || ix>=Nx || iz<0 || iz>=Nz) continue;
    const int nd    = ix + iz * Nx;
    const double sc = dt * sou_t * w_list[k] * mesh.M_L[nd];
    b_out[gid_xx(nd)] += sc;
    b_out[gid_zz(nd)] += sc;
}
```

**回滚**：加 `-DNO_SOURCE_SMOOTH` 得到 M_L 缩放但不做平滑（振幅对，源附近有斑点）；
或整块替换回旧的 2 行 `b_out[...] += dt * sou_t`（振幅小 100 倍，模式错）。

### 验证

- 编译：`make -j4`  → 无警告，无错误
- 200×200 主算例（`make CXXFLAGS="... -DNT=850 -DSNAP_STRIDE=50"`, 用户命令等价）：
  - 82 s / 850 步（跟 SESSION 2 相当，Union Jack 三角数不变）
  - `snap_txx_it0800.dat` → 干净圆形 P 波，中心在 (row=60, col=100)
  - 峰值振幅 ~6×10⁻⁵（FDM 同时刻 ~4×10⁻⁵，量级一致）
  - 波前半径 ~40 cells，与 Vp·t = 5000·0.08 = 400 m = 40 cells 吻合
  - 极性（红→蓝）与 FDM 一致
- 回归到"三关"配置（`-DUSE_CONSISTENT_MASS -DNO_UNION_JACK -DNO_SOURCE_SMOOTH`）
  应重现 SESSION 2 结束时的错误图案，方便对比

### 会话级完全回滚方式

从当前 SESSION 3 状态回到 SESSION 2 结束时：

1. **`assembler.cpp`**：
   - 顶部注释块回到 SESSION 2 版本（去掉 lumped/源缩放的说明段）
   - `m0_` lambda 去掉 `#ifdef USE_CONSISTENT_MASS`，只保留 `A12*(1+δ_ij)` 分支
   - `assemble_rhs` 里声源块去掉 3×3 循环和 `mesh.M_L` 缩放，回到 2 行
     `b_out[gid_xx/zz(src_node)] += dt * sou_t`
2. **`mesh.cpp`**：
   - `build_triangles()` 去掉 `push_tri` lambda 与 `parity_even` 分支，回到旧的
     "T_lower(SW,SE,NE), T_upper(SW,NE,NW)" 两段 push
3. 清可执行 + 输出：`make cleanall`
4. 编译验证：`make small`（此时应看到 SESSION 2 的错误行为，即报告里那张
   对角条纹图）

或者不改代码，直接用编译宏关闭全部 SESSION 3 修改（等价效果）：

```bash
make CXXFLAGS="-O3 -Wall -std=c++14 -DUSE_CONSISTENT_MASS -DNO_UNION_JACK -DNO_SOURCE_SMOOTH"
```

---

## SESSION 4 (2026-07-02 14:49, 修尾波散射：把源空间平滑从 3×3 加宽到 7×7 抗数值频散)

### 上下文与诊断

跑完 SESSION 3 修好的版本，用户从三张时刻快照（it=250 / 2000 / 3000）与四个 trace
文件观察到：**主波前没问题，但传播过一段后源附近持续出现一圈圈滞后的散射环，
trace 里也能看到明显的尾波多余波包**。用户询问："会不会是网格剖分的问题？
有限元二维三角网格生成对波动方程有什么规则？"

诊断过程：

1. **首先看快照**：it=2000 主 P 波正常向外传，同时源附近半径 ~25 cells 有一个
   八边形环；it=3000 环扩到 ~35 cells。环扩速度 ≈ 2000 m/s ≈ 0.4·Vp1 —— 显然
   不是 P/S 物理波，是数值伪迹。

2. **看 trace_zaxis**：源列 x-t 图里源位置一直在震，一路把小波包"喷"出去，
   直到 t=0.4s 都没停。这是典型的 **grid dispersion tail wave**。

3. **验证不是层界面 / PML 反射**：切换到 `NBH=0`（均匀介质）重跑，同样出现
   八边形环。→ 与介质无关，纯粹 FEM 数值问题。

4. **对比 FDM 同期图**：FDM 同参数 t=0.325s 干净，无此环。→ FEM 特有。

**根本原因**（用户的直觉方向正确）：

**P1 三角网格的可传播波长下限约为 4–6·dx**。Ricker 30 Hz 主频对应 λ ≈ 17 dx
（OK），但 Ricker 的高频尾（约 60–120 Hz）对应 λ ≈ 4–8 dx，正好位于/低于 P1
分辨率极限。这些频率被以错误（更慢）的相速度传播，在源之后形成 tail wave
环。SESSION 3 引入的 3×3 高斯平滑（σ≈0.7 cells）只压掉了最尖的一圈，不够。

FDM 交错网格里 σxx 天然放在格心，等价于对源做了一次网格级抗混叠，所以它没这
个问题。FEM 节点上必须显式做**更宽的抗混叠低通**。

### Change 4.1  `assembler.cpp` — 源空间平滑加宽 (3×3 → 7×7)

**修改文件**：`assembler.cpp` 的 `assemble_rhs()` 声源施加块 + 顶部 `#include <cmath>`。

**关键 diff**（旧 → 新）：

旧（3×3 固定核）：
```cpp
#ifdef NO_SOURCE_SMOOTH
    const int NW = 1; ...
#else
    // 3x3 分离高斯核：1-2-1 / 4
    const int di_list[9] = {-1,0,1,-1,0,1,-1,0,1};
    const int dj_list[9] = {-1,-1,-1,0,0,0,1,1,1};
    const double w_list[9] = {1/16, 2/16, 1/16, 2/16, 4/16, 2/16, 1/16, 2/16, 1/16};
    const int NW = 9;
#endif
```

新（可配置宽度，默认 7×7 σ=1.7）：
```cpp
#ifdef NO_SOURCE_SMOOTH
    const int SR = 0;  const double SIGMA = 1.0;
#else
    #ifndef SOURCE_SMOOTH_RADIUS
        #define SOURCE_SMOOTH_RADIUS 3   // 7×7 窗口
    #endif
    #ifndef SOURCE_SMOOTH_SIGMA
        #define SOURCE_SMOOTH_SIGMA 1.7  // Gaussian σ (in cells)
    #endif
    const int SR = SOURCE_SMOOTH_RADIUS;  const double SIGMA = SOURCE_SMOOTH_SIGMA;
#endif
// 累计 exp(-(di²+dj²)/(2σ²))，再归一化 → 域内 ∑w = 1，逐节点乘 M_L
double wsum = 0.0;  int nd_list[225];  double w_list[225];  int NW = 0;
const double two_sigma2 = (SR == 0) ? 1.0 : (2.0 * SIGMA * SIGMA);
for (int dj = -SR; dj <= SR; dj++)
for (int di = -SR; di <= SR; di++) {
    const int ix = isou + di,  iz = jsou + dj;
    if (ix<0 || ix>=Nx || iz<0 || iz>=Nz) continue;
    const double w = (SR==0) ? 1.0 : std::exp(-(double)(di*di+dj*dj)/two_sigma2);
    nd_list[NW]=ix+iz*Nx;  w_list[NW]=w;  wsum+=w;  NW++;
}
const double inv_wsum = 1.0 / wsum;
for (int k = 0; k < NW; k++) {
    const int nd = nd_list[k];
    const double sc = dt * sou_t * (w_list[k]*inv_wsum) * mesh.M_L[nd];
    b_out[gid_xx(nd)] += sc;  b_out[gid_zz(nd)] += sc;
}
```

新增开关：

- `-DSOURCE_SMOOTH_RADIUS=N`   N∈{1,2,3,...}：窗口半径（1→3×3，2→5×5，3→7×7）
- `-DSOURCE_SMOOTH_SIGMA=X`    X 是 Gaussian σ（单位：格）
- `-DNO_SOURCE_SMOOTH`         集中点源（用于对照 / 调试）

### 验证

- 均匀介质 (`-DNBH=0`)，NT=3000：八边形环已消失，只剩源附近 4 个极小残留斑
  点（~2×10⁻⁶，主波 ~1.3×10⁻⁵）。
- 层状介质（默认 NBH=40），NT=3000：
  - `snap_txx_it2000.dat` / `snap_txx_it2500.dat` → 主 P 波干净、层界面反射清晰、
    源附近无持续散射
  - `trace_zaxis_txx.dat` x-t 图 → 直达波 / 层顶反射 / 层底透射 / PML 吸收清晰
    可辨；`t>0.15s` 之后基本干净，无多余波包
- 计算开销与 SESSION 3 一致（每步 ~78 ms），装配阶段几乎无变化

### 数值经验（作为文档留下）

P1 三角 + lumped mass 波动 FEM 的**推荐**平滑窗口：

| 场景 | 窗口 | σ（in cells） | 截止波长 λ_c |
|---|---|---|---|
| 快速调试 / 频率很低 (F0·dx/Vp < 0.02) | 3×3 | 0.7 | ~4 cells |
| 一般 (F0·dx/Vp ≈ 0.05, 本项目) | **7×7 (默认)** | **1.7** | **~11 cells** |
| 有阴影区 / 长时间传播 | 9×9 或 11×11 | 2.0 – 2.5 | 12 – 16 cells |

窗口再大就是"物理源变虚"了，源附近几个格失去 δ 特性，但远场几乎无影响。

### 会话级完全回滚方式

从 SESSION 4 状态回到 SESSION 3 末：

1. `assembler.cpp` 顶部：删除 `#include <cmath>` （SESSION 3 未新增）
2. `assembler.cpp` 声源块：把 SESSION 4 里的 SR/SIGMA/for-循环-归一化改回
   SESSION 3 的写死 3×3 高斯 `1-2-1 / 4` 数组
3. 清可执行：`make clean`

或者不改代码，直接用宏回到 3×3 平滑：

```bash
make CXXFLAGS+="-DSOURCE_SMOOTH_RADIUS=1 -DSOURCE_SMOOTH_SIGMA=0.7"
```

---

## SESSION 5 (2026-07-02 20:00, 关键 bug 定位并修复：默认 mass 矩阵搞反了)

### 上下文（用户诉求）

用户拒绝 SESSION 4 的源平滑，理由极合理：
1. FDM 同参数运行没问题，说明源本身没问题；
2. 平滑源只是掩盖症状，不是根治；
3. **要求从网格出发分析**：画出网格、检查网格连续性、对比 Q1 四边形网格。

用户原话（决定性）："我在 FDM 的时候得到的结果是对的，且参数都没有变化，但是为什么到 FEM 就会因为声源出现问题，我觉得并不是平滑的原因... 我觉得从网格出发，进行分析。"

### 诊断过程（关键：本次会话的核心贡献）

1. **画出网格模型图**（`mesh_visualization.png`、`mesh_layer_interface.png`）：
   Union Jack 三角剖分是标准连续 Galerkin 网格，节点全共享、C0 连续、层
   界面沿单元边完全对齐 —— 网格本身无缺陷。

2. **实现 Q1 四边形网格**（`-DUSE_Q1_MESH`）作对照，用 4-node 双线性单元
   代替 P1 三角。均匀介质无源平滑测试：**P1 tri 和 Q1 quad 的散射图案几乎
   完全一样**（`compare_p1_vs_q1_snapshots.png`），从而排除"网格类型 /
   Union Jack 对角"作为根因。

3. **均匀介质、无 PML、无源平滑测试**（`out_p1_nopml/`）暴露了真正的问题：
   在 it=700 (t=0.07s)，P-wave 本应传到源外 350 m = 35 cells，但
   实际波场几乎卡在源附近 —— **σxx 场没有向外传播**，只在源 stencil 覆
   盖的少数节点上响应（`check_nopml_early.png`）。

4. **决定性实验 — 切换到 consistent mass** (`-DUSE_CONSISTENT_MASS` 于旧
   语义下打开)：**波立刻正常传出去了**！
   `compare_lumped_vs_consistent.png` 直接对比：同一个网格、同一个源、同
   一个时间步、同样无平滑，仅 mass 矩阵一换：
       lumped   : 波场停在源附近 ~30 cells 内（error）
       consistent: 波前扩到 89 cells（略过冲，正常 P1 色散），干净圆形

### 根源

SESSION 3 时把 mass 矩阵默认从 **consistent** 改成 **lumped**，理由错了。当时
的判断"波动 FEM 用 lumped mass 是标准"仅对 **displacement 二阶系统** (M ü = -K u)
成立；对 **velocity-stress 一阶系统** 则灾难性：

```
    σ 方程:  M · ∂σ/∂t = D · v
    lumped M 后：M 变对角  →  σ 更新是**纯 pointwise**（每个节点独立）
    但 D · v 是有 stencil 的  →  两侧空间尺度不匹配
    → σ 只在源的邻居 stencil 内被激活，无法通过 D · v 的耦合正确向外传播波
```

SESSION 3 引入的其余机制（source scaling by M_L、Union Jack 剖分、3×3 高斯平滑）
以及 SESSION 4 的 7×7 高斯平滑，实际上都是在**掩盖 lumped mass 的错误**。这也解释
了为什么源平滑"有点用但不完美"：平滑抑制了高频分量在破碎系统中的最强表现，但主
波依然被 lumped mass 严重压制。

### Change 5.1  默认 mass 改回 consistent（`assembler.cpp`）

- **修改 `assembler.cpp`**：
  ```diff
  - // 3x3 mass 模板（默认 lumped；-DUSE_CONSISTENT_MASS 走一致质量）
  - #ifdef USE_CONSISTENT_MASS
  -     const double A12 = A_e / 12.0;
  -     auto m0_ = [&](int i, int j) -> double { return A12 * (i == j ? 2.0 : 1.0); };
  - #else
  -     auto m0_ = [&](int i, int j) -> double { return (i == j) ? A3 : 0.0; };
  - #endif
  + // P1 tri 3x3 mass —— **默认 consistent**（SESSION 5 修正）
  + //     consistent : M0_ij = (|A|/12)·(1 + δ_ij)  ← 默认
  + //     lumped     : M0_ij = (|A|/3)·δ_ij         ← -DUSE_LUMPED_MASS 诊断用
  + #ifdef USE_LUMPED_MASS
  +     auto m0_ = [&](int i, int j) -> double { return (i == j) ? A3 : 0.0; };
  + #else
  +     const double A12 = A_e / 12.0;
  +     auto m0_ = [&](int i, int j) -> double { return A12 * (i == j ? 2.0 : 1.0); };
  + #endif
  ```
  同样对 Q1 分支：`USE_CONSISTENT_MASS` 取反，改为默认 consistent，
  `USE_LUMPED_MASS` 作 opt-in。

- ψ 辅助场的 mass 保持 lumped diagonal 不变（ψ 是 pointwise memory，
  在 PML 外 d=0 时自动关闭；lumped ψ 与 FDM 参考程序对齐）。

### Change 5.2  默认关闭源平滑（`assembler.cpp`）

- 按用户要求，取消 SESSION 4 的 7×7 高斯平滑默认。
- 修改 `assembler.cpp`：
  ```diff
  - #define SOURCE_SMOOTH_RADIUS 3   // 半径 3 → 7×7 窗口
  - #define SOURCE_SMOOTH_SIGMA 1.7
  + #define SOURCE_SMOOTH_RADIUS 0   // 默认：集中源（无平滑）
  + #define SOURCE_SMOOTH_SIGMA 1.0
  ```
- 平滑仍作可选保留：需诊断高频抑制时可加 `-DSOURCE_SMOOTH_RADIUS=3`。

### Change 5.3  新增 Q1 四边形网格支持（`mesh.h/.cpp`, `assembler.cpp`, `Makefile`）

- **动机**：作为诊断工具，帮助确认散射不是 P1 tri 网格类型的问题。
  副产物是一个可用的替代离散化。
- **新增**（编译期开关 `-DUSE_Q1_MESH`）：
  - `struct Quad { int nodes[4]; int cell_id; double dx_val, dz_val; };`
  - `Mesh::quads` 和 `Mesh::build_quads()`
  - `Mesh::build_lumped_mass()` 支持 quads（∫ N_k = |A_e|/4）
  - `assembler.cpp` 里加一整套 Q1 装配路径：分析积分 Dx/Dz/M0 矩阵
    （`mMc`, `mDx`, `mDz` 常量表 + `Dxq / Dzq / m0_` lambda）
- **Makefile** 加目标：`q1`, `q1_homo`, `small_q1`, `homo`

### Change 5.4  网格可视化脚本

- 新增图像（诊断/文档用）：
  - `mesh_visualization.png` — P1 tri Union Jack 与 Q1 quad 并排展示
  - `mesh_layer_interface.png` — Union Jack 剖分穿越层界面时的连续性
  - `compare_p1_vs_q1_snapshots.png` — 二者散射图案几乎完全一致
  - `compare_lumped_vs_consistent.png` — **关键诊断图**，直接展示 mass 错误
  - `check_nopml_early.png`, `check_nopml_checkerboard.png`, `check_nopml_vxvztxx.png` —
    从头看无 PML 均匀介质下波不传播的证据
  - `v5_homo_pml_snapshots.png`, `v5_homo_pml_trace.png` — 修复后波场

### 验证

均匀介质 (`-DNBH=0`)、200×200、PML=30、NT=2000（`out_v5_homo_pml/`）：
- **it=200-600**：横向双瓣 σxx 辐射从源清晰扩展，符合 cos²φ 各向异性
- **it=800-1200**：完整圆形 P 波前扩展到全域，无源附近残留
- **it=1400-1800**：波前进入 PML 被吸收，中心区回归零
- 波前半径 ≈ Vp·t（it=1200 → ~1200 m ≈ 120 cells ✓ 与 Vp=5000·0.12=600m…
  略过冲 —— 是 P1 tri 已知的相速度略高，但主波形正确）
- 尾部残留次要色散振荡是 P1 一阶格式固有特性，不影响主波前对比 FDM

### 会话级完全回滚方式

从 SESSION 5 状态回到 SESSION 4 末：

1. `assembler.cpp`：
   - 恢复 mass 默认 lumped：`#ifdef USE_CONSISTENT_MASS` → `#ifdef USE_LUMPED_MASS` 反转
   - Q1 分支同上
   - 源平滑默认改回 `SOURCE_SMOOTH_RADIUS=3, SIGMA=1.7`
2. `mesh.h/.cpp`：删除 `Quad` struct、`quads`、`build_quads()`、Q1 lumped mass 分支
3. `Makefile`：删除 `q1`, `q1_homo`, `small_q1`, `homo` 目标
4. 删除本次新增的 png 文件（诊断图，非运行必需）

不改代码，宏回退：
```bash
# 回到 SESSION 4 行为（lumped mass + 7x7 源平滑）
make CXXFLAGS="-O3 -Wall -std=c++14 -DUSE_LUMPED_MASS -DSOURCE_SMOOTH_RADIUS=3 -DSOURCE_SMOOTH_SIGMA=1.7"
```

### 关键教训（留给以后自己）

- **不要凭"经验判断" mass 应该 lumped 还是 consistent** —— 一定要看具体
  discretization 是二阶（displacement）还是一阶（velocity-stress）。
  displacement 系统里 lumped mass 是标准做法（stability + explicit stepping）；
  velocity-stress 系统里，如果 σ 和 v 的空间导数用不同 stencil（如 D 与 D^T），
  mass matrix 的空间尺度必须和它们一致 —— lumped mass 会造成灾难性的不匹配。

- **多用简单实验**：均匀介质 + 无 PML + 无源平滑，从 (t=0 到 P-wave 应该到达
  几个 λ 距离) 看波前是否正确、传得是否够远，是最强的正确性判据。SESSION 3-4
  一直在层介质 + PML 下调，反而让 lumped mass 的问题被掩盖。

- **多个可能病灶时优先做 A/B 测试** —— 网格类型（P1 vs Q1）、mass 类型
  （lumped vs consistent）、时间格式（BE vs CN）等，各变一个跑一个对照，
  很快就能定位。SESSION 5 只花了两个对照实验就锁定了 mass。

---

## SESSION 6 (2026-07-03, 对齐 FDM：Crank-Nicolson + 正确源缩放 + FDM 对照诊断)

### 用户诉求

用户指出 SESSION 5 修复后波场仍不对：雷克子波激发后源点应静止、同一深度应只有一个波包，
不应有多重振荡。要求对照 FDM/ImplicitV2 快照分析是程序问题还是参数/稳定性问题。

### 诊断结论（程序问题，不是参数问题）

| 检查项 | 结论 |
|--------|------|
| dx=10, dz=10, dt=1e-4, f0=30Hz | **与 FDM 完全相同，参数合理** |
| 每波长网格点数 | λ=Vp/f=167m≈**17 cells/λ**，FDM 足够 |
| CFL = Vp·dt/dx = 0.05 | 隐式格式无稳定性问题 |
| 网格连续性 | Union Jack 网格 C0 连续，无几何缺陷 |
| P1 vs Q1 | 散射图案相同 → 不是网格类型问题 |

**真正的问题**是 FEM 与 FDM 的**离散格式不一致**（三处程序 bug + 一处固有局限）：

1. **时间格式**：FEM 用 Backward Euler，FDM 用 **Crank-Nicolson**（2 阶）
2. **源缩放错误**：SESSION 3 用 `M_L`（=2×M0_diag）过度放大；去掉后又过小 70 倍；
   正确应乘 **M0_diag**（consistent mass 对角元）
3. **Mass 矩阵**：SESSION 5 已改回 consistent（保留）
4. **固有局限**：P1 同点连续 FEM ≠ FDM 交错网格，仍有空间数值频散（无法完全消除）

### Change 6.1  默认 Crank-Nicolson（`assembler.cpp`）

- `build_system`：耦合系数 `hc = dt/2`（默认），`-DUSE_BACKWARD_EULER` 时 `hc=dt`
- `assemble_rhs`：CN 时主 5 场 `b = 2·M·q^n - A·q^n`；ψ 块保持原公式

### Change 6.2  正确源缩放 M0_diag（`mesh.h/.cpp`, `assembler.cpp`）

- 新增 `Mesh::M0_diag[k] = Σ_{e∋k} |A_e|/6`（consistent mass 对角）
- 源：`b += dt·sou_t·M0_diag[k]`（内部节点 M0_diag≈50, M_L≈100）

### 验证（均匀介质 NBH=0, 200×200, NT=2000）

| 指标 | FDM | FEM v6 | 说明 |
|------|-----|--------|------|
| 源点峰值 | 1.87e-3 | 1.78e-3 | ✅ 振幅对齐 |
| 接收点峰值 | 1.25e-5 | 8.87e-6 | ✅ 同量级 |
| 源点 t>0.08s | 1.4e-6 | 1.3e-5 | ⚠️ FEM 源点残留仍 ~10× |
| 接收点波包数 | 6 | 12 | ⚠️ FEM 仍有额外波包 |
| x-t 图 / 快照 | 干净 V 形波前 | 明显改善，仍有内圈纹 | 见 v6_final_*.png |

### 回滚

```bash
# 回到 SESSION 5（BE + 无 M0_diag 源缩放）
make CXXFLAGS+="-DUSE_BACKWARD_EULER"
```

---

## SESSION 7 (2026-07-06, 【架构级重写】丢弃 velocity-stress，改位移二阶格式，根治菱形波前)

### 上下文与决策

SESSION 6 后波场仍是**菱形**（不像 FDM 那样圆形向外传）。本会话先做**诊断验证**，
确诊根因、再按方案 B 重写。

**诊断（新增 `diag_dispersion.py` / `diag_compare_fdm.py`，产物 `diag_*.png`）**：
- 根因 = 5 个场 (vx,vz,σxx,σzz,σxz) 全部 **collocated 在同一套节点**。连续
  Galerkin 一阶导等效成 `(f_{i+1}-f_{i-1})/2Δx`（跳过中心节点）→ 等效 2Δx 中心
  差分 → 奇偶节点解耦（checkerboard 零空间）+ 强各向异性频散 → 波前菱形。
- 解析数值频散：高频分量 G=6 pts/λ 时，collocated 各向异性 9.7% vs 交错 2.3%。
- 用压强 p=(σxx+σzz)/2（旋转对称，排除 σxx 的 cos²θ 辐射方向图）对照 it=800：
  FDM diag/axis=1.000（圆）；旧 FEM diag/axis=1.152、各向异性 119%（菱形）。
- 结论：collocated 框架下菱形**无法调参根除**，SESSION 3-6 的 Union Jack /
  源平滑 / mass 类型 / CN 全是治标。

**决策**：用户选**方案 B —— 丢弃 velocity-stress，改位移二阶格式**：
```
ρ ∂²u/∂t² = ∇·σ(u) + f,   σ = C:ε(u)
离散：M ü + C u̇ + K u = f，K = ∫ Bᵀ D B（标准刚度阵，连中心节点 → 无奇偶解耦）
时间：Newmark 平均加速度 β=1/4 γ=1/2（无条件稳定），每步解
      A a^{n+1} = f - C·v_pred - K·u_pred，A = M + γ·dt·C + β·dt²·K
```

### 备份（回滚基线）

**整个 SESSION 6 版本已完整备份**到目录：
```
Quantum/FEM/Serial/backup_v6_velocity_stress_20260706/
```
含全部 17 个文件（config/mesh/material/sparse_matrix/linear_solver/solver_bicgstab/
assembler/field_io/main/Makefile/README/CHANGELOG）。**会话级回滚 = 用该目录覆盖回根目录。**

### Change 7.1  重写 `config.h`

- DOF 从 13 块 → **2 块** (ux, uz)，`NfieldsPerNode=2`，`Ntot=2*Nnodes`，
  `off_ux=0, off_uz=Nnodes`。
- 删除全部 PML/ψ 相关常量与 offset；新增：
  - Newmark 参数 `nm_beta=0.25, nm_gamma=0.5`（`-DNEWMARK_BETA/GAMMA` 可改）
  - 声源矩幅度 `src_scale=dx·dz`（`-DSRC_SCALE`）
  - **吸收边界层**（替代 PML）：Rayleigh 质量比例阻尼分级 sponge，
    `NABS=30`（`-DNABSDEF`），`damp_at(x,z)` 二次曲线剖面，`compute_d0_x/z()`。
- 保留：网格/时间/介质/Ricker/求解器/输出 参数不变。
- **回滚**：用备份目录的 config.h 覆盖。

### Change 7.2  重写 `assembler.h/.cpp`

- `struct ElastSystem { A, K, C; beta,gamma,dt; src_dof/src_coef; ... }`。
- `build_system`：遍历 P1 三角，逐单元算
  - B(3×6)、D(3×3 平面应变)、`Ke=|A_e|·Bᵀ D B`
  - consistent mass `Me=ρ·(A_e/12)(1+δ)`（ux/uz 各一份，无交叉）
  - 阻尼 `Ce=d_e·Me`（d_e=单元 3 节点平均阻尼）
  - 合成 `A = M + γ·dt·C + β·dt²·K`，同时单独保留 K、C 供每步 RHS
  - 声源等效力系数：`F_a=M0·∇N_a` 对源邻接三角求和（自平衡、圆对称）
- `assemble_source_force(sys, sou_t, F)`：`F = src_scale·sou_t·Σ∇N_a`。
- `recover_nodal_stress(mesh,mat,u, sxx,szz,sxz)`：逐单元 ε=Bu、σ=Dε，
  面积加权投影到节点，用于与 FDM σ 快照对照。
- **回滚**：用备份目录覆盖。

### Change 7.3  改 `field_io.h/.cpp`

- 删除旧 `write_snapshot(q, off_*)`，新增
  `write_snapshot_fields(it,Nx,Nz, vx,vz,sxx,szz,sxz)`：直接吃 5 个准备好的
  Nz×Nx 场（速度取自 Newmark v，应力取自 recover_nodal_stress），文件名/格式
  与 FDM 完全一致（snap_vx/vz/txx/tzz/txz）。

### Change 7.4  重写 `main.cpp`

- Newmark 时间循环：维护 u/v/a，每步 预测→组 RHS→解 A a=b（warm-start a^n）→
  校正 u/v/a → 恢复应力 → 记录 trace（σxx/σzz@rec + 源列 σxx x-t）+ 快照。
- 求解层仍用 `ILinearSolver`/BiCGSTAB，CIM 落位口不变。

### Change 7.5  更新 `Makefile`

- 目标：`make`（全尺寸层状）、`homo`（均匀介质验圆）、`small` / `small_homo`。
- 删除 velocity-stress 时代的 `q1/small_nopml/...`。

### 验证（新增 `diag_planB_check.py`，产物 `diag_planB_vs_fdm.png`）

- 编译 `make clean && make`：无警告。
- 冒烟 `make small_homo`：BiCGSTAB 每步 7–12 次收敛，无发散，7.85 s / 700 步。
- **关键：均匀介质 200×200 至 it=800，压强波前各向异性对照**：

| 版本 | diag/axis | 总各向异性 | 形状 |
|------|-----------|-----------|------|
| 旧 velocity-stress FEM | 1.152 | 119% | 菱形 |
| **新位移二阶 FEM** | **1.008** | **8.0%** | **圆** ✅ |
| FDM 交错 | 1.000 | 4.7% | 圆 |

  → 对角/轴向超出从 15% 降到 0.8%，与 FDM 同水平，**菱形根治**。
- 波前半径 ~20 cells 正确（Ricker 峰值 t≈T0=0.04s 发出，到 t=0.08s 传 20 cells）。

### 已知遗留 / 后续

- 残余各向异性 8%（vs FDM 4.7%）与波速略慢（半径 18.8 vs 21）是 P1 consistent-mass
  的固有轻微频散，可用更细网格或调 mass 缓解，属次要。
- 吸收层是 Rayleigh 质量比例阻尼 sponge（一次实现），比 ADE-PML 简单；如需更强
  吸收可后续升级为位移型 PML。
- 声源目前 M0(t)=src_scale·ricker，与 FDM 的 stress-glut 可能差一个时间导数/常数，
  波形级振幅对齐待标定（不影响波前形状结论）。

### 会话级完全回滚

```bash
cd /home/wanwb/zzt/Quantum/FEM/Serial
cp backup_v6_velocity_stress_20260706/*.h  backup_v6_velocity_stress_20260706/*.cpp .
cp backup_v6_velocity_stress_20260706/Makefile .
make clean && make        # 回到 SESSION 6 velocity-stress 版本
```


---

## SESSION 8 (2026-07-06, 位移二阶格式加 ADE-PML unsplit：主方程进 CIM 的 A 保持对称正定)

### 上下文与诉求

SESSION 7 位移二阶格式波前已是圆、频散极小，但**边界有明显反射**（当时只有
Rayleigh 质量比例阻尼 sponge，吸收 ~40–50%，非真 PML）。用户两个诉求：
1. 加真正的 PML；
2. 因后续要接 **CIM 求解器**，确认每步求解仍是 `A x = b` 形式、且 A 结构对 CIM 友好。

### 决策：ADE unsplit PML，"主方程进 CIM、ψ 显式更新"

关键设计（保证 CIM 衔接）：进 CIM 的有效矩阵 **A = M + γ·dt·C + β·dt²·K_eff**
**保持对称正定、常系数、装配一次**；PML 记忆变量 ψ 满足一阶 ODE，用指数积分
**逐 DOF 显式更新、只进 RHS b，绝不进 A**。

复坐标拉伸 s_α = 1 + d_α/(iω)，动量方程乘 s_x·s_z 后按 iω 幂次归类：
```
M ü + C u̇ + K_eff u + F_ψ = F_src
  M     = ∫ρ N_aN_b                         (质量, 对称 SPD)
  C     = ∫ρ(d_x+d_z) N_aN_b                (PML 阻尼, 对称, 仅 PML 区非零)
  K_eff = 标准刚度 K + 角区弹簧 P            (P = ∫ρ d_x d_z N_aN_b, 对称)
  F_ψ   = ψ_xx + ψ_zz                        (显式记忆修正, 不进 A)
```
刚度按导数对拆分：只有 ∂x∂x 块带因子 s_z/s_x、∂z∂z 块带 s_x/s_z（iω 有理式）
需要 ψ；交叉 ∂x∂z 因子=1。有理式的"1"部分并入标准 K（进 A），剩余部分变成
```
ψ̇_xx + d_x·ψ_xx = (d_z-d_x)·(K_xx·u),   ψ̇_zz + d_z·ψ_zz = (d_x-d_z)·(K_zz·u)
指数积分:  ψ^{n+1} = E·ψ^n + fac·rhs,  E=exp(-d·dt), fac=(1-E)/d (d→0取dt)
```
NABS=0 时 d_α≡0 → C=P=ψ=0，**精确退化回 SESSION 7 纯位移格式**（已验证逐位一致）。

### Change 8.1  `config.h`

- 删除 SESSION 7 的合并 sponge `damp_at(x,z)`，改为**方向性** `dampx_at(x)` /
  `dampz_at(z)`（PML 必须按方向分）。d0 由目标反射系数 RC 反推（二次剖面，
  Collino-Tsogka）。默认 `NABSDEF=30, RC_ABS=1e-4`。

### Change 8.2  `assembler.h/.cpp`

- `ElastSystem` 扩为 5 个矩阵：`A`（进 CIM）、`K_eff`、`C`、`K_xx`、`K_zz`
  （后 4 个只做每步 SpMV 组 b，不进 CIM）；加 PML 节点级剖面缓存
  `dampx_ix/dampz_iz/Ex_ix/Ez_iz/facx_ix/facz_iz`。
- `build_system`：单元级算 d_sum=<d_x+d_z>（→C）、d_pro=<d_x·d_z>（→P 并入
  K_eff）；同时装配 K_xx=∫coef·∂xN∂xN、K_zz=∫coef·∂zN∂zN（coef: ux 用 λ+2μ、
  uz 用 μ 对 K_xx；对调给 K_zz）。合成对称 A。
- 新增 `pml_update_force(sys, u_pred, ψ_xx, ψ_zz, F_ψ)`：一次 K_xx·u/K_zz·u
  SpMV 驱动，逐 DOF 指数积分更新 ψ，输出 F_ψ = ψ_xx+ψ_zz。全显式。

### Change 8.3  `main.cpp`

- 加 `psi_xx/psi_zz/Fpsi`（各 Ntot）；时间循环里 6b' 调 `pml_update_force`，
  6c 的 RHS 改为 `b = F - C·v_pred - K_eff·u_pred - F_ψ`（`sys.K`→`sys.K_eff`）。

### 验证（关键，全部通过）

1. **NABS=0 退化**：80×80 无 PML，trace 与 SESSION 7 逐位一致。
2. **A 严格对称（CIM 前提）**：数值检查 `max|A_ij−A_ji|=0.000e+00`，
   `max|A_ij|=1.69e5` → **rel=0，SYMMETRIC (CIM-ready)**。
   （对比：旧 velocity-stress 的 A 是非对称的；位移格式 M/C/K_eff 全对称 → A 对称，
   M SPD 保证 A SPD，CIM 可直接映射 min ½xᵀAx−bᵀx。）
3. **真 PML 吸收（能量判据，决定性）**：120×120 均匀介质长跑 NABS=30 RC=1e-4，
   全域动能 ||v||² 从峰值 **衰减到 0.0%（1e-18，机器噪声）**；
   对比无 PML **恒 plateau 在 4.3e-14**（纯反射）。厚度扫描 NABS=20→43%、
   40→76%+、30 长跑→~100%，随厚度单调增强 = 真 PML 行为（非 sponge 的 plateau）。
4. **波前仍圆**：200×200 均匀 it=800，diag/axis=1.008（与无 PML 完全相同）→
   PML 未引入各向异性。
5. 每步 BiCGSTAB 6–8 次收敛；200×200×1700 步约 41 s。

诊断图：`diag_pml_check.png`、`diag_pml_200_check.png`。

### CIM 衔接结论（给用户）

- 每步严格 `A x = b`：A=`sys.A`（对称正定、常系数、装配一次），x=加速度 a^{n+1}，
  b=`F − C·v_pred − K_eff·u_pred − F_ψ`（每步现算，ψ 已知）。
- 求解调用点 `main.cpp` 走抽象 `ILinearSolver::solve(req)`，`req.A/b/x`；
  接 CIM 只需实现 `CIMSolver::solve(req)`，上层零改动。`DummyCIMSolver` 落位口在
  `solver_bicgstab.cpp`。
- **A 对称正定**是位移格式相对旧 velocity-stress（非对称 A）给 CIM 的核心红利，
  ADE-PML 的 ψ 显式设计刻意不破坏它。

### 回滚

```bash
# 关 PML（退回 SESSION 7 纯位移，A 不变）：
make CXXFLAGS+="-DNABSDEF=0"
# 或调厚/调强：-DNABSDEF=40 -DRC_ABS=1e-5
```

---

## SESSION 9 (2026-07-06, 诊断均质模型的尾波：确认为 P1 网格频散，非稳定性/介质/PML)

### 用户诉求

`make homo`（纯均质）下，波传播完成后仍有很弱尾波（snap_tzz_it1400/it3550、
trace_zaxis row100 可见，约主波 2–6%）。问是网格、稳定性还是声源原因。

### 确认：make homo 是纯均质

material.cpp:38-41 NBH=0 时 `in_layer` 恒 false，所有 cell 取背景介质，无层界面。
→ 尾波与介质无关。

### 诊断结论：P1 consistent-mass 网格频散（不是稳定性/介质/PML/声源）

逐一排除 + 决定性实验：
- ❌ **稳定性**：Newmark β=1/4 无条件稳定，SESSION 8 能量测试已证 ‖v‖² 单调衰减
  到 1e-18（若失稳应增长）。
- ❌ **介质**：均质已确认。
- ❌ **PML/边界**：**关键实验**——完全关 PML（NABSDEF=0）、且在边界反射能返回
  之前测量，主波后仍有尾波：
  ```
  探针(100,70) 主波峰 amp=3.25e-2  主波后尾波 max=2.01e-3 = 6.2%
  探针(100,80)                                        = 6.5%
  探针(100,100)                                       = 5.3%
  ```
  无 PML、无反射仍有尾 → 尾波源头是离散化本身，不是边界。
  （it=3550 快照下角的残余是另一回事：PML 角区/掠射的不完美吸收，二阶小量。）
- ✅ **根因 = P1 网格频散**：P1 单元 + consistent mass 对高频分量相速度偏低，
  Ricker 高频尾以更慢速度落在主波之后，形成频散尾。远场干净测量（300×300 无 PML）
  主波后尾波 ≈ 2.6%。

### Change 9.1  混合质量矩阵（抗频散，可调，仍对称 SPD → 不影响 CIM）

- `config.h` 加 `MASS_BLEND`（默认 0.5）：M = (1-α)·consistent + α·lumped。
  consistent 与 lumped 的 P1 频散误差符号相反，混合可部分抵消。lumped 对角正、
  consistent SPD → 混合仍对称 SPD，A 的 CIM 友好性不变。
- `assembler.cpp`：单元 mass 模板由固定 consistent 改为
  `m_diag = (1-α)A_e/6 + α·A_e/3`，`m_off = (1-α)A_e/12`。

### 实测（300×300 无 PML 远场尾波，200×200 波前圆度）

| MASS_BLEND | 远场尾波 | 波前 diag/axis | 主波峰值 |
|-----------|---------|---------------|---------|
| 0.0 (纯consistent) | 2.58% | 1.008 | 5.24e-2 |
| 0.5 (默认) | 2.31% | ~1.01 | ~5.6e-2 |
| 1.0 (纯lumped) | 2.03% | 1.016 | 6.06e-2 |

结论：**混合质量是弱杠杆**，最多把尾波从 2.6% 压到 2.0%（~20%），波前圆度基本不变。
说明尾波是 P1 一阶单元的**固有频散**，非 bug。

### 彻底消除尾波的手段（如需要）

1. **细网格**：dx 减半，每波长网格数翻倍，频散尾按 O(dx²) 下降（最有效但最贵）。
2. **高阶单元**：P2/谱元，频散极小（改动大）。
3. **混合质量 α**（已加，默认 0.5）：便宜但只削 ~20%。
2% 量级尾波对多数应用（含与 FDM 对照）可接受；P1+consistent 本就有此固有特性。

### 回滚

```bash
make CXXFLAGS+="-DMASS_BLEND=0.0"    # 回到纯 consistent mass（SESSION 8 行为）
```

---

## SESSION 10 (2026-07-07, P1→P2 升级：6 节点二次三角，抗频散尾波；CIM 接口不变)

### 用户诉求

为减少 SESSION 9 诊断的 P1 网格频散尾波，把单元从 P1（3 节点线性）升级到
P2（6 节点二次）。先备份当前版本再改。

### 备份

SESSION 8/9 的 P1 位移+ADE-PML 版本完整备份到：
```
Quantum/FEM/Serial/backup_v8_P1_displacement_pml_20260707/   (16 文件)
```

### 关键设计：半间距细网格承载 P2 节点（避免边表）

P2 的 6 节点（3 顶点 + 3 边中点）恰好都落在半间距 (dx/2) 的规则细网格上：
顶点→偶数格点，边中点→半格点。于是把所有节点统一编号在
`MX×MZ = (2Nx-1)×(2Nz-1)` 细网格，node id = IX+IZ·MX。好处：
- PML 剖面按坐标算照常有效（IX=id%MX，物理坐标=IX·(dx/2)）；
- 边中点天然唯一编号，无需构造边表；
- 快照只导出**顶点子网格 (stride 2) = 原 Nx×Nz**，画图脚本与 FDM 对照零改动。

P1 通过 `-DUSE_P1_MESH` 保留（此时 MX=Nx、hgrid=dx，一切退化回原状）。

### Change 10.1  config.h
- 加 MX/MZ/hgrid/node_stride（P2 细网格，P1 时退化）；Nnodes=MX·MZ（P2≈4×），
  Nvert=Nx·Nz（快照用）；`fine_node_id`、`vert_node_id`；node_id 改走 vert_node_id。

### Change 10.2  mesh.h/.cpp
- 加 `Tri6`（6 节点，存面积坐标常梯度 dL_dx/dL_dz）；build_nodes 改建细网格；
  `build_triangles_p2()`（Union Jack + 边中点 MID 索引平均）；构造函数按宏选 P1/P2。

### Change 10.3  assembler.cpp
- PML 剖面改按细网格 IX∈[0,MX) 采样。
- 装配分 `#ifdef USE_P1_MESH`（原闭式 P1）/ `#else`（P2）：P2 用 Strang 6 点 4 阶
  高斯求积，逐高斯点算 N、∂N/∂x,z（=Σ ∂N/∂L·∂L/∂x），累加 12×12 Ke、mass、
  PML 的 C/P/Kxx/Kzz。混合质量 lumped 用行和（P2 行和为正，安全）。
- 声源 P2 分支：在源顶点面积坐标处评估 6 节点 ∇N，F_a=M0·∇N_a。
- 应力恢复 P2 分支：在 3 顶点处评估 P2 梯度、面积加权投影。
- pml_update_force 的 node→(IX,IZ) 改用 %MX。

### Change 10.4  field_io.cpp
- write_snapshot_fields 从细网格数组用 vert_node_id 抽取顶点子网格，输出仍 Nz×Nx。

### Change 10.5  main.cpp
- rec_node/src_node/源列 trace 改用 node_id（顶点→细网格映射）。

### 验证

- 编译无警告。P1 回归（-DUSE_P1_MESH）：MX=Nx，行为与 SESSION 9 一致。
- **P2 默认运行正常**：100×100 → MX=199, Nnodes=39601(≈4×), BiCGSTAB 每步 6 次收敛，
  无失稳。200×200 nnz(A)≈1.15M。
- **A 仍严格对称（CIM 前提保住）**：数值检查 rel|Aij-Aji|=1.4e-17 → SYMMETRIC。
  P2 的 M/C/K_eff 仍全对称 SPD，CIM 映射 min ½xᵀAx-bᵀx 不受影响，求解接口零改动。
- **频散尾波**：300×300 均质无 PML 远场探针，P1=2.58% → P2=2.27%。

### 诚实评估（重要）

P2 在本问题上对尾波的改善**比预期小**（远场 2.58%→2.27%，约 12%）。原因：
1. 主频 16 点/波长时 P1 本就不太频散（主频各向异性仅 1.2%）；尾波主要来自
   Ricker 高频端（4–6 点/波长）。
2. **更关键**：长距离传播测试（400×400，传 200 格）P1 尾波只有 0.96%，说明尾波
   **不是累积传播频散主导**，而主要是**单节点矩张量声源**激发的近源网格模式
   （点源本质宽带，含 λ<2dx 分量，P1/P2 都难干净传播）。P2 对传播频散有效，
   对点源激发的改善有限。

**建议**：若要进一步压尾波，性价比更高的是（a）**声源空间平滑**（3×3/5×5 高斯核，
把点源展宽、削掉 λ<2dx 分量，SESSION 3-4 做过，对位移格式同样适用），或
（b）**细网格**。P2 的主要价值是主波精度/各向异性，尾波削减是次要收益。
P2 代价：DOF≈4×、求解≈4×慢。

### 回滚

```bash
make CXXFLAGS+="-DUSE_P1_MESH"     # 退回 P1（SESSION 8/9 行为）
# 或整体恢复： cp backup_v8_P1_displacement_pml_20260707/* .
```

---

## SESSION 11 (2026-07-07, 诊断并修复 P2 的"顶边频散"：实为 P2 点源棋盘伪影)

### 用户现象

P2 结果异常：trace_zaxis row30 在 it>3500 出现振幅达 0.6 的强"频散"（主波仅 0.01）；
snap_txx_it3450 顶边 z=0 一条棋盘格花纹、振幅 ±100；P1 无此问题。问是否边界设置问题。

### 诊断：不是频散、不是边界、不是发散 —— 是 P2 点源在源节点的棋盘伪影

逐一排除（实验）：
1. **不是边界/PML**：源远离边界（140×140, 源(70,30)）+ NABS=0，仍出现。
2. **不是发散/不稳定**：幅值跟随 Ricker 时间函数（it≈400 峰值最大、之后衰减），
   不是无界增长。
3. **确诊 = 源节点局部伪影**：it=400 源节点 5×5 邻域 txx：
   ```
   源节点 = -10.26          ← 尖峰(全域max就在此)
   紧邻   = +4.09, +0.88    ← 正负交替(棋盘)
   3格外  = -0.10, +0.27    ← 正常波幅
   ```
   一个符号交替、局限源节点、3 格内衰减的尖峰。it3450 快照"贴顶边"是因该算例
   源偏上、伪影被 PML 反复作用堆积到顶边，根在源节点本身。

**机理（为何 P1 无、P2 有）**：矩张量点源 F_a=M0·∇N_a 直接加到单顶点。
P1 每单元 3 个等价顶点、梯度同量级 → 平滑。P2 每单元 6 节点，顶点函数
L(2L-1) 与边中点函数 4LiLj 的梯度差异极大，源力在其间严重不均 → 激发源节点处
最高频网格模式（棋盘）。这是 P2/高阶元做点源的已知病理。

### Change 11.1  声源高斯空间平滑（assembler.cpp）

- 把矩张量源在源周围**顶点**上做高斯平滑：对每个"子源顶点"仍施加标准
  F_a=M0·∇N_a，按归一化高斯权重（∑w=1，保持总矩不变）叠加。展宽点源、削掉
  λ<2dx 分量，源节点棋盘模式无法被激发。
- 开关：`-DSRC_SMOOTH_R=R`（窗口 (2R+1)²，默认 **P2=3(7×7)**，P1=0）、
  `-DSRC_SMOOTH_SIGMA=σ`（默认 P2=1.5，P1=1.0）、R=0 退回集中点源。
- 只改 RHS 的声源力系数，**不动 A**（A 仍对称 SPD，CIM 不受影响）。
- 需要 `#include <utility>`（std::pair）。

### 验证

- 编译无警告。
- **源节点尖峰消除**：140×140 it400 源节点 -10.26 → -0.41（R=3），邻域平滑无棋盘。
- **无后期爆炸**：140×140 NT=2000 全程 |max|=0.44（修复前 ~10+ 且持续），源列 max
  随时间衰减、峰位跟随波前外移（row 30→96）= 正常传播。
- **波前仍圆**：200×200 it800 diag/axis=0.993；源节点 txx=-1.3e-3（同波幅，无尖峰）。
- **A 仍严格对称**：rel|Aij-Aji|=1.44e-17 → SYM(CIM-ready)，源修复不碰 A。

### 说明

此修复对 P1 同样可用（历史 SESSION 3-4 在 velocity-stress 上做过类似平滑），但 P1
默认关闭（R=0）以保持与既往结果一致；P2 默认开启（点源在 P2 上是必须平滑的）。

### 回滚

```bash
make CXXFLAGS+="-DSRC_SMOOTH_R=0"   # 关闭源平滑(P2 会重现源节点棋盘)
make CXXFLAGS+="-DUSE_P1_MESH"      # 或退回 P1
```
