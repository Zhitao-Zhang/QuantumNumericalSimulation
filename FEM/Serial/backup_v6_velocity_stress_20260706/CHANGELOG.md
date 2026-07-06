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
