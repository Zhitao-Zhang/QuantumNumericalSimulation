# 2D 速度–应力弹性波有限元（P1 三角，串行 CPU，第一版）

本目录实现文档 `CIM_FEM_Ax_b_guide_velocity_stress_v2.md` 中「阶段 A」的
要求，用于**与 `Quantum/FDM/ImplicitV2/main.cpp` 的输出直接对照**。

---

## 1. 数学模型速览

物理方程（二维速度–应力一阶弹性波，未加 PML）：

$$
\rho \frac{\partial v_x}{\partial t} = \frac{\partial \sigma_{xx}}{\partial x} + \frac{\partial \sigma_{xz}}{\partial z} + f_x
$$

$$
\rho \frac{\partial v_z}{\partial t} = \frac{\partial \sigma_{xz}}{\partial x} + \frac{\partial \sigma_{zz}}{\partial z} + f_z
$$

$$
\frac{\partial \sigma_{xx}}{\partial t} = (\lambda+2\mu) \frac{\partial v_x}{\partial x} + \lambda \frac{\partial v_z}{\partial z}
$$

$$
\frac{\partial \sigma_{zz}}{\partial t} = \lambda \frac{\partial v_x}{\partial x} + (\lambda+2\mu) \frac{\partial v_z}{\partial z}
$$

$$
\frac{\partial \sigma_{xz}}{\partial t} = \mu \left( \frac{\partial v_x}{\partial z} + \frac{\partial v_z}{\partial x} \right)
$$

FEM 空间离散：**5 个未知量全部用连续 P1 三角形单元**（`README §7.2`）。

对速度方程做分部积分（自由面边界 → 边界项自然为零），应力方程直接测试。
最终得到半离散系统

$$
M \dot q + K q = f,\qquad q = (v_x, v_z, \sigma_{xx}, \sigma_{zz}, \sigma_{xz})^T
$$

时间用 **Backward Euler**：

$$
(M + \Delta t\, K)\, q^{n+1} = M q^n + \Delta t\, f^{n+1}
$$

于是每步就是显式的 $Ax=b$，供任意后端（BiCGSTAB / CIM / QUBO）求解。

---

## 2. 目录结构

```
Serial/
├── config.h                # 全部参数（对齐 FDM/ImplicitV2）
├── mesh.h / mesh.cpp       # 结构化 P1 三角网格 + 单元几何
├── material.h / .cpp       # 逐 cell 的 (ρ, λ, μ)（含高速夹层）
├── sparse_matrix.h / .cpp  # CSR + RowMap 装配容器
├── linear_solver.h         # ILinearSolver 抽象接口（CIM 替换口）
├── solver_bicgstab.cpp     # Jacobi 预条件 BiCGSTAB + DummyCIM 落位
├── assembler.h / .cpp      # 从 mesh + material → (A, M)
├── field_io.h / .cpp       # wfile2d 与 FDM 同格式的 ASCII 输出
├── main.cpp                # 时间推进主循环
├── Makefile
└── README.md
```

**离散层 ↔ 求解层严格解耦**（文档 §3.2）：
- `assembler::build_system(...)` 只填 `A` 和 `M`
- `assembler::assemble_rhs(...)` 只填 `b`
- `ILinearSolver::solve(SolveRequest)` 只求 `x`

未来把 BiCGSTAB 换成 CIM，只需实现一个新的 `ILinearSolver`。

---

## 3. 与 FDM 的参数一一对齐

| 参数 | FDM（ImplicitV2/main.cpp） | FEM 本目录 |
|---|---|---|
| 网格 | `NZ=200, NX=200` | 一致 |
| 步长 | `dx=dz=10 m` | 一致 |
| 时间 | `dt=1.0e-4 s, NT=4500` | 一致 |
| 声源 | Ricker 30 Hz，`jsou=NZ/2-40, isou=NX/2` | 同节点 `(isou=100, jsou=60)` |
| 接收 | `jrec=NZ/2+10, irec=NX/2` | 同节点 `(irec=100, jrec=110)` |
| 固体背景 | `Vp=5000, Vs=3300, ρ=2450` | 一致 |
| 高速夹层 | `Vp=6000, Vs=3400, ρ=2770, iz∈[100,120)` | 一致 |
| PML | ADE-PML（8 个辅助场） | **不启用**（文档 §14） |
| 求解 | Jacobi-BiCGSTAB | 一致 |
| PML | ADE-PML（8 个 ψ 辅助场）NP=30 | ADE-PML（8 个 ψ 辅助场）NP=30 |

---

## 3.5 PML 实现（本项目 vs FDM）

**共同点**：都用 ADE-PML，都引入 8 个 ψ 辅助场并入 `Ax=b`，ψ 采用指数-梯形离散：

$$
\psi^{n+1} + (d\Delta t/2) g^{n+1} = E \cdot [ \psi^n - (d\Delta t/2) g^n ], \quad E = e^{-d\Delta t}
$$

**FEM 具体做法**（本目录）：

- 所有 ψ 场用 P1 连续，与主 5 场共节点，共 `13 * Nnodes` DOF；
- ψ 行的“质量”用 **lumped mass**（`M_L[k] = ∫ N_k dΩ`），保证 `d=0` 处 `ψ^{n+1}=ψ^n`；
- ψ 与主场的耦合矩阵是 `G_S,α_ij = ∫ N_i ∂N_j/∂α dΩ = dN_dα[j]·|A|/3`（P1 三角闭式）；
- `d_α(x)` 只按节点采样一维数组 `dampx[ix]`, `dampz[iz]`（比 FDM 简单，因为无交错网格）；
- `NPMLDEF=0` 时全部 ψ 行退化为 `M_L·ψ^{n+1}=M_L·ψ^n`，程序自动关掉 PML，与原
  「无 PML」版本 trace **逐位一致**。

---

## 4. 编译与运行

```bash
cd Quantum/FEM/Serial
make                # 默认 200x200 x 4500 步
./fem_elastic_2d
```

冒烟测试（≈几秒钟即可跑通全流程）：

```bash
make small          # 50x50, 800 步, NP=10 (PML 开)
make small_nopml    # 50x50, 800 步, NP=0  (PML 关，用于回归对比)
```

`small_nopml` 用来做**回归验证**：新版（13 块 + PML 关闭）应与老版
（5 块无 PML）trace 逐位相同——已实证。

也可以按需覆盖任意参数：

```bash
make CXXFLAGS+="-DNT=2000 -DBICG_TOL=1e-8"
```

### 4.1 数值伪迹修正的编译开关

`velocity-stress` FEM 的 4 个坑，按重要性排序（重要教训见 CHANGELOG SESSION 5）：

| 宏 | 默认 | 打开时行为 |
|---|---|---|
| `USE_BACKWARD_EULER`    | 未定义（→ **Crank-Nicolson**，SESSION 6，对齐 FDM） | 退回 Backward Euler |
| `USE_LUMPED_MASS`       | 未定义（→ **consistent mass**，SESSION 5 修正） | 强制 lumped mass（破坏波传播，仅诊断） |
| `USE_Q1_MESH`           | 未定义（→ **P1 三角 Union Jack**） | 换成 Q1 双线性四边形网格（每 cell 一个单元）；DOF 布局不变，可与 P1 直接对比 |
| `NO_UNION_JACK`         | 未定义（→ **Union Jack 三角剖分**） | 回退到"全 SW-NE 对角"，波场沿对角方向偏置（P1 tri 才生效） |
| `SOURCE_SMOOTH_RADIUS`  | **0**（→ 集中源） | ≥1 时开高斯空间平滑；SESSION 4 曾默认 3，SESSION 5 关掉 |
| `SOURCE_SMOOTH_SIGMA`   | 1.0 cells | Gaussian σ（仅在 RADIUS≥1 时使用） |
| `NO_SOURCE_SMOOTH`      | 未定义 | 显式关闭平滑（与默认等价，冗余保留兼容旧脚本） |

`USE_CONSISTENT_MASS` 现在是**噪声名**（默认就已经是 consistent），但保留可用不
影响任何行为。旧脚本里带这个宏也照常工作。

快速对照命令：

```bash
# 默认（推荐）
make

# Q1 四边形替代 P1 三角，同时也是默认 consistent mass
make CXXFLAGS+="-DUSE_Q1_MESH"

# 复现 SESSION 5 之前的（有 bug 的）行为，用于对照实验
make CXXFLAGS+="-DUSE_LUMPED_MASS -DSOURCE_SMOOTH_RADIUS=3 -DSOURCE_SMOOTH_SIGMA=1.7"

# 无 PML 均匀介质做正确性验证：应该看到清晰圆形 P 波
make CXXFLAGS+="-DNBH=0 -DNPMLDEF=0 -DNT=800"
```

**关键教训（SESSION 5）**：
- **velocity-stress 一阶系统必须用 consistent mass**。lumped mass 只对 displacement 二阶
  系统 (M ü = -K u) 是安全的。
- P1 tri 与 Q1 quad 在均匀介质下**行为几乎相同**（都对 = 都错），说明网格类型不是
  波动 FEM 精度的关键。真正影响精度的是：mass 矩阵、时间格式、多项式阶数、staggering。
- 剩余的 P1-P1 数值色散（尾部小振荡）是**低阶格式固有**的，主波前已正确；进一步
  减小需要更高阶元（P2/spectral）或换 CN 时间格式。

---

## 5. 输出与 FDM 对照

FEM 输出文件名与 FDM 完全相同：

- `snap_vx_it####.dat` / `snap_vz_it####.dat`
- `snap_txx_it####.dat` / `snap_tzz_it####.dat` / `snap_txz_it####.dat`
- `source_time.dat`
- `trace_txx.dat` / `trace_tzz.dat` / `trace_pressure.dat`
- `trace_zaxis_txx.dat`（形状 `Nz × NT`）

**尺寸差异需要留意**：

| 场 | FDM 尺寸（交错） | FEM 尺寸（节点） | 是否可直接对比 |
|---|---|---|---|
| σxx | Nz × Nx | Nz × Nx | ✅ 逐点 |
| σzz | Nz × Nx | Nz × Nx | ✅ 逐点 |
| vx  | Nz × (Nx+1) | Nz × Nx | ➖ 位置有半格偏差，可插值后对比 |
| vz  | (Nz+1) × Nx | Nz × Nx | ➖ 同上 |
| σxz | (Nz+1) × (Nx+1) | Nz × Nx | ➖ 同上 |
| trace_txx.dat | 1 × NT | 1 × NT | ✅ 接收点直接对比 |
| trace_zaxis_txx.dat | Nz × NT | Nz × NT | ✅ 沿源 z 轴 σxx 时空图直接对比 |

推荐用 `trace_txx.dat` / `trace_pressure.dat` / `trace_zaxis_txx.dat` 做定量对比。

---

## 6. 面向 CIM 的预留接口（文档 §17）

已到位的架构位：

1. **`LinearSystem { A, b, x }`** — 显式暴露的中间层表示
2. **`ILinearSolver::solve(SolveRequest)`** — 统一后端接口
3. **`DummyCIMSolver`** — 已经出现在工厂里，可作为落位入口
4. **`SolverKind`** 枚举 —— 后续加 `CIM`, `BlockJacobi`, `Schwarz` 等

**下一版**（文档阶段 C-D）计划：
- 引入 `Subdomain { id, elem_ids, interior_dofs, interface_dofs, A_local }`
- `BlockJacobiSolver` / `AdditiveSchwarzSolver`
- 最后把 `A_i x_i = b_i` 交给 `CIMSolver::solve(...)`

---

## 7. 一些数值细节 / 注意事项

1. **消耗**：`Ntot = 5 * 200 * 200 = 200 000` 个未知量。CSR 中 nnz ≈ 200k × 25 ≈ 5M。
   Windows/桌面 g++ -O3 组装耗时约 10–30 s，时间循环每步 ~ 数十次 BiCGSTAB 迭代。
2. **警告**：若使用 `-DNBH=0`，则整域退化为均匀固体，与 FDM 的 `NBH=0` 均匀情形一致。
3. **Backward Euler 会引入数值耗散**（文档 §8.3 明确接受这一点）；后续升级到 Crank–Nicolson 时，只要把 `assembler.cpp` 的时间因子改成 `dt/2`，并在 `assemble_rhs` 中额外加 `-dt/2 K q^n` 即可，其余代码无需变动。
4. **首次跑通建议**：先 `make small`（50×50 × 800 步），看到 `res` 稳定、`sou` 的 Ricker 波形正常再上大算例。
