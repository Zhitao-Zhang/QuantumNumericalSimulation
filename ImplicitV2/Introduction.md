# ImplicitV2：`main.cpp` 技术说明

本文档对应 `Quantum/ImplicitV2/main.cpp`：二维速度–应力形式弹性波方程，采用**交错网格有限差分（SGFD）**空间离散、**Crank–Nicolson（CN）**时间离散，并联合 **ADE 型 PML**；每个时间步将“半隐式离散后的所有未知量”组成一个大型稀疏线性方程组 **`A x = b`**，用 **Jacobi 预处理 BiCGSTAB** 迭代求解。

---

## 1. 项目概览

| 项目 | 说明 |
|------|------|
| 方程 | 2D 各向同性弹性介质中的速度–应力一阶双曲系统 |
| 空间 | 交错网格：`v_x`、`v_z`、`σ_{xx}`、`σ_{zz}`、`σ_{xz}` 位于不同网格点 |
| 时间 | Crank–Nicolson：`n` 与 `n+1` 各取一半 |
| 吸收边界 | PML：用 8 个辅助场 `ψ` 把“拉伸坐标导数”写成 ADE，并与主变量**联立**进 `A x = b` |
| 线性代数 | `A`：CSR 稀疏矩阵（`CSRMatrix`），每时间步组装 `b`，用 `bicgstab` 解 `x` |
| 介质 | 背景为一套固体参数；可选在代码中划定区域使用**第二套固体**参数（实现与文件头部分中文注释所述“含水井流体”不一致，见 §5.3） |

---

## 2. 控制方程（连续形式）

记 `v_x, v_z` 为速度分量，`σ_{xx}, σ_{zz}, σ_{xz}` 为应力，`ρ` 密度，`λ, μ` Lamé 常数，`C = λ+2μ`，`L = λ`，`M = μ`。

### 2.1 动量方程

\[
\rho \frac{\partial v_x}{\partial t} = \frac{\partial \sigma_{xx}}{\partial x} + \frac{\partial \sigma_{xz}}{\partial z},
\qquad
\rho \frac{\partial v_z}{\partial t} = \frac{\partial \sigma_{xz}}{\partial x} + \frac{\partial \sigma_{zz}}{\partial z}.
\]

程序在 `v_x` 位置用 `ρ` 在 `x` 方向半点上的值 `ρ_x`，在 `v_z` 位置用 `ρ` 在 `z` 方向半点上的值 `ρ_z`（见 `build_medium`）。

### 2.2 本构（Hooke）方程

\[
\frac{\partial \sigma_{xx}}{\partial t} = C \frac{\partial v_x}{\partial x} + L \frac{\partial v_z}{\partial z},
\qquad
\frac{\partial \sigma_{zz}}{\partial t} = L \frac{\partial v_x}{\partial x} + C \frac{\partial v_z}{\partial z},
\qquad
\frac{\partial \sigma_{xz}}{\partial t} = M \left( \frac{\partial v_z}{\partial x} + \frac{\partial v_x}{\partial z} \right).
\]

`C, L` 定义在应力整点 `(i_x, i_z)`；`M` 定义在 `σ_{xz}` 的交错点，由周围四格 `μ` **算术平均**（`g_M`）。

### 2.3 PML：ADE 与对导数的替换

在 PML 内，对每个相关空间导数采用“**原导数 + 辅助项**”：

\[
\tilde{\partial}_\alpha = \partial_\alpha + \psi_\alpha.
\]

辅助量满足一阶方程（`d_α` 为阻尼，与位置、方向有关）：

\[
\frac{\partial \psi_\alpha}{\partial t} + d_\alpha \psi_\alpha = - d_\alpha \, (\partial_\alpha \text{ field}).
\]

时间方向用与 `Quantum/Implicit/main.cu` 一致的**指数因子 + 梯形规则**离散，整理为（见源文件顶部注释）：

\[
\psi_\alpha^{n+1} + \frac{d_\alpha \Delta t}{2}\, (\partial_\alpha \text{field})^{n+1}
= E_\alpha \left[ \psi_\alpha^{n} - \frac{d_\alpha \Delta t}{2}\, (\partial_\alpha \text{field})^{n} \right],
\quad E_\alpha = e^{-d_\alpha \Delta t}.
\]

在 `d_α=0` 区域：`E_α=1`，若初值 `ψ≡0` 则 `ψ` 恒为 0，主方程退化为无 PML 版本。

---

## 3. 时间离散：Crank–Nicolson 主方程

定义

\[
c_{ax} = \frac{\Delta t}{2\,\Delta x}, \qquad c_{az} = \frac{\Delta t}{2\,\Delta z}.
\]

对任一量 `f`，CN 写法为：时间导数一侧在 `n` 与 `n+1` 上各取一半。以动量方程为例，加入 PML 后左端（未知在 `n+1`）形如（与源文件注释一致，略去下标网格）：

\[
\rho_x v_x^{n+1}
- \frac{\Delta t}{2}\bigl( \partial_x \sigma_{xx} + \psi_{x,\sigma xx} \bigr)^{n+1}
- \frac{\Delta t}{2}\bigl( \partial_z \sigma_{xz} + \psi_{z,\sigma xz} \bigr)^{n+1}
=
\rho_x v_x^{n}
+ \frac{\Delta t}{2}\bigl( \cdots \bigr)^{n}.
\]

应力方程同理，`C,L,M` 乘在相应 `∂v` 与 `ψ` 项上。

**代码中的对应关系**

- `assemble_A`：把上式**左端**所有与 `^{n+1}` 有关的项写成 `A` 的列系数（含 `σ`、`v`、各类 `ψ` 的耦合）。
- `assemble_b`：把**右端**全部用 `n` 时刻已知场与 `ψ^n` 算出来。
- 因而每个时间步确实是解联立线性系统 **`A x^{n+1} = b^n`**。

---

## 4. 空间离散：有限差分（交错网格）

程序使用**两点中心差分**逼近一阶空间导数，形式为

\[
(\partial_x f)_{i+1/2} \approx \frac{f_{i+1} - f_i}{\Delta x},
\qquad
(\partial_z f)_{j+1/2} \approx \frac{f_{j+1} - f_j}{\Delta z},
\]

具体哪两个指标参与差分由场分量的交错位置决定（见 `assemble_A` / `assemble_b` 中 `if` 边界条件：`ix ≤ Nx-2`、`ix ≥ 1`、`iz ≥ 1` 等），在域外“缺邻居”时**不写该差分项**（与注释“缺邻居 → 跳过”一致）。

常数引用（与代码一致）：

- `cax = dt/(2*dx)`，`caz = dt/(2*dz)`；
- `hd = 0.5*dt`：主方程中 `ψ` 项的系数前缀（再乘 `C/L/M`）。

---

## 5. PML 阻尼曲线与代码中的采样

### 5.1 连续定义

\[
d_0 = \frac{3 V_0 \ln(1/R_c)}{2 D_p}, \quad D_p = N_p\, h,
\]

其中 `V0_pml = Vp_val`，`RC_pml = 1e-6`，`NP` 为 PML 层数，`h` 为 `dx` 或 `dz`。在左、右吸收层内 `d` 为**二次**分布，区外为 0（`dampx_at` / `dampz_at`）。

### 5.2 五条 1D 采样曲线与 \(E = e^{-d\Delta t}\)

为避免在每个方程里重复调用连续函数，程序预先填充：

- `g_dx_int`, `g_ex_int`：`x = i_x \Delta x`，`i_x = 0..N_x-1`
- `g_dx_half`, `g_ex_half`：`x = (i_x + 1/2)\Delta x`，`i_x = 0..N_x`
- `g_dz_int`, `g_ez_int`：`z = i_z \Delta z`，`i_z = 0..N_z-1`
- `g_dz_hv`, `g_ez_hv`：`z = (i_z + 1/2)\Delta z`（用于部分 `ψ`）
- `g_dz_hxz`, `g_ez_hxz`：`z = (i_z - 1/2)\Delta z`（用于 `σ_{xz}` 相关 `ψ`）

八个 `ψ` 各自绑定到哪条 `(d, E)` 曲线，源文件顶部注释表已有完整对应。

---

## 5.3 介质建模：`build_medium`（实现细节）

1. 在每个应力整点单元 `(iz, ix)` 存储 `rho_cell, lam_cell, mu_cell`，其中  
   \(\mu = \rho V_s^2\)，\(\lambda + 2\mu = \rho V_p^2\)，\(\lambda = (\lambda+2\mu) - 2\mu\)。
2. **异常区域**：当 `NBH > 0` 且 `100 ≤ iz < 120` 时，使用第二套固体 `Vp_val2, Vs_val2, rho_val2`；否则使用默认固体 `Vp_val, Vs_val, rho_val`。  
   **注意**：此处并未把井内置为 `Vp_fluid` 等流体参数；文件中打印的流体参数主要用于参考或与别处版本对齐，`g_*` 数组实际由上述两套固体填充。
3. `g_rho_x`：`v_x` 位置，`ρ` 为左右两单元算术平均。  
4. `g_rho_z`：`v_z` 位置，`ρ` 为上下两单元算术平均。  
5. `g_C`, `g_L`：应力整点直接取单元值。  
6. `g_M`：`σ_{xz}` 点四周四个单元的 `μ` 算术平均（边界外索引 `clamp`）。

---

## 6. 未知向量 `x` 的布局（13 块）

块顺序与偏移在 `main` 开头计算：

| 块索引 | 物理量 | 尺寸 |
|--------|--------|------|
| 0 | `v_x` | `Nz × (Nx+1)` |
| 1 | `v_z` | `(Nz+1) × Nx` |
| 2 | `σ_{xx}` | `Nz × Nx` |
| 3 | `σ_{zz}` | `Nz × Nx` |
| 4 | `σ_{xz}` | `(Nz+1) × (Nx+1)` |
| 5–12 | 八个 `ψ`（`pxSxx`, `pzSxz`, …） | 与各自主变量同形 |

全局索引函数：`id_vx`, `id_vz`, `id_xx`, `id_zz`, `id_xz`, `id_pxSxx`, …（均等于 **块偏移 + 二维线性下标**）。  
仅供裸向量用的本地线性下标：`lx_vx`, `lx_vz`, …。

总未知数：`Ntot = off_pzVx + Nxz`。

---

## 7. 是否真的写成 \(A x = b\)？

**是的。** 论证如下：

1. **`assemble_A`** 只依赖网格、介质 `g_*`、PML 曲线 `g_d*`（及固定 `dt,dx,dz`），**不含**当前时刻波场数值 → `A` **与时间步无关**，只需组装一次。
2. 对每个网格点的每条离散方程，`assemble_A` 向 `rows[row]` 的 `std::map<int,double>` 累加：**该行未知量的线性组合系数**，对应 \(A_{row,col}\)。
3. **`assemble_b`** 使用 **\(n\) 时刻** 的 `vx_old, vz_old, σ_old, ψ_old` 及震源项，计算出每条方程的右端 **`b[row]`**，其中：
   - 主方程右端 = CN 右端形式（含 `+ (Δt/2)·ψ^n` 或 `+ (Δt·弹性模量/2)·ψ^n`）；
   - `ψ` 行右端 = \(E\cdot(\psi^n - \frac{d\Delta t}{2}(\partial_\alpha \text{field})^n)\)，其中旧时刻空间导数用**已知场**的有限差分算出。
4. **`CSRMatrix::build_from_rows`** 把 `rows` 转为 CSR；**`bicgstab`** 解 \(A x = b\)（\(x\) 为 \(n+1\) 时刻全场）。
5. **`unpack_solution`** 把 \(x\) 写回 13 个二维数组。

因此：**有限差分 + CN + ADE-PML 离散后，每一步确实是稀疏线性代数意义下的 `A x = b`，而非算子分裂显式推进。**  
矩阵 \(A\) 一般**不对称**（BiCGSTAB 适用于非对称方程）。

---

## 8. 结构体与函数说明（逐函数）

### 8.1 `dampx_at(x)` / `dampz_at(z)`

给定物理坐标，返回 PML 阻尼 `d(x)` 或 `d(z)`（二次剖面，`NP≤0` 时恒为 0）。

### 8.2 `build_damping_profiles()`

遍历离散位置填充 `g_dx_*`, `g_dz_*` 及对应的 `g_ex_* = exp(-d·Δt)`。

### 8.3 `build_medium()`

构造 `g_rho_x`, `g_rho_z`, `g_C`, `g_L`, `g_M`（§5.3）。

### 8.4 `CSRMatrix`

| 成员 | 作用 |
|------|------|
| `build_from_rows(n, rows)` | 由 `vector<map<int,double>>` 生成 CSR；提取对角线 `diag`；若对角为 0 则置 1（避免预处理除零） |
| `multiply(x, y)` | \(y = A x\) |

### 8.5 `bicgstab(A, b, x, maxiter, tol, final_res)`

带 **Jacobi 预处理**（用 `diag[i]`）的双共轭梯度稳定化方法迭代求解 \(A x = b\)。

- 初值 `x` 由调用方提供；时间循环里用上一步全场作为 **warm start**。
- 返回迭代次数；相对残差 `\|r\|/\|b\|` 与 `tol` 比较。
- `rho` 过小等病态情形有限次重置搜索方向。

### 8.6 `assemble_A(rows)`

核心：**组装 \(A\)**。

- **(A) 五条主方程行**：对每个 `(iz, ix)`，`rows[row]` 累加：
  - `v`：`ρ` 在对角；
  - `∂σ`：`±cax` 或 `±caz` 作用在相邻 `σ` 未知数上（边界处按代码条件跳过）；
  - `ψ`：`−hd` 或 `−hd·C/L/M` 作用在对应的 `ψ` 未知数上。
- **(B) 八条 `ψ` 方程行**：对角线 `+1`；若该处 `d>0` 且差分合法，则把 \(\frac{d\Delta t}{2}(\partial_\alpha \text{field})^{n+1}\) 展开为对 \(σ\) 或 \(v\) 的 **±系数**，系数为 `(half_dt * d)/dx` 或 `/dz`。

边界合法条件必须与主方程中同一导数完全一致（源文件注释强调否则会破坏闭合性）。

### 8.7 `assemble_b(b, …_old, sou_t)`

核心：**组装 \(b^n\)**。

- 主方程：逐项镜像 `assemble_A` 的 CN 右端（用 `_old` 场计算差分），并加上 `+ hd * ψ_old`（或乘弹性模量）。
- 震源：`b[id_xx(jsou,isou)] += dt * sou_t`，`σ_{zz}` 同理（应力张量各向同性膨胀型点源）。
- 八个 `ψ`：`b = E * (ψ_old - half_dt * d * (∂field)_old)`，其中 `(∂field)_old` 与 `assemble_A` 中模板一致。

### 8.8 `unpack_solution(x, …)`

按 §6 块偏移把全局向量 \(x\) 拷贝回各 `std::vector<double>`。

### 8.9 `wfile2d` / `wfile2d_d`

将二维场写成 ASCII 文本（行优先，`%e`）；`wfile2d_d` 先转为 `float` 再调用 `wfile2d`。

### 8.10 `main`

1. 计算各块 `Nvx`…、`off_*`、`Ntot`。  
2. `build_medium()`、`build_damping_profiles()`。  
3. 分配 13 个场向量、`b`、`xsol`。  
4. **一次性** `assemble_A` → `A.build_from_rows`。可选 `Nz==Nx==4` 时打印第一行验证。  
5. 时间循环 `it = 0 .. NT-1`：  
   - Ricker：`a = π F_0 (t - T_0)`，`sou_t = -(1 - 2a²)exp(-a²)`（导数形式因子）。  
   - 将上一时刻全场写入 `xsol` 作初值；`assemble_b`；`bicgstab`；`unpack_solution`。  
   - 记录道集、快照、`trace_zaxis_txx`。  
6. 写出 `source_time.dat`、`trace_*.dat`、`trace_zaxis_txx.dat`。

---

## 9. 震源（Ricker）

主频 `F0`，`T0 = 1.2/F0`，时间 \(t = i_t \Delta t\)：

\[
a = \pi F_0 (t - T_0), \qquad
s(t) = -(1 - 2 a^2)\, e^{-a^2}.
\]

加到 \(\sigma_{xx}\) 与 \(\sigma_{zz}\) 的右端：**各加 `dt * sou_t`**（即时间积分意义下的应力源项缩放）。

---

## 10. 公式汇总清单

以下为文档中出现的离散与物理公式索引：

| 编号 | 公式含义 |
|------|----------|
| P1 | 动量方程（\(ρ \partial_t v_x\)，\(ρ \partial_t v_z\)） |
| P2 | Hooke 方程（\(σ_{xx}, σ_{zz}, σ_{xz}\)） |
| P3 | PML：\(\tilde\partial_\alpha = \partial_\alpha + \psi_\alpha\) |
| P4 | ADE：\(\partial_t ψ_\alpha + d_\alpha ψ_\alpha = -d_\alpha \partial_\alpha \text{field}\) |
| P5 | \(ψ\) 的 CN 离散（含 \(E_\alpha = e^{-d_\alpha \Delta t}\)） |
| P6 | CN 系数 \(c_{ax}=\Delta t/(2\Delta x)\)，\(c_{az}=\Delta t/(2\Delta z)\) |
| P7 | PML：`d_0` 与二次剖面 |
| P8 | 弹性参数：\(\mu=\rho V_s^2\)，\(\lambda+2\mu=\rho V_p^2\) |
| P9 | 交错网格中心差分模板 |
| P10 | Ricker：`s(t)=-(1-2a^2)e^{-a^2}` |

---

## 11. 输出文件说明

| 文件模式 | 含义 |
|----------|------|
| `snap_*_itXXXX.dat` | 每 50 步输出 \(v_x, v_z, σ_{xx}, σ_{zz}, σ_{xz}\) 快照 |
| `source_time.dat` | `1 × NT`，震源时间序列 |
| `trace_txx.dat` / `trace_tzz.dat` / `trace_pressure.dat` | 接收点 `(jrec, irec)` 上的应力道，`p = (σ_{xx}+σ_{zz})/2` |
| `trace_zaxis_txx.dat` | `Nz × NT`，\(x=i_{sou}\) 一整列 \(\sigma_{xx}\) 的 z–t 记录 |

---

## 12. 编译与尺度提示

- 网格尺度、`NT`、`NP`、`NBH`、`DX_VAL`、`DZ_VAL`、`DT_VAL`、`F0_HZ` 等可通过 `-D` 宏在编译期覆盖（见 `main.cpp` 顶部 `#ifndef`）。
- **矩阵规模**：未知数约 \(13\) 倍应力网格量级；`NNZ` 随 `Nz,Nx` 增长，大步长网格时注意内存与时间。

---

*文档生成依据：`Quantum/ImplicitV2/main.cpp` 源码结构与注释；若注释与 `build_medium` 中介质赋值不一致，以 §5.3 所述代码路径为准。*
