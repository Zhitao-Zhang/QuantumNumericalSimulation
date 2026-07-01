# CIM 相干伊辛机模拟与地震速度反演 — 项目总结

## 目录

1. [项目概览](#1-项目概览)
2. [物理原理](#2-物理原理)
3. [CIM 求解算法](#3-cim-求解算法)
4. [地震反演应用](#4-地震反演应用)
5. [与 FDM 正演程序的关系](#5-与-fdm-正演程序的关系)
6. [代码结构与模块分析](#6-代码结构与模块分析)
7. [Benchmark 测试框架](#7-benchmark-测试框架)
8. [关键公式推导](#8-关键公式推导)

---

## 1. 项目概览

本项目实现了一个 **相干伊辛机 (Coherent Ising Machine, CIM)** 的数值模拟器，并将其应用于**地震速度反演**问题。核心思想是：

1. 将地下速度模型的反演问题**映射**为 Ising 自旋优化 (QUBO) 问题
2. 利用 CIM 模拟器**求解**该组合优化问题，获得最优速度模型
3. 与穷举搜索、贪心搜索、模拟退火等经典方法**对比**性能

### 技术栈

| 组件 | 技术 |
|------|------|
| 编程语言 | Python 3 |
| 核心依赖 | NumPy ≥1.21, PyTorch ≥1.12, SciPy ≥1.7, Matplotlib ≥3.5 |
| 计算后端 | PyTorch (支持 CPU/GPU 批量并行) |

### 目录结构

```
CIM/
├── requirements.txt          # Python 依赖
├── cimsim/                   # 自主开发的 CIM 求解器核心
│   ├── ising.py              # Ising 问题定义与能量计算
│   ├── ahc.py                # AHC 算法实现
│   ├── cac.py                # CAC 算法实现
│   └── solver.py             # 多配置集成求解器
├── baseline/                 # 参考文献中的基线实现
│   └── cim_optimizer/
│       ├── AHC.py            # 原版 AHC (Leleu et al. 2019)
│       ├── CAC.py            # 原版 CAC (Leleu et al. 2021)
│       ├── extAHC.py         # 扩展 AHC (含外场 h)
│       └── CIM_helper.py     # 工具函数 (文件加载, 暴力搜索)
├── benchmark/                # 性能评估框架
│   ├── instances.py          # 测试实例生成器
│   ├── evaluate.py           # 对比评估框架
│   └── scientific.py         # 科学计算 Benchmark (5 个领域)
├── instances/                # 预生成的测试实例 (.npz)
│   ├── random_N{20..1000}_d*.npz
│   ├── maxcut_N{20..500}_d*.npz
│   └── sk_N{50..500}.npz
└── seismic/                  # 地震反演应用
    ├── forward.py            # 2D 声波正演模拟器
    ├── inversion.py          # QUBO 构建与反演求解
    ├── run_inversion.py      # 阻挫反演对比实验
    └── run_scale.py          # 规模缩放性能实验
```

---

## 2. 物理原理

### 2.1 Ising 模型

Ising 模型是统计物理中描述自旋系统的经典模型。系统由 \(N\) 个自旋变量 \(s_i \in \{-1, +1\}\) 组成，其能量 (哈密顿量) 定义为：

$$
H(\mathbf{s}) = -\frac{1}{2} \mathbf{s}^T J \mathbf{s} - \mathbf{h}^T \mathbf{s}
$$

其中：
- \(J\) 是 \(N \times N\) 对称耦合矩阵，\(J_{ij}\) 描述自旋 \(i\) 和 \(j\) 之间的相互作用
- \(\mathbf{h}\) 是外部磁场向量
- \(J_{ij} > 0\)：铁磁耦合（倾向于自旋对齐）
- \(J_{ij} < 0\)：反铁磁耦合（倾向于自旋反向）

**求解 Ising 模型的基态**（最小化 \(H\)）是一个 NP-hard 的组合优化问题。

### 2.2 相干伊辛机 (CIM) 原理

CIM 是一种基于光学参量振荡器 (OPO) 网络的物理计算系统。其核心思想是：

1. **每个 OPO 脉冲对应一个 Ising 自旋**：脉冲的相位 (0 或 π) 对应 \(s_i = +1\) 或 \(-1\)
2. **OPO 振幅 \(x_i\) 是连续变量**：通过参数泵浦从零态经历分岔，最终二值化为 \(\pm 1\)
3. **自旋间耦合通过反馈实现**：计算 \(J \cdot \mathbf{x}\) 并反馈回 OPO，驱动系统演化到低能态

数值模拟通过常微分方程 (ODE) 迭代求解：

$$
\frac{dx_i}{dt} = (p - 1)x_i - \mu x_i^3 + \epsilon \sum_j J_{ij} x_j e_j
$$

其中：
- \(p\)：泵浦参数，控制增益
- \(\mu x_i^3\)：饱和非线性项
- \(\epsilon\)：反馈强度
- \(e_j\)：误差变量（用于振幅校正）

### 2.3 为什么 CIM 能求解组合优化

CIM 的物理演化过程本质上是在连续空间中进行的**全局搜索**：

- 初始状态：所有 \(x_i \approx 0\)（量子噪声级别）
- 泵浦增加 → 参数振荡 → 自旋逐渐分岔到 \(\pm 1\)
- 反馈耦合使系统偏向低能态构型
- 多个 OPO 同时演化 → **天然并行性**
- 批量运行 (batch) → 统计采样多个初始条件

---

## 3. CIM 求解算法

### 3.1 AHC 算法 (Amplitude Heterogeneity Correction)

**文献来源**: T. Leleu, Y. Yamamoto, P.L. McMahon, K. Aihara, *Phys. Rev. Lett.* **122**, 040607 (2019)

AHC 解决了基础 CIM 中振幅不均匀性导致的亚最优问题。核心更新方程：

```
自旋振幅更新:
  x += dt * x * ((r - 1) - μ * x²)              # 泵浦 + 饱和
  x += dt * ε * (x @ J) * e                      # 反馈 (带误差校正)

误差变量更新:
  target_a = a₀ + ε * mean(sig @ J * sig)        # 自适应目标振幅
  e += dt * (-β * (x² - target_a) * e)           # 误差校正
  e = clamp(e, max=32)                            # 防止发散
```

**关键参数**:

| 参数 | 默认值 | 物理含义 |
|------|--------|----------|
| `r` | 0.2 | 泵浦增益率 |
| `β` | 0.05 | 误差变量学习率 |
| `ε` | 0.07 | 反馈耦合强度 |
| `μ` | 1.0 | 饱和非线性系数 |
| `dt` | 0.05 | 时间步长 |
| `a₀` | 0.2 | 目标振幅基线 |

### 3.2 CAC 算法 (Chaotic Amplitude Control)

**文献来源**: T. Leleu, F. Khoyratee, T. Levi, R. Hamerly, T. Kohno, K. Aihara, *Commun. Phys.* **4**, 266 (2021)

CAC 在 AHC 基础上引入**混沌动力学**，利用混沌逃逸机制避免局部最优。关键创新：

```
混沌振幅控制:
  ξ += γ * dt                                     # 阻尼参数缓慢增加
  ΔH = H_current - H_optimal                      # 与最优能量的差异
  a = α + ρ * tanh(δ * ΔH)                        # 自适应目标振幅

误差变量:
  e += dt * (-ξ * (x² - a) * e)                   # ξ 控制校正强度

重置机制:
  当 ΔH < 0 (发现更优解): 重置 t_c, 更新 H_opt
  当 t - t_c > τ/dt (长时间无改善): 重置 t_c
```

**关键参数**:

| 参数 | 默认值 | 物理含义 |
|------|--------|----------|
| `α` | 3.0 | 基线目标振幅 |
| `β` | 0.25 | 反馈强度 |
| `γ` | 1.1e-4 | 阻尼增长率 |
| `δ` | 10.0 | 混沌非线性灵敏度 |
| `ρ` | 3.0 | 混沌控制幅度 |
| `τ` | 1000 | 重置时间窗口 |

### 3.3 多配置集成求解器 (Smart Solver)

`solver.py` 实现了一个集成策略，同时运行多组超参数配置：

- **AHC**: `ε ∈ {0.04, 0.07, 0.12, 0.20}` → 4 组
- **CAC**: `γ ∈ {5e-3, 1e-2, 5e-2}` → 3 组
- 每组使用 `batch_size` 个独立运行
- 总有效试验数 = `batch_size × 7`
- 返回所有试验中能量最低的解

这种集成策略通过超参数多样性补偿了单一配置可能的盲区。

### 3.4 与基线实现的差异

本项目包含两套实现——`baseline/` (原始论文代码) 和 `cimsim/` (优化重写)：

| 特性 | `baseline/` | `cimsim/` |
|------|------------|-----------|
| 外部磁场 h | 不支持 (AHC/CAC) | 完整支持 |
| 内存效率 | 预分配轨迹存储 | 可选轨迹存储 |
| 矩阵-向量乘法 | `torch.bmm` 逐批次 | `torch.mm` 批量并行 |
| Scratch 张量 | 每步重分配 | 预分配复用 |
| 调度器支持 | Callable 函数 | 预构建张量向量化 |
| 扩展算法 | `extAHC` (带 h 场) | 原生集成在核心中 |

---

## 4. 地震反演应用

### 4.1 问题描述

地震速度反演的目标是从**观测到的地震走时数据**推断地下**速度模型**。这是一个经典的**反问题**：

- **正问题**: 给定速度模型 → 计算地震波走时
- **反问题**: 给定观测走时 → 反推速度模型

### 4.2 二值速度层析 (Binary Velocity Tomography)

本项目采用**二值参数化**方案：

- 将地下模型划分为 \(N\) 个区域
- 每个区域 \(i\) 只有两种速度选择：\(v_{\text{slow}}[i]\) 或 \(v_{\text{fast}}[i]\)
- 使用 Ising 自旋编码：\(s_i = -1 \to v_{\text{slow}}\)，\(s_i = +1 \to v_{\text{fast}}\)

### 4.3 无惩罚项 QUBO 公式 (inversion.py)

**慢度参数化**：

$$
u_i(s) = \frac{u_{\text{slow},i} + u_{\text{fast},i}}{2} + \frac{u_{\text{fast},i} - u_{\text{slow},i}}{2} \cdot s_i = \bar{u}_i + \Delta u_i \cdot s_i
$$

其中 \(u = 1/v\) 为慢度。

**走时计算**（直线射线近似）：

$$
T_{\text{syn}} = L \cdot \mathbf{u}(s) = L \cdot \bar{\mathbf{u}} + L_\Delta \cdot \mathbf{s}
$$

其中 \(L\) 是射线路径长度矩阵 (\(n_{\text{rec}} \times n_{\text{regions}}\))，\(L_\Delta\) 的第 \(i\) 列为 \(L_{:,i} \cdot \Delta u_i\)。

**目标函数** (走时残差最小二乘)：

$$
f(\mathbf{s}) = \|T_{\text{obs}} - T_{\text{syn}}\|^2 = \|\mathbf{r} - A\mathbf{s}\|^2
$$

其中 \(\mathbf{r} = T_{\text{obs}} - L\bar{\mathbf{u}}\)，\(A = L_\Delta\)。

**精确映射到 Ising 哈密顿量**：

$$
J_{ij} = -2(A^T A)_{ij} \quad (i \neq j), \qquad h_i = 2(A^T \mathbf{r})_i
$$

这是一个**无惩罚项**的精确映射——不需要任何约束违反惩罚项。

### 4.4 阻挫反演 (Frustrated Binary Tomography, run_inversion.py)

纯走时反演对贪心方法"太容易"。实际反演需要**正则化**来处理不适定性，这引入了**阻挫** (frustration)：

$$
E(\mathbf{s}) = \alpha \|T_{\text{obs}} - T_{\text{syn}}(\mathbf{s})\|^2 + \beta \sum_{\langle ij \rangle} (1 - s_i s_j)
$$

- **数据保真项** (\(\alpha\))：反铁磁耦合 (\(J < 0\))，驱动异质性
- **空间平滑项** (\(\beta\), Potts 模型)：铁磁耦合 (\(J > 0\))，驱动均匀性
- **阻挫 = 两者竞争**：\(J\) 矩阵同时含正负元素，形成典型的 NP-hard 优化问题

映射到 Ising 后：

$$
J_{ij} = -2\alpha (A^T A)_{ij} + \beta \cdot \text{adj}(i,j), \qquad h_i = 2\alpha (A^T \mathbf{r})_i
$$

### 4.5 射线追踪

`compute_ray_lengths()` 使用直线射线近似计算每条射线在各区域的路径长度：

- 将射线离散为 1000 段等距采样点
- 每个采样点定位到所属区域
- 累计每个区域内的总路径长度
- 支持多源多检波器配置

### 4.6 实验设计

`run_inversion.py` 设计了 4 组测试：

| 测试 | 网格 | 自旋数 N | 速度范围 | α | β |
|------|------|---------|---------|---|---|
| TEST 1 | 4×5 | 20 | 2500-4500 m/s | 100 | 0.05 |
| TEST 2 | 6×6 | 36 | 2000-5000 m/s | 50 | 0.08 |
| TEST 3 | 8×8 | 64 | 2500-4500 m/s | 30 | 0.1 |
| TEST 4 | 10×10 | 100 | 2000-5000 m/s | 20 | 0.15 |

每组测试对比 CIM、Greedy (多启动贪心)、SA (模拟退火)、Random (随机搜索) 四种方法。

### 4.7 规模缩放实验 (run_scale.py)

`run_scale.py` 测试 CIM 在不同问题规模下的缩放性能：

- 问题规模：N = 50, 100, 200, 400, 800
- CIM 使用固定参数 (timesteps=2000, batch=30)
- 两种对比模式：
  1. **固定预算**：各方法使用各自默认配置
  2. **等时对比**：给 Greedy/SA 与 CIM 相同的计算时间

计算复杂度分析：
- CIM：O(N² × timesteps × batch) — 矩阵-向量乘法主导
- Greedy：O(N² × sweeps × restarts) — 局部场更新主导
- SA：O(N × steps) — 单自旋翻转

---

## 5. 与 FDM 正演程序的关系

### 5.1 FDM 程序概述 (ImplicitV2/main.cpp)

`/Quantum/FDM/ImplicitV2/main.cpp` 实现了一个 **2D 弹性波方程的 Crank-Nicolson 隐式时间推进有限差分求解器**，带有 ADE-PML 吸收边界。

**物理方程** — 速度-应力弹性波方程：

| 方程 | 内容 |
|------|------|
| \(v_x\) 方程 | \(\rho \partial_t v_x = \partial_x \sigma_{xx} + \partial_z \sigma_{xz}\) |
| \(v_z\) 方程 | \(\rho \partial_t v_z = \partial_x \sigma_{xz} + \partial_z \sigma_{zz}\) |
| \(\sigma_{xx}\) 方程 | \(\partial_t \sigma_{xx} = (\lambda+2\mu)\partial_x v_x + \lambda \partial_z v_z\) |
| \(\sigma_{zz}\) 方程 | \(\partial_t \sigma_{zz} = \lambda \partial_x v_x + (\lambda+2\mu) \partial_z v_z\) |
| \(\sigma_{xz}\) 方程 | \(\partial_t \sigma_{xz} = \mu(\partial_x v_z + \partial_z v_x)\) |

**关键特性**：
- **Crank-Nicolson 隐式格式**：无条件稳定，允许大时间步长
- **ADE-PML 吸收边界**：8 个辅助 ψ 场实现完美匹配层
- **13 块未知向量**：5 个物理场 + 8 个 PML 辅助场
- **网格参数**：Nz=200, Nx=200, dx=dz=10m, dt=1e-4s, NT=4500
- **非均匀介质**：含水井孔模型 (固体 Vp=5000 + 流体 Vp=1500)

### 5.2 两个程序的互补关系

```
┌─────────────────────────────────────────────────────────┐
│                    完整反演工作流                          │
│                                                         │
│  ┌──────────────┐        ┌──────────────────┐          │
│  │   FDM 正演    │        │    CIM 反演       │          │
│  │ (ImplicitV2)  │        │   (seismic/)     │          │
│  │              │        │                  │          │
│  │ 输入: 速度模型 │──观测──→│ 输入: 观测走时     │          │
│  │ 输出: 波场/走时│  数据   │ 输出: 反演速度模型  │          │
│  │              │        │                  │          │
│  │ 弹性波方程    │        │ Ising 优化        │          │
│  │ CN 隐式格式   │        │ AHC/CAC 算法      │          │
│  │ ADE-PML 边界  │        │ QUBO 映射         │          │
│  └──────────────┘        └──────────────────┘          │
│                                                         │
│  C++, 高精度正演                Python, 量子启发反演       │
└─────────────────────────────────────────────────────────┘
```

**具体联系**：

1. **FDM 是正演引擎**：给定速度模型计算波场和走时数据
2. **CIM 是反演引擎**：从走时数据反推速度模型
3. **CIM 的 `forward.py` 是简化版正演器**：2D 声波标量方程，2 阶 FD 显式格式，仅用于快速走时计算
4. **FDM 的 `main.cpp` 是高精度正演器**：2D 弹性波矢量方程，CN 隐式格式 + PML，用于精确波场模拟

**协同使用场景**：
- 先用 FDM 程序（或真实地震实验）获取高精度观测数据
- 然后用 CIM 反演程序从观测数据中恢复速度模型
- 反演结果可用 FDM 再次正演验证

### 5.3 参数对应关系

| 参数 | FDM (main.cpp) | CIM (seismic/) |
|------|----------------|----------------|
| 网格步长 | dx=dz=10m | dx=dz=10m |
| 介质速度 | Vp=5000/1500 m/s | v_slow=2000-2500, v_fast=4500-5000 m/s |
| 模型尺寸 | 200×200 | 80-200 (可变) |
| 波方程 | 弹性波 (5 场分量) | 声波 (标量压力) |
| 时间推进 | CN 隐式 (无条件稳定) | 2 阶 FD 显式 (CFL 约束) |
| 边界条件 | ADE-PML (8 辅助场) | 指数衰减吸收层 |

---

## 6. 代码结构与模块分析

### 6.1 cimsim/ising.py — Ising 问题定义

```python
# 核心数据结构
@dataclass
class IsingResult:
    best_spins: np.ndarray          # 最优自旋构型 (N,)
    best_energy: float              # 最低 Ising 能量
    spin_trajectories: np.ndarray   # 自旋振幅演化 (batch, N, T)
    energy_evolution: np.ndarray    # 能量演化 (batch, T)
    all_best_spins: np.ndarray      # 各批次最优 (batch, N)
    all_best_energies: np.ndarray   # 各批次能量 (batch,)
    wall_time: float                # 计算耗时 (秒)

# 能量计算
def ising_energy(J, spins, h=None):
    """H = -1/2 * s^T J s - h^T s"""

def ising_energy_batch(J, sig, h=None):
    """批量能量计算 (PyTorch)"""
```

### 6.2 cimsim/ahc.py — AHC 求解器

完整的 AHC 迭代流程：

```
初始化: x ~ Uniform(-5e-4, 5e-4), e = 1
for t = 0 to ticks-1:
    1. 二值化: sig = sign(x)
    2. 计算 Ising 能量
    3. 矩阵-向量乘: MVM = x @ J
    4. 泵浦 + 饱和: x += dt * x * ((r-1) - μ*x²)
    5. 反馈:         x += dt * ε * MVM * e
    6. 噪声 (可选):  x += ε * noise * rand
    7. 目标振幅:     a = a₀ + ε * mean(sJ * sig)
    8. 误差更新:     e += dt * (-β * (x² - a) * e)
    9. 追踪最优
返回 IsingResult
```

### 6.3 cimsim/cac.py — CAC 求解器

CAC 在 AHC 基础上增加混沌控制：

```
额外状态: ξ=0 (阻尼参数), a=α (目标振幅), t_c=0 (最优时刻)

主循环中额外步骤:
    混沌控制:
      ξ += γ * dt
      ΔH = H - H_opt
      a = α + ρ * tanh(δ * ΔH)
    
    重置:
      if t - t_c > τ/dt: 重置 t_c
    
    更新最优:
      if H < H_opt: 记录并重置
```

### 6.4 seismic/forward.py — 声波正演器

2D 声波方程有限差分求解器，用于地震反演中的快速正演计算：

```python
# 标量波方程: 1/c²·∂²p/∂t² = ∂²p/∂x² + ∂²p/∂z² + s(t)δ(x-xs,z-zs)
# 2 阶 FD 空间差分, 2 阶 FD 时间推进
# 简单指数衰减吸收边界

# 主要类:
class AcousticModel:    # 速度模型 (vp, dx, dz, nx, nz)
class Source:           # 震源 (位置, 子波, dt)
class Receiver:         # 检波器 (位置)

def ricker_wavelet(f0, dt, nt, t0=None):   # Ricker 子波生成
def solve_acoustic_2d(model, source, receivers, nt, abc_width=20):  # 正演求解
```

### 6.5 seismic/inversion.py — QUBO 构建与反演

核心函数链：

```
compute_ray_lengths()      # 射线追踪 → L 矩阵
build_binary_qubo()        # 构建 Ising: J, h = f(L, T_obs, v_slow, v_fast)
cim_inversion()            # CIM 求解 → 最优速度模型
decode_spins()             # s → 速度
compute_travel_time()      # L @ (1/v) → 走时

# 对比求解器:
exhaustive_search()        # 穷举 2^N (验证用)
random_search()            # 随机采样
greedy_search()            # 多启动贪心坐标下降
simulated_annealing()      # 模拟退火
```

### 6.6 seismic/run_inversion.py — 阻挫反演实验

完整的实验框架：

```python
# 核心函数
make_grid_regions()        # 创建区域划分 (n_vert × n_horiz 网格)
get_adjacency()            # 计算 4 连通邻接关系
build_frustrated_ising()   # 构建阻挫 Ising: J_data + J_smooth
compute_objective()        # 数据项 + 平滑项的总目标函数

# 实验流程
run_test():
  1. 生成真实模型 (随机二值)
  2. 设置观测几何 (多源多检)
  3. 计算射线长度矩阵
  4. 生成观测走时 (+ 可选噪声)
  5. 构建阻挫 Ising 哈密顿量
  6. CIM / Greedy / SA / Random 四法求解
  7. 对比 Ising 能量、目标函数、正确率、耗时
```

---

## 7. Benchmark 测试框架

### 7.1 测试实例类型 (instances.py)

| 类型 | 说明 | 尺寸范围 |
|------|------|---------|
| Random Ising | 均匀分布 \(J_{ij} \in [-1,1]\) | N=20~1000 |
| MAX-CUT | 随机图最大割 | N=20~500 |
| SK 模型 | Sherrington-Kirkpatrick 全连接 \(J_{ij} \sim N(0, 1/\sqrt{N})\) | N=50~500 |

### 7.2 科学计算 Benchmark (scientific.py)

覆盖 5 个科学领域：

| 领域 | 问题 | 特点 |
|------|------|------|
| 凝聚态物理 | 2D Edwards-Anderson 自旋玻璃 | ±1 随机最近邻耦合 |
| 图论 | Petersen/Cube/Dodecahedron 最大割 | 经典图结构 |
| 组合优化 | 数划分问题 | 全连接 \(J=-2a_ia_j\) |
| 统计力学 | 三角格子反铁磁 | 几何阻挫 |
| 运筹学 | 加权 MAX-2-SAT | 混合正负耦合 |

### 7.3 评估框架 (evaluate.py)

对比 5 种求解器：

| 求解器 | 代码 | 描述 |
|--------|------|------|
| baseline_CAC | `baseline/` | 原始 CAC 论文实现 |
| baseline_AHC | `baseline/` | 原始 AHC 论文实现 |
| ours_CAC | `cimsim/` | 优化 CAC 实现 |
| ours_AHC | `cimsim/` | 优化 AHC 实现 |
| ours_SMART | `cimsim/` | 多配置集成求解器 |

---

## 8. 关键公式推导

### 8.1 走时反演 → Ising 映射的完整推导

**步骤 1**: 慢度参数化

$$
u_i(s_i) = \bar{u}_i + \Delta u_i \cdot s_i
$$

其中 \(\bar{u}_i = (u_{\text{slow},i} + u_{\text{fast},i})/2\)，\(\Delta u_i = (u_{\text{fast},i} - u_{\text{slow},i})/2\)

**步骤 2**: 走时线性化

$$
T_{\text{syn}} = L \mathbf{u}(\mathbf{s}) = L\bar{\mathbf{u}} + A\mathbf{s}
$$

其中 \(A_{ki} = L_{ki} \cdot \Delta u_i\)

**步骤 3**: 最小二乘展开

$$
f(\mathbf{s}) = \|\mathbf{r} - A\mathbf{s}\|^2 = \mathbf{r}^T\mathbf{r} - 2\mathbf{r}^TA\mathbf{s} + \mathbf{s}^TA^TA\mathbf{s}
$$

**步骤 4**: 利用 \(s_i^2 = 1\) 的约束分离对角项

$$
\mathbf{s}^T(A^TA)\mathbf{s} = \sum_{i \neq j}(A^TA)_{ij}s_is_j + \sum_i(A^TA)_{ii} = \sum_{i \neq j}(A^TA)_{ij}s_is_j + \text{const}
$$

**步骤 5**: 与 Ising 哈密顿量匹配

$$
H = -\frac{1}{2}\sum_{i \neq j}J_{ij}s_is_j - \sum_i h_i s_i
$$

对比得：

$$
-\frac{1}{2}J_{ij} = (A^TA)_{ij} \implies J_{ij} = -2(A^TA)_{ij}, \quad i \neq j
$$

$$
-h_i = -2(A^T\mathbf{r})_i \implies h_i = 2(A^T\mathbf{r})_i
$$

### 8.2 阻挫机制

加入平滑正则化后：

$$
J_{ij}^{\text{total}} = \underbrace{-2\alpha(A^TA)_{ij}}_{\text{反铁磁 (J<0)}} + \underbrace{\beta \cdot \text{adj}(i,j)}_{\text{铁磁 (J>0)}}
$$

- 数据项要求异质性（邻居自旋反向 → 速度对比度 → 走时差异）
- 平滑项要求均匀性（邻居自旋同向 → 空间连续性）
- **两者竞争 = 阻挫 = NP-hard** → CIM 的优势场景

---

## 附录：运行指南

```bash
# 安装依赖
cd /path/to/CIM
pip install -r requirements.txt

# 生成 benchmark 实例
python benchmark/instances.py

# 运行 CIM vs 基线对比
python benchmark/evaluate.py

# 运行科学计算 benchmark
python benchmark/scientific.py

# 运行地震反演实验 (阻挫版)
python seismic/run_inversion.py

# 运行规模缩放实验
python seismic/run_scale.py
```
