# CIM-Based Seismic Velocity Inversion

## Overview

将相干伊辛机（CIM）应用于地震速度反演问题。通过将地下速度模型参数化为 Ising 变量，将反演问题映射为 QUBO 优化。

## 方法

### 物理模型

- **正演器**: 2D 声波方程有限差分求解器 (`forward.py`)
- **反演目标**: 旅行时层析成像 + 空间平滑正则化 (Potts model)

### QUBO 公式（无惩罚项）

每个区域 i 选择快/慢两种速度之一：

```
s_i = +1 → v_fast[i]
s_i = -1 → v_slow[i]
```

目标函数:

```
E(s) = α||T_obs - T_syn(s)||² + β Σ_{<ij>} (1 - s_i·s_j)
```

精确映射为 Ising Hamiltonian H = -(1/2)s^T J s - h^T s:

- **数据项** → 反铁磁耦合 (J < 0)
- **平滑项** → 铁磁耦合 (J > 0)
- **阻挫** = 两者竞争 → NP-hard

### 文件结构

- `forward.py`: 2D 声波正演器 (显式有限差分)
- `inversion.py`: QUBO 构建、射线追踪、各种求解器
- `run_inversion.py`: 完整对比实验

## 实验结果

| 问题规模 | CIM gap | Greedy gap | SA gap | Random gap |
|---------|---------|-----------|--------|-----------|
| N=20    | -0.5%   | -2.6%     | -2.6%  | +1.5%     |
| N=36    | +5.7%   | -2.2%     | -2.4%  | +9.6%     |
| N=64    | -19.8%  | -30.3%    | -30.7% | +11.3%    |
| N=100   | -19.7%  | -33.0%    | -33.0% | +22.0%    |

(负数 = 找到比参考解更低的能量)

## 分析

1. **CIM 优于随机搜索**: 在所有测例中 CIM 显著优于随机搜索
2. **局部搜索更适合稀疏问题**: Greedy/SA 在此类网格结构问题上更高效，因为 J 的连通性以局部邻居为主
3. **CIM 优势场景**: 稠密耦合系统（如 SK 模型、稠密 MAX-CUT），这是已有 benchmark 证实的

## 运行

```bash
cd /path/to/CIMSim
python seismic/run_inversion.py
```
