# README: 2D 速度–应力弹性波方程的 Crank–Nicolson 隐式交错网格有限差分串行代码实现说明

本文档用于指导 AI 或开发者编写一个**串行版本**的二维弹性波模拟程序。程序要求使用**速度–应力一阶弹性波方程**，空间上采用**交错网格有限差分法**，时间上采用 **Crank–Nicolson 隐式格式**，并在每一个时间步显式构建线性系统：

\[
A x = b
\]

然后通过直接线性求解器或基础稠密/稀疏求解器求出新时间层的五个未知场：

\[
v_x^{n+1},\quad v_z^{n+1},\quad \sigma_{xx}^{n+1},\quad \sigma_{zz}^{n+1},\quad \sigma_{xz}^{n+1}.
\]

本 README 的目标是让 AI 在写代码时不要写成显式 leapfrog，也不要写成 Picard 近似，而是**真正组装 Crank–Nicolson 隐式格式对应的 \(Ax=b\) 线性系统**。

---

## 1. 程序目标

编写一个二维弹性波模拟程序，满足以下要求：

1. 求解二维 \(x-z\) 平面内的各向同性弹性波速度–应力方程；
2. 使用一阶速度–应力形式，而不是二阶位移形式；
3. 空间采用交错网格有限差分；
4. 时间采用 Crank–Nicolson 隐式格式；
5. 每一个时间步组装一次线性系统 \(Ax=b\)，然后求解；
6. 初始版本写成串行 CPU 代码；
7. 暂时不设置 PML，不设置吸收边界；
8. 声源采用应力点源，同时加到 \(\sigma_{xx}\) 和 \(\sigma_{zz}\) 上；
9. 介质为均匀介质：

\[
V_p=6000,\qquad V_s=3000,\qquad \rho=2700.
\]

10. 声源函数记为 `sou(t)`，代码中可先使用 Ricker 子波。

---

## 2. 连续二维速度–应力弹性波方程

二维 \(x-z\) 平面内，各向同性弹性介质的速度–应力方程为：

\[
\rho \frac{\partial v_x}{\partial t}
=
\frac{\partial \sigma_{xx}}{\partial x}
+
\frac{\partial \sigma_{xz}}{\partial z},
\]

\[
\rho \frac{\partial v_z}{\partial t}
=
\frac{\partial \sigma_{xz}}{\partial x}
+
\frac{\partial \sigma_{zz}}{\partial z},
\]

\[
\frac{\partial \sigma_{xx}}{\partial t}
=
(\lambda+2\mu)\frac{\partial v_x}{\partial x}
+
\lambda\frac{\partial v_z}{\partial z},
\]

\[
\frac{\partial \sigma_{zz}}{\partial t}
=
\lambda\frac{\partial v_x}{\partial x}
+
(\lambda+2\mu)\frac{\partial v_z}{\partial z},
\]

\[
\frac{\partial \sigma_{xz}}{\partial t}
=
\mu
\left(
\frac{\partial v_x}{\partial z}
+
\frac{\partial v_z}{\partial x}
\right).
\]

Lamé 参数为：

\[
\mu=\rho V_s^2,
\]

\[
\lambda=\rho(V_p^2-2V_s^2),
\]

\[
\lambda+2\mu=\rho V_p^2.
\]

对于本程序的均匀介质：

\[
\mu=2700\times 3000^2,
\]

\[
\lambda=2700\times(6000^2-2\times3000^2),
\]

\[
\lambda+2\mu=2700\times6000^2.
\]

---

## 3. 网格定义与交错网格变量位置

设整点应力网格大小为：

\[
N_z\times N_x.
\]

其中行方向为 \(z\)，列方向为 \(x\)。

整点坐标为：

\[
x_i=i\Delta x,\qquad z_j=j\Delta z.
\]

五个变量的位置如下：

| 变量 | 物理位置 | 数组尺寸 |
|---|---|---|
| \(\sigma_{xx}\) | \((i,j)\) | \(N_z\times N_x\) |
| \(\sigma_{zz}\) | \((i,j)\) | \(N_z\times N_x\) |
| \(v_x\) | \((i+\frac12,j)\) | \(N_z\times (N_x+1)\) |
| \(v_z\) | \((i,j+\frac12)\) | \((N_z+1)\times N_x\) |
| \(\sigma_{xz}\) | \((i+\frac12,j+\frac12)\) | \((N_z+1)\times (N_x+1)\) |

在代码中建议使用如下数组名：

```c
vx[Nz][Nx+1]
vz[Nz+1][Nx]
txx[Nz][Nx]
tzz[Nz][Nx]
txz[Nz+1][Nx+1]
```

注意：数组第一维是 \(z\)，第二维是 \(x\)。如果用一维数组存储，应采用行优先展开：

```c
id = ix + iz * nx_local;
```

---

## 4. Crank–Nicolson 时间离散思想

对一般方程：

\[
\frac{\partial u}{\partial t}=F(u),
\]

Crank–Nicolson 格式为：

\[
\frac{u^{n+1}-u^n}{\Delta t}
=
\frac12\left[F(u^{n+1})+F(u^n)\right].
\]

对于速度–应力弹性波系统，这意味着：

1. 更新速度时，应力导数取 \(n\) 和 \(n+1\) 两个时间层的平均；
2. 更新应力时，速度导数取 \(n\) 和 \(n+1\) 两个时间层的平均；
3. 因此新时间层的速度和应力互相耦合，必须联立求解；
4. 不能先单独更新速度，再更新应力；
5. 必须构建整体线性系统 \(Ax=b\)。

---

## 5. Crank–Nicolson 交错网格空间差分方程

下面给出代码必须实现的离散方程。注意：这些方程中的左端含有新时间层未知量，右端全是旧时间层已知量。

### 5.1 \(v_x\) 方程

\[
\begin{aligned}
\rho_x v_x^{n+1}(i+\frac12,j)
&-
\frac{\Delta t}{2}
\frac{
\sigma_{xx}^{n+1}(i+1,j)-\sigma_{xx}^{n+1}(i,j)
}{\Delta x}
\\
&-
\frac{\Delta t}{2}
\frac{
\sigma_{xz}^{n+1}(i+\frac12,j+\frac12)
-
\sigma_{xz}^{n+1}(i+\frac12,j-\frac12)
}{\Delta z}
\\
&=
\rho_x v_x^{n}(i+\frac12,j)
+
\frac{\Delta t}{2}
\frac{
\sigma_{xx}^{n}(i+1,j)-\sigma_{xx}^{n}(i,j)
}{\Delta x}
\\
&+
\frac{\Delta t}{2}
\frac{
\sigma_{xz}^{n}(i+\frac12,j+\frac12)
-
\sigma_{xz}^{n}(i+\frac12,j-\frac12)
}{\Delta z}.
\end{aligned}
\]

### 5.2 \(v_z\) 方程

\[
\begin{aligned}
\rho_z v_z^{n+1}(i,j+\frac12)
&-
\frac{\Delta t}{2}
\frac{
\sigma_{xz}^{n+1}(i+\frac12,j+\frac12)
-
\sigma_{xz}^{n+1}(i-\frac12,j+\frac12)
}{\Delta x}
\\
&-
\frac{\Delta t}{2}
\frac{
\sigma_{zz}^{n+1}(i,j+1)
-
\sigma_{zz}^{n+1}(i,j)
}{\Delta z}
\\
&=
\rho_z v_z^{n}(i,j+\frac12)
+
\frac{\Delta t}{2}
\frac{
\sigma_{xz}^{n}(i+\frac12,j+\frac12)
-
\sigma_{xz}^{n}(i-\frac12,j+\frac12)
}{\Delta x}
\\
&+
\frac{\Delta t}{2}
\frac{
\sigma_{zz}^{n}(i,j+1)-\sigma_{zz}^{n}(i,j)
}{\Delta z}.
\end{aligned}
\]

### 5.3 \(\sigma_{xx}\) 方程

\[
\begin{aligned}
\sigma_{xx}^{n+1}(i,j)
&-
\frac{\Delta t}{2}
(\lambda+2\mu)_{i,j}
\frac{
v_x^{n+1}(i+\frac12,j)
-
v_x^{n+1}(i-\frac12,j)
}{\Delta x}
\\
&-
\frac{\Delta t}{2}
\lambda_{i,j}
\frac{
v_z^{n+1}(i,j+\frac12)
-
v_z^{n+1}(i,j-\frac12)
}{\Delta z}
\\
&=
\sigma_{xx}^{n}(i,j)
+
\frac{\Delta t}{2}
(\lambda+2\mu)_{i,j}
\frac{
v_x^{n}(i+\frac12,j)
-
v_x^{n}(i-\frac12,j)
}{\Delta x}
\\
&+
\frac{\Delta t}{2}
\lambda_{i,j}
\frac{
v_z^{n}(i,j+\frac12)
-
v_z^{n}(i,j-\frac12)
}{\Delta z}.
\end{aligned}
\]

### 5.4 \(\sigma_{zz}\) 方程

\[
\begin{aligned}
\sigma_{zz}^{n+1}(i,j)
&-
\frac{\Delta t}{2}
\lambda_{i,j}
\frac{
v_x^{n+1}(i+\frac12,j)
-
v_x^{n+1}(i-\frac12,j)
}{\Delta x}
\\
&-
\frac{\Delta t}{2}
(\lambda+2\mu)_{i,j}
\frac{
v_z^{n+1}(i,j+\frac12)
-
v_z^{n+1}(i,j-\frac12)
}{\Delta z}
\\
&=
\sigma_{zz}^{n}(i,j)
+
\frac{\Delta t}{2}
\lambda_{i,j}
\frac{
v_x^{n}(i+\frac12,j)
-
v_x^{n}(i-\frac12,j)
}{\Delta x}
\\
&+
\frac{\Delta t}{2}
(\lambda+2\mu)_{i,j}
\frac{
v_z^{n}(i,j+\frac12)
-
v_z^{n}(i,j-\frac12)
}{\Delta z}.
\end{aligned}
\]

### 5.5 \(\sigma_{xz}\) 方程

\[
\begin{aligned}
\sigma_{xz}^{n+1}(i+\frac12,j+\frac12)
&-
\frac{\Delta t}{2}
\mu_{xz}
\frac{
v_z^{n+1}(i+1,j+\frac12)
-
v_z^{n+1}(i,j+\frac12)
}{\Delta x}
\\
&-
\frac{\Delta t}{2}
\mu_{xz}
\frac{
v_x^{n+1}(i+\frac12,j+1)
-
v_x^{n+1}(i+\frac12,j)
}{\Delta z}
\\
&=
\sigma_{xz}^{n}(i+\frac12,j+\frac12)
+
\frac{\Delta t}{2}
\mu_{xz}
\frac{
v_z^{n}(i+1,j+\frac12)
-
v_z^{n}(i,j+\frac12)
}{\Delta x}
\\
&+
\frac{\Delta t}{2}
\mu_{xz}
\frac{
v_x^{n}(i+\frac12,j+1)
-
v_x^{n}(i+\frac12,j)
}{\Delta z}.
\end{aligned}
\]

Important correction: the \(\sigma_{xz}\) equation uses \(\mu\), not \(\lambda\), for both terms.

---

## 6. Boundary policy for the first implementation

The first implementation should be simple and verifiable. Do not implement PML.

Use one of the following boundary strategies consistently:

### Recommended for first version: zero-derivative boundary by skipping boundary equations

Only assemble equations for points whose required stencil neighbors exist. Boundary variables may be fixed to zero or copied from the previous time layer.

However, this makes the unknown vector irregular unless carefully handled.

### Easier for direct \(Ax=b\) assembly: include ghost-like boundary half grids and use zero rows for invalid derivative stencils

For a small proof-of-concept, keep all staggered arrays in the unknown vector and define boundary derivative rows as zero where a neighbor is missing.

For example, if \(v_x(0,0)=v_x(\frac12,0)\) is treated as an interior half position between \(\sigma_{xx}(0,0)\) and \(\sigma_{xx}(1,0)\), then its \(x\)-derivative term must be:

\[
\frac{\sigma_{xx}(1,0)-\sigma_{xx}(0,0)}{\Delta x}.
\]

The \(z\)-derivative term requires:

\[
\sigma_{xz}\left(\frac12,\frac12\right)-\sigma_{xz}\left(\frac12,-\frac12\right).
\]

If the lower ghost value is not stored, choose a clear boundary rule. For the first simple code, set the missing outside value to zero or set the corresponding derivative term to zero. Do not mix conventions.

The most important requirement: the indexing convention must be explicit and consistent.

---

## 7. Unknown vector layout for \(Ax=b\)

The unknown vector \(x\) contains all new-time-level variables:

\[
x=
\begin{bmatrix}
x_{vx}\\
x_{vz}\\
x_{xx}\\
x_{zz}\\
x_{xz}
\end{bmatrix}
=
\begin{bmatrix}
\mathrm{vec}_r(V_x^{n+1})\\
\mathrm{vec}_r(V_z^{n+1})\\
\mathrm{vec}_r(T_{xx}^{n+1})\\
\mathrm{vec}_r(T_{zz}^{n+1})\\
\mathrm{vec}_r(T_{xz}^{n+1})
\end{bmatrix}.
\]

Here:

\[
T_{xx}=\sigma_{xx},\qquad T_{zz}=\sigma_{zz},\qquad T_{xz}=\sigma_{xz}.
\]

Use row-major vectorization:

\[
\mathrm{vec}_r(U)=[U(0,0),U(0,1),\ldots,U(0,N_x-1),U(1,0),\ldots]^T.
\]

In code:

```c
local_id = ix + iz * nx_local;
```

For general \(N_z\times N_x\) stress grid:

\[
N_{vx}=N_z(N_x+1),
\]

\[
N_{vz}=(N_z+1)N_x,
\]

\[
N_{xx}=N_zN_x,
\]

\[
N_{zz}=N_zN_x,
\]

\[
N_{xz}=(N_z+1)(N_x+1).
\]

Total unknown count:

\[
N_{total}=N_{vx}+N_{vz}+N_{xx}+N_{zz}+N_{xz}.
\]

Block offsets:

```c
off_vx = 0;
off_vz = off_vx + Nz * (Nx + 1);
off_xx = off_vz + (Nz + 1) * Nx;
off_zz = off_xx + Nz * Nx;
off_xz = off_zz + Nz * Nx;
Ntot   = off_xz + (Nz + 1) * (Nx + 1);
```

Index functions:

```c
id_vx(iz, ix) = off_vx + ix + iz * (Nx + 1);      // iz=0..Nz-1, ix=0..Nx
id_vz(iz, ix) = off_vz + ix + iz * Nx;            // iz=0..Nz,   ix=0..Nx-1
id_xx(iz, ix) = off_xx + ix + iz * Nx;            // iz=0..Nz-1, ix=0..Nx-1
id_zz(iz, ix) = off_zz + ix + iz * Nx;            // iz=0..Nz-1, ix=0..Nx-1
id_xz(iz, ix) = off_xz + ix + iz * (Nx + 1);      // iz=0..Nz,   ix=0..Nx
```

---

## 8. Matrix row assembly principle

Each discrete equation contributes one row to \(A\) and one entry to \(b\).

For each unknown at \(n+1\), put its coefficient into matrix \(A\). For each known quantity at \(n\), accumulate it into the right-hand side \(b\).

The row order should match the unknown block order:

1. all \(v_x\) equations;
2. all \(v_z\) equations;
3. all \(\sigma_{xx}\) equations;
4. all \(\sigma_{zz}\) equations;
5. all \(\sigma_{xz}\) equations.

That means the row index of each equation can use the same ordering as its corresponding unknown variable.

---

## 9. Example: row for \(v_x(0,0)=v_x(\frac12,0)\)

Assume a \(4\times4\) stress grid. Then:

\[
V_x\in \mathbb{R}^{4\times5},\quad
V_z\in \mathbb{R}^{5\times4},\quad
T_{xx},T_{zz}\in \mathbb{R}^{4\times4},\quad
T_{xz}\in \mathbb{R}^{5\times5}.
\]

The total unknown count is:

\[
20+20+16+16+25=97.
\]

The block offsets in 0-based indexing are:

```text
x_vx:  0  ... 19
x_vz:  20 ... 39
x_xx:  40 ... 55
x_zz:  56 ... 71
x_xz:  72 ... 96
```

For \(v_x(0,0)=v_x(\frac12,0)\), the CN equation is:

\[
\begin{aligned}
\rho^x_{0,0}(v^x_{0,0})^{n+1}
&-
\frac{\Delta t}{2\Delta x}
\left[(\sigma^{xx}_{0,1})^{n+1}-(\sigma^{xx}_{0,0})^{n+1}\right]
\\
&-
\frac{\Delta t}{2\Delta z}
\left[(\sigma^{xz}_{1,0})^{n+1}-(\sigma^{xz}_{0,0})^{n+1}\right]
=
 b_0.
\end{aligned}
\]

Therefore the nonzero entries of row 0 of \(A\) are:

\[
A_{0,0}=\rho^x_{0,0},
\]

\[
A_{0,40}=+\frac{\Delta t}{2\Delta x},
\]

\[
A_{0,41}=-\frac{\Delta t}{2\Delta x},
\]

\[
A_{0,72}=+\frac{\Delta t}{2\Delta z},
\]

\[
A_{0,77}=-\frac{\Delta t}{2\Delta z}.
\]

The right-hand side is:

\[
\begin{aligned}
b_0
&=
\rho^x_{0,0}(v^x_{0,0})^n
+
\frac{\Delta t}{2\Delta x}
\left[(\sigma^{xx}_{0,1})^{n}-(\sigma^{xx}_{0,0})^{n}\right]
\\
&+
\frac{\Delta t}{2\Delta z}
\left[(\sigma^{xz}_{1,0})^{n}-(\sigma^{xz}_{0,0})^{n}\right].
\end{aligned}
\]

Important: if a term references a value outside the chosen stored grid, apply the chosen boundary rule consistently. Do not silently drop only one direction unless the boundary rule says so.

---

## 10. General row assembly formulas

The following formulas are the most direct way to write the code.

### 10.1 Assemble \(v_x\) equation rows

Loop:

```c
for iz = 0 .. Nz-1
  for ix = 0 .. Nx
```

The equation corresponds to \(v_x(ix,iz)=v_x(i+1/2,j)\).

Always add:

```c
row = id_vx(iz, ix);
A[row][id_vx(iz, ix)] += rho_x[iz][ix];
b[row] += rho_x[iz][ix] * vx_old[iz][ix];
```

For the \(\partial_x \sigma_{xx}\) term, valid when the stress neighbors exist. If \(ix\) maps to the half point between \(Txx(ix-1)\) and \(Txx(ix)\), use that convention. If using the convention that \(v_x(ix,iz)\) lies between \(Txx(ix,iz)\) and \(Txx(ix+1,iz)\), then valid for `ix=0..Nx-2`:

```c
A[row][id_xx(iz, ix)]     += +dt / (2*dx);
A[row][id_xx(iz, ix + 1)] += -dt / (2*dx);

b[row] += dt/(2*dx) * (txx_old[iz][ix+1] - txx_old[iz][ix]);
```

This matches:

\[
-\frac{\Delta t}{2}\frac{\sigma_{xx}^{n+1}(i+1,j)-\sigma_{xx}^{n+1}(i,j)}{\Delta x}
\]

moved to the left-hand side.

For the \(\partial_z \sigma_{xz}\) term, use the stored \(T_{xz}\) convention. If `txz(iz,ix)` and `txz(iz+1,ix)` represent the two vertical neighboring half points, then valid for `iz=0..Nz-1`:

```c
A[row][id_xz(iz, ix)]     += +dt / (2*dz);
A[row][id_xz(iz + 1, ix)] += -dt / (2*dz);

b[row] += dt/(2*dz) * (txz_old[iz+1][ix] - txz_old[iz][ix]);
```

### 10.2 Assemble \(v_z\) equation rows

Loop:

```c
for iz = 0 .. Nz
  for ix = 0 .. Nx-1
```

Add the mass term:

```c
row = id_vz(iz, ix);
A[row][id_vz(iz, ix)] += rho_z[iz][ix];
b[row] += rho_z[iz][ix] * vz_old[iz][ix];
```

For \(\partial_z \sigma_{zz}\), valid when stress neighbors exist. If \(v_z(iz,ix)\) lies between \(Tzz(iz,ix)\) and \(Tzz(iz+1,ix)\), then:

```c
A[row][id_zz(iz, ix)]     += +dt / (2*dz);
A[row][id_zz(iz + 1, ix)] += -dt / (2*dz);

b[row] += dt/(2*dz) * (tzz_old[iz+1][ix] - tzz_old[iz][ix]);
```

For \(\partial_x \sigma_{xz}\):

```c
A[row][id_xz(iz, ix)]     += +dt / (2*dx);
A[row][id_xz(iz, ix + 1)] += -dt / (2*dx);

b[row] += dt/(2*dx) * (txz_old[iz][ix+1] - txz_old[iz][ix]);
```

### 10.3 Assemble \(\sigma_{xx}\) equation rows

Loop:

```c
for iz = 0 .. Nz-1
  for ix = 0 .. Nx-1
```

Add identity:

```c
row = id_xx(iz, ix);
A[row][id_xx(iz, ix)] += 1.0;
b[row] += txx_old[iz][ix];
```

Let:

```c
C = lambda[iz][ix] + 2.0 * mu[iz][ix];
L = lambda[iz][ix];
```

For \(\partial_x v_x\):

```c
A[row][id_vx(iz, ix)]     += +dt * C / (2*dx);
A[row][id_vx(iz, ix + 1)] += -dt * C / (2*dx);

b[row] += dt*C/(2*dx) * (vx_old[iz][ix+1] - vx_old[iz][ix]);
```

For \(\partial_z v_z\):

```c
A[row][id_vz(iz, ix)]     += +dt * L / (2*dz);
A[row][id_vz(iz + 1, ix)] += -dt * L / (2*dz);

b[row] += dt*L/(2*dz) * (vz_old[iz+1][ix] - vz_old[iz][ix]);
```

### 10.4 Assemble \(\sigma_{zz}\) equation rows

Loop:

```c
for iz = 0 .. Nz-1
  for ix = 0 .. Nx-1
```

Add identity:

```c
row = id_zz(iz, ix);
A[row][id_zz(iz, ix)] += 1.0;
b[row] += tzz_old[iz][ix];
```

Let:

```c
C = lambda[iz][ix] + 2.0 * mu[iz][ix];
L = lambda[iz][ix];
```

For \(\partial_x v_x\):

```c
A[row][id_vx(iz, ix)]     += +dt * L / (2*dx);
A[row][id_vx(iz, ix + 1)] += -dt * L / (2*dx);

b[row] += dt*L/(2*dx) * (vx_old[iz][ix+1] - vx_old[iz][ix]);
```

For \(\partial_z v_z\):

```c
A[row][id_vz(iz, ix)]     += +dt * C / (2*dz);
A[row][id_vz(iz + 1, ix)] += -dt * C / (2*dz);

b[row] += dt*C/(2*dz) * (vz_old[iz+1][ix] - vz_old[iz][ix]);
```

### 10.5 Assemble \(\sigma_{xz}\) equation rows

Loop:

```c
for iz = 0 .. Nz
  for ix = 0 .. Nx
```

Add identity:

```c
row = id_xz(iz, ix);
A[row][id_xz(iz, ix)] += 1.0;
b[row] += txz_old[iz][ix];
```

Let:

```c
M = mu_xz[iz][ix];
```

For \(\partial_x v_z\), valid when `ix < Nx` and the neighbor exists:

```c
A[row][id_vz(iz, ix)]     += +dt * M / (2*dx);
A[row][id_vz(iz, ix + 1)] += -dt * M / (2*dx);

b[row] += dt*M/(2*dx) * (vz_old[iz][ix+1] - vz_old[iz][ix]);
```

For \(\partial_z v_x\), valid when `iz < Nz` and the neighbor exists:

```c
A[row][id_vx(iz, ix)]     += +dt * M / (2*dz);
A[row][id_vx(iz + 1, ix)] += -dt * M / (2*dz);

b[row] += dt*M/(2*dz) * (vx_old[iz+1][ix] - vx_old[iz][ix]);
```

Boundary rows require a clear rule. For a first version, if a neighbor does not exist, skip that derivative term or impose a fixed boundary value. But document the choice.

---

## 11. Block matrix form

The full system has the block form:

\[
\begin{bmatrix}
R_x
&0
&-\frac{\Delta t}{2}K_{x,xx}^{c\rightarrow h}
&0
&-\frac{\Delta t}{2}K_{z,xz}^{h\rightarrow c}
\\
0
&R_z
&0
&-\frac{\Delta t}{2}K_{z,zz}^{c\rightarrow h}
&-\frac{\Delta t}{2}K_{x,xz}^{h\rightarrow c}
\\
-\frac{\Delta t}{2}C K_{x,vx}^{h\rightarrow c}
&-\frac{\Delta t}{2}L K_{z,vz}^{h\rightarrow c}
&I
&0
&0
\\
-\frac{\Delta t}{2}L K_{x,vx}^{h\rightarrow c}
&-\frac{\Delta t}{2}C K_{z,vz}^{h\rightarrow c}
&0
&I
&0
\\
-\frac{\Delta t}{2}M_{xz}K_{z,vx}^{c\rightarrow h}
&-\frac{\Delta t}{2}M_{xz}K_{x,vz}^{c\rightarrow h}
&0
&0
&I
\end{bmatrix}
\begin{bmatrix}
x_{vx}\\
x_{vz}\\
x_{xx}\\
x_{zz}\\
x_{xz}
\end{bmatrix}
=
\begin{bmatrix}
b_{vx}\\
b_{vz}\\
b_{xx}\\
b_{zz}\\
b_{xz}
\end{bmatrix}.
\]

Here:

- \(R_x=\mathrm{diag}(\mathrm{vec}_r(\rho_x))\);
- \(R_z=\mathrm{diag}(\mathrm{vec}_r(\rho_z))\);
- \(L=\mathrm{diag}(\mathrm{vec}_r(\lambda))\);
- \(C=\mathrm{diag}(\mathrm{vec}_r(\lambda+2\mu))\);
- \(M_{xz}=\mathrm{diag}(\mathrm{vec}_r(\mu_{xz}))\).

Do not confuse this with an explicit matrix-multiplication stencil. This code must assemble \(A\) and solve \(Ax=b\) directly.

---

## 12. Source term

Use an isotropic stress point source. At source grid point \((i_s,j_s)\), add the source to the right-hand side of the \(T_{xx}\) and \(T_{zz}\) equations.

For example, at each time step:

```c
sou = ricker(t);
b[id_xx(js, is)] += dt * sou;
b[id_zz(js, is)] += dt * sou;
```

or equivalently add the source to the old stress fields before assembling \(b\). Prefer adding directly to \(b\) so the linear system remains clear.

A Ricker wavelet can be:

\[
sou(t)=\left(1-2\pi^2 f_0^2(t-t_0)^2\right)
\exp\left[-\pi^2 f_0^2(t-t_0)^2\right].
\]

If using the negative sign convention from an older code, document it clearly. The sign only changes source polarity.

---

## 13. Material interpolation

For the first uniform-medium version, all parameters can be constants:

\[
\rho_x=\rho_z=2700.
\]

\[
\mu=2700\times 3000^2.
\]

\[
\lambda=2700\times(6000^2-2\times3000^2).
\]

\[
\lambda+2\mu=2700\times6000^2.
\]

\[
\mu_{xz}=\mu.
\]

For later heterogeneous models:

- \(\rho_x\) should be interpolated to the \(v_x\) location;
- \(\rho_z\) should be interpolated to the \(v_z\) location;
- \(\mu_{xz}\) should be interpolated to the \(\sigma_{xz}\) location, preferably by harmonic averaging near material interfaces.

---

## 14. Suggested code structure

A clear serial implementation should contain:

```text
main.c / main.cpp
  - define grid and time parameters
  - allocate arrays
  - initialize model
  - initialize wavefields
  - time loop
      - assemble_A_and_b(...)
      - add_source_to_b(...)
      - solve_linear_system(A, b, x)
      - unpack x into vx, vz, txx, tzz, txz
      - output snapshots or traces
```

Recommended functions:

```c
int id_vx(int iz, int ix);
int id_vz(int iz, int ix);
int id_xx(int iz, int ix);
int id_zz(int iz, int ix);
int id_xz(int iz, int ix);

void initialize_model(...);
void initialize_fields(...);
void assemble_A_and_b(...);
void add_source_to_b(...);
void solve_dense_or_sparse(...);
void unpack_solution(...);
void write_snapshot(...);
```

For the first implementation, a dense matrix is acceptable for very small grids such as \(4\times4\), \(10\times10\), or \(20\times20\). For larger grids, use sparse storage such as CSR.

---

## 15. Solver recommendation

For a simple serial prototype:

1. Small grid: use dense Gaussian elimination or LU;
2. Medium grid: use sparse CSR + BiCGSTAB or GMRES;
3. The matrix is generally nonsymmetric because of variable coefficients, staggered layout, and boundary handling;
4. Do not use conjugate gradient unless you have proven the matrix is symmetric positive definite.

For initial testing, choose a very small grid and dense LU to avoid solver complications.

---

## 16. Verification checklist

The AI-generated code should be checked using the following tests:

### 16.1 Matrix assembly sanity

For a \(4\times4\) stress grid:

\[
N_{total}=97.
\]

Check that:

```text
A is 97 x 97
x is length 97
b is length 97
```

### 16.2 Row check for \(v_x(0,0)\)

The row for \(v_x(0,0)=v_x(\frac12,0)\) should include coefficients for:

- \(v_x(0,0)\);
- \(\sigma_{xx}(0,0)\);
- \(\sigma_{xx}(0,1)\), depending on row-major notation this means same z row, next x column;
- \(\sigma_{xz}(0,0)\);
- \(\sigma_{xz}(1,0)\).

The signs should be:

\[
+\rho_x,
\quad
+\frac{\Delta t}{2\Delta x},
\quad
-\frac{\Delta t}{2\Delta x},
\quad
+\frac{\Delta t}{2\Delta z},
\quad
-\frac{\Delta t}{2\Delta z}.
\]

### 16.3 CN residual

After solving \(Ax=b\), compute:

\[
r=A x-b.
\]

The relative residual should be small:

\[
\frac{\|Ax-b\|_2}{\|b\|_2}\ll 1.
\]

### 16.4 Physical sanity

For homogeneous medium:

- P-wave speed should be approximately \(6000\);
- S-wave speed should be approximately \(3000\);
- wavefronts should be physically reasonable;
- if no absorbing boundary is used, reflections from boundaries are expected.

### 16.5 Difference from explicit method

The CN implicit result should not be expected to match explicit leapfrog exactly. They solve the same continuous PDE but use different time discretizations. As \(\Delta t,\Delta x,\Delta z\) decrease, both should converge toward the same physical solution.

---

## 17. Important implementation warnings

1. Do not implement Picard iteration as a substitute for \(Ax=b\) unless explicitly requested. This project requires direct matrix assembly and solution.
2. Do not use \(\lambda\) in the \(\sigma_{xz}\) equation. Use \(\mu\) for both shear terms.
3. Keep the variable block order fixed:

```text
Vx, Vz, Txx, Tzz, Txz
```

4. Keep row-major vectorization fixed.
5. Carefully handle boundary derivative terms. The first version may skip invalid terms or impose zero derivative, but the choice must be documented.
6. Do not add PML in the first implementation.
7. The source should be added to both \(T_{xx}\) and \(T_{zz}\) equations.
8. Use double precision for the first serial verification if possible. Single precision can be used later.

---

## 18. Minimal expected output

The code should output at least:

1. source time function;
2. receiver trace of \(T_{xx}\) or pressure-like quantity \((T_{xx}+T_{zz})/2\);
3. wavefield snapshots at selected time steps;
4. matrix residual \(\|Ax-b\|/\|b\|\) for each time step or selected time steps.

---

## 19. Summary for AI code generation

Write a serial 2D elastic-wave solver using the velocity–stress formulation. Use staggered-grid finite differences in space and Crank–Nicolson in time. At every time step, assemble the coupled linear system \(Ax=b\) for the five new-time-level unknown fields:

\[
V_x^{n+1},\quad V_z^{n+1},\quad T_{xx}^{n+1},\quad T_{zz}^{n+1},\quad T_{xz}^{n+1}.
\]

Use the exact block order:

\[
x=[\mathrm{vec}_r(V_x),\mathrm{vec}_r(V_z),\mathrm{vec}_r(T_{xx}),\mathrm{vec}_r(T_{zz}),\mathrm{vec}_r(T_{xz})]^T.
\]

Use the CN discrete equations listed in Section 5. Assemble all new-time-layer terms into \(A\), all old-time-layer terms into \(b\). Solve \(Ax=b\) directly for small grids. Do not use explicit leapfrog and do not use Picard iteration. The first version should be homogeneous, serial, and without PML.
