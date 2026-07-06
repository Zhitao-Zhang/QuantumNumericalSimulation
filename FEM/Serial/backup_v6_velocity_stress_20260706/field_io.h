/* =============================================================================
 *  field_io.h
 *  ---------------------------------------------------------------------------
 *  ASCII 二维矩阵输出。签名与 Quantum/FDM/ImplicitV2/main.cpp 中的
 *  wfile2d(...) 一致，便于两个程序的输出直接对照。
 *
 *  文件格式（每行 ncols 个 %e，共 nrows 行）：
 *      v[0,0] v[0,1] ... v[0,ncols-1]
 *      v[1,0] v[1,1] ...
 *      ...
 *
 *  与 FDM 的关键差异：
 *      FDM 的 vx/vz/σxz 是交错网格（尺寸 Nz×(Nx+1) 等），FEM 里 5 个场都在
 *      同一 Nz×Nx 节点网格上。σxx / σzz 尺寸完全相同（Nz×Nx），可直接逐点对照。
 * ===========================================================================*/
#ifndef FEM_FIELD_IO_H
#define FEM_FIELD_IO_H

#include <vector>

namespace femio {

// float 版本 (与 FDM 完全同名同签名)
void wfile2d(const char* fn, const float* data, int nrows, int ncols);

// double vector 版本（内部转 float 写盘，保持 %e 格式一致）
void wfile2d_d(const char* fn, const std::vector<double> &v, int nrows, int ncols);

// 从块向量 q 中抽出 5 个场并写盘 (每个都是 Nz x Nx 的节点网格)
void write_snapshot(int it,
                    const std::vector<double> &q,
                    int Nx, int Nz,
                    int off_vx, int off_vz, int off_xx, int off_zz, int off_xz);

} // namespace femio

#endif // FEM_FIELD_IO_H
