/* =============================================================================
 *  field_io.cpp
 * ===========================================================================*/
#include "field_io.h"
#include "config.h"

#include <cstdio>
#include <cstddef>

namespace femio {

void wfile2d(const char* fn, const float* data, int nrows, int ncols)
{
    FILE* fp = std::fopen(fn, "wt");
    if (!fp) return;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++)
            std::fprintf(fp, "%e ", data[(size_t)i * (size_t)ncols + (size_t)j]);
        std::fprintf(fp, "\n");
    }
    std::fclose(fp);
}

void wfile2d_d(const char* fn, const std::vector<double> &v, int nrows, int ncols)
{
    std::vector<float> tmp((size_t)nrows * (size_t)ncols);
    for (size_t k = 0; k < tmp.size(); k++) tmp[k] = (float)v[k];
    wfile2d(fn, tmp.data(), nrows, ncols);
}

// 输入场按细网格 Nnodes(=MX*MZ) 存储；快照只导出【顶点子网格】(Nz×Nx)，
// 用 vert_node_id 抽取，保持文件格式与 P1/FDM 完全一致（画图脚本不用改）。
void write_snapshot_fields(int it, int Nx, int Nz,
                           const std::vector<double> &vx,
                           const std::vector<double> &vz,
                           const std::vector<double> &sxx,
                           const std::vector<double> &szz,
                           const std::vector<double> &sxz)
{
    char fn[128];
    const size_t Nout = (size_t)Nx * (size_t)Nz;
    std::vector<float> tmp(Nout);
    auto dump = [&](const std::vector<double> &fld, const char *tag){
        for (int iz = 0; iz < Nz; iz++)
            for (int ix = 0; ix < Nx; ix++)
                tmp[(size_t)iz * (size_t)Nx + (size_t)ix]
                    = (float)fld[femcfg::vert_node_id(ix, iz)];
        std::snprintf(fn, sizeof(fn), "snap_%s_it%04d.dat", tag, it);
        wfile2d(fn, tmp.data(), Nz, Nx);
    };
    dump(vx,  "vx" );
    dump(vz,  "vz" );
    dump(sxx, "txx");
    dump(szz, "tzz");
    dump(sxz, "txz");
}

} // namespace femio
