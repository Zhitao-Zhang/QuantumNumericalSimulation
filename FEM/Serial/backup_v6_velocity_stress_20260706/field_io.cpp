/* =============================================================================
 *  field_io.cpp
 * ===========================================================================*/
#include "field_io.h"

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

void write_snapshot(int it,
                    const std::vector<double> &q,
                    int Nx, int Nz,
                    int off_vx, int off_vz, int off_xx, int off_zz, int off_xz)
{
    char fn[128];
    const size_t Nnodes = (size_t)Nx * (size_t)Nz;

    std::vector<float> tmp(Nnodes);
    auto dump = [&](int base, const char *tag){
        for (size_t k = 0; k < Nnodes; k++) tmp[k] = (float)q[base + k];
        std::snprintf(fn, sizeof(fn), "snap_%s_it%04d.dat", tag, it);
        wfile2d(fn, tmp.data(), Nz, Nx);
    };
    dump(off_vx, "vx" );
    dump(off_vz, "vz" );
    dump(off_xx, "txx");
    dump(off_zz, "tzz");
    dump(off_xz, "txz");
}

} // namespace femio
