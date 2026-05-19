/*
 * 2D 弹性波 SGFD + SPML，时间方向：Crank–Nicolson 的 Picard 固定点近似。
 *
 * 记一步显式子步（Source 之后仅 FD_V + FD_T）为映射 R(u)。CN 离散为
 *   u^{n+1} = u^s + (1/2)(R(u^s)-u^s) + (1/2)(R(u^{n+1})-u^{n+1})
 * 用 Picard 迭代逼近 u^{n+1}（初值 u^{(0)}=R(u^s)）。
 *
 * 网格与模型与 SGFDM2D.cu 首段 main 一致：NY=128, NX=128, NP=32。
 *
 * 为何低频一致、高频隐式更“光滑”：
 * Picard 逼近 CN；若迭代已收敛，再增加次数不改变结果。与显式差异还来自时间离散。
 *
 * 若改 N_PICARD 结果不变：多半是 Picard 在前几步已在 float 下收敛，多出来的
 * 循环不改变状态。另：若曾开启「按主频加迭代」，不同 N_PICARD 可能被算到同一
 * n_picard（都顶到上限），看起来像改宏无效——现默认已关闭该加步，见下。
 *
 * 调试：export QN_PICARD_TRACE=1；export QN_N_PICARD=0 为跳过 Picard（与每步仅
 * R(u^s) 的显式子步对齐，波形应与多轮 Picard 明显不同）。
 * 若需按 F0 自动加 Picard 次数：编译加 -DPICARD_FREQ_BOOST。
 */
 #include "cuda_runtime.h"
 #include "device_launch_parameters.h"
 
 #include <stdio.h>
 #include <stdlib.h>
 #include <string.h>
 #include <math.h>
 #include <time.h>
 
 #define PI 3.1415926535f
 #define NY 328
 #define NX 128
 #define NP 32
 #define BLOCK_SIZE_X 8
 #define BLOCK_SIZE_Y 16
 
 /* Picard 迭代次数：默认即宏值（不再默认按 F0 加大，避免与宏不一致、顶到上限后“改宏无效”） */
 #ifndef N_PICARD
 #define N_PICARD 3
 #endif
 #ifndef N_PICARD_MAX
 #define N_PICARD_MAX 256
 #endif
 
 #define NX_EXT (NX + 2 * NP)
 #define NY_EXT (NY + 2 * NP)
 #define IDX2(ix, iy, nxe) ((ix) + (iy) * (nxe))
 
 int divUp(int a, int b) { return (a - 1) / b + 1; }
 
 void wfile2d(const char* fn, const float* data, int nrows, int ncols);
 void create_model_2d(float* vp, float* vs, float* rho, float* lamda, float* lamda2u, float* mu, int ny, int nx);
 void extmodel2d_pad(const float* src, float* dst, int ny, int nx, int np);
 
 __global__ void Source2D(float* txx, float* tyy, float I_sou, int sn)
 {
     int ix = threadIdx.x + blockIdx.x * blockDim.x;
     int iy = threadIdx.y + blockIdx.y * blockDim.y;
     int o = IDX2(ix, iy, NX_EXT);
     if (o == sn)
     {
         txx[o] += I_sou;
         tyy[o] += I_sou;
     }
 }
 
 template <int BLOCK_X, int BLOCK_Y>
 __global__ void FD_V_2D(
     float* vux, float* vuy,
     float* rho_tempx, float* rho_tempy,
     float* txx, float* tyy, float* txy,
     float* pmlxSxx, float* pmlySxy, float* pmlxSxy, float* pmlySyy,
     float* SXxx, float* SXxy, float* SYxy, float* SYyy,
     float* e_dxi, float* dxi, float* e_dxi2, float* dxi2, float* e_dyj, float* dyj, float* e_dyj2, float* dyj2,
     float DT, float H)
 {
     const int nxe = NX_EXT;
     const int nye = NY_EXT;
     int ix = blockIdx.x * blockDim.x + threadIdx.x;
     int iy = blockIdx.y * blockDim.y + threadIdx.y;
     if (ix <= 0 || ix >= nxe - 1 || iy <= 0 || iy >= nye - 1)
         return;
     int o = IDX2(ix, iy, nxe);
     int offset_r = IDX2(ix + 1, iy, nxe);
     int offset_h = IDX2(ix, iy - 1, nxe);
     int offset_l = IDX2(ix - 1, iy, nxe);
     int offset_q = IDX2(ix, iy + 1, nxe);
 
     float x1 = (txx[offset_r] - txx[o]) / H;
     float y1 = (tyy[offset_q] - tyy[o]) / H;
     float x2 = (txy[o] - txy[offset_h]) / H;
     float y2 = (txy[o] - txy[offset_l]) / H;
 
     pmlxSxx[o] = pmlxSxx[o] * e_dxi2[o] + (-DT * dxi2[o] * 0.5f) * (e_dxi2[o] * SXxx[o] + x1);
     pmlySxy[o] = pmlySxy[o] * e_dyj[o] + (-DT * dyj[o] * 0.5f) * (e_dyj[o] * SXxy[o] + x2);
     SXxx[o] = x1; SXxy[o] = x2;
 
     pmlxSxy[o] = pmlxSxy[o] * e_dxi[o] + (-DT * dxi[o] * 0.5f) * (e_dxi[o] * SYxy[o] + y2);
     pmlySyy[o] = pmlySyy[o] * e_dyj2[o] + (-DT * dyj2[o] * 0.5f) * (e_dyj2[o] * SYyy[o] + y1);
     SYxy[o] = y2; SYyy[o] = y1;
 
     y2 = y2 + pmlxSxy[o];
     y1 = y1 + pmlySyy[o];
     x1 = x1 + pmlxSxx[o];
     x2 = x2 + pmlySxy[o];
 
     vux[o] += DT * rho_tempx[o] * (x1 + x2);
     vuy[o] += DT * rho_tempy[o] * (y1 + y2);
 }
 
 template <int BLOCK_X, int BLOCK_Y>
 __global__ void FD_T_2D(
     float* vux, float* vuy,
     float* txx, float* tyy, float* txy,
     float* lamda2u_ext, float* lamda_ext,
     float* pmlxVux, float* muxy,
     float* pmlyVuy, float* pmlxVuy, float* pmlyVux,
     float* Vuxx, float* Vuxy, float* Vuyx, float* Vuyy,
     float* e_dxi, float* dxi, float* e_dxi2, float* dxi2, float* e_dyj2, float* dyj2, float* e_dyj, float* dyj,
     float DT, float H)
 {
     const int nxe = NX_EXT;
     const int nye = NY_EXT;
     int ix = threadIdx.x + blockIdx.x * blockDim.x;
     int iy = threadIdx.y + blockIdx.y * blockDim.y;
     if (ix <= 0 || ix >= nxe - 1 || iy <= 0 || iy >= nye - 1)
         return;
     int o = IDX2(ix, iy, nxe);
 
     float uxx = (vux[o] - vux[IDX2(ix - 1, iy, nxe)]) / H;
     float uxy = (vux[IDX2(ix, iy + 1, nxe)] - vux[o]) / H;
     float uyy = (vuy[o] - vuy[IDX2(ix, iy - 1, nxe)]) / H;
     float uyx = (vuy[IDX2(ix + 1, iy, nxe)] - vuy[o]) / H;
 
     pmlxVux[o] = pmlxVux[o] * e_dxi[o] + (-DT * dxi[o] * 0.5f) * (e_dxi[o] * Vuxx[o] + uxx);
     pmlxVuy[o] = pmlxVuy[o] * e_dyj2[o] + (-DT * dyj2[o] * 0.5f) * (e_dyj2[o] * Vuxy[o] + uxy);
     Vuxx[o] = uxx; Vuxy[o] = uxy;
 
     pmlyVuy[o] = pmlyVuy[o] * e_dyj[o] + (-DT * dyj[o] * 0.5f) * (e_dyj[o] * Vuyy[o] + uyy);
     pmlyVux[o] = pmlyVux[o] * e_dxi2[o] + (-DT * dxi2[o] * 0.5f) * (e_dxi2[o] * Vuyx[o] + uyx);
     Vuyy[o] = uyy; Vuyx[o] = uyx;
 
     uxx += pmlxVux[o];
     uyy += pmlyVuy[o];
     uxy += pmlxVuy[o];
     uyx += pmlyVux[o];
 
     txx[o] += DT * (lamda_ext[o] * (uyy)+lamda2u_ext[o] * uxx);
     tyy[o] += DT * (lamda_ext[o] * (uxx)+lamda2u_ext[o] * uyy);
     txy[o] += muxy[o] * DT * (uxy + uyx);
 }
 
 __global__ void CrankBlend1(const float* bak, const float* rs, const float* cur, const float* wrk, float* out, int n2)
 {
     int o = blockIdx.x * blockDim.x + threadIdx.x;
     if (o >= n2) return;
     out[o] = bak[o] + 0.5f * ((rs[o] - bak[o]) + (wrk[o] - cur[o]));
 }
 
 #define N_STATE_FIELDS 21
 
 static void copy_state_fields(float** dst, float** src, size_t mem2)
 {
     for (int i = 0; i < N_STATE_FIELDS; i++)
         cudaMemcpy(dst[i], src[i], mem2, cudaMemcpyDeviceToDevice);
 }
 
 static void blend_all_fields(
     float** bak, float** rs, float** cur, float** wrk, float** out, int n2, int threads)
 {
     int blocks = (n2 + threads - 1) / threads;
     for (int i = 0; i < N_STATE_FIELDS; i++)
         CrankBlend1 <<<blocks, threads>>> (bak[i], rs[i], cur[i], wrk[i], out[i], n2);
     cudaDeviceSynchronize();
 }

 /* 默认 n_picard == N_PICARD。若编译加 -DPICARD_FREQ_BOOST，则再按 F0 加若干步（旧行为） */
 static int picard_iters_for_source(float F0)
 {
     int n = N_PICARD;
 #ifdef PICARD_FREQ_BOOST
     int boost = (int)(0.0025f * F0 + 0.5f);
     if (boost > 36)
         boost = 36;
     n += boost;
 #endif
     if (n < 0)
         n = 0;
     if (n > N_PICARD_MAX)
         n = N_PICARD_MAX;
     return n;
 }

 /* QN_N_PICARD：未设置则 -1；>=0 则覆盖（0=完全跳过 Picard，仅 R(u^s)，便于对照） */
 static int picard_iters_from_env_override(void)
 {
     const char* e = getenv("QN_N_PICARD");
     if (!e || !e[0])
         return -1;
     int v = atoi(e);
     if (v < 0)
         return -1;
     if (v > N_PICARD_MAX)
         v = N_PICARD_MAX;
     return v;
 }
 
 template <int BLOCK_X, int BLOCK_Y>
 static void launch_VT(
     float* vux, float* vuy, float* rtx, float* rty, float* txx, float* tyy, float* txy,
     float* pmlxSxx, float* pmlySxy, float* pmlxSxy, float* pmlySyy,
     float* SXxx, float* SXxy, float* SYxy, float* SYyy,
     float* pmlxVux, float* muxy, float* pmlyVuy, float* pmlxVuy, float* pmlyVux,
     float* Vuxx, float* Vuxy, float* Vuyx, float* Vuyy,
     float* d_lam2, float* d_lam,
     float* e_dxi, float* dxi, float* e_dxi2, float* dxi2, float* e_dyj, float* dyj, float* e_dyj2, float* dyj2,
     float DT, float H, dim3 Grid1, dim3 Block1)
 {
     FD_V_2D<BLOCK_X, BLOCK_Y> <<<Grid1, Block1>>> (
         vux, vuy, rtx, rty, txx, tyy, txy,
         pmlxSxx, pmlySxy, pmlxSxy, pmlySyy,
         SXxx, SXxy, SYxy, SYyy,
         e_dxi, dxi, e_dxi2, dxi2, e_dyj, dyj, e_dyj2, dyj2, DT, H);
     FD_T_2D<BLOCK_X, BLOCK_Y> <<<Grid1, Block1>>> (
         vux, vuy, txx, tyy, txy, d_lam2, d_lam,
         pmlxVux, muxy, pmlyVuy, pmlxVuy, pmlyVux,
         Vuxx, Vuxy, Vuyx, Vuyy,
         e_dxi, dxi, e_dxi2, dxi2, e_dyj2, dyj2, e_dyj, dyj, DT, H);
 }
 
 static void fill_state_ptrs(
     float** arr,
     float* vux, float* vuy, float* txx, float* tyy, float* txy,
     float* pmlxSxx, float* pmlySxy, float* pmlxSxy, float* pmlySyy,
     float* SXxx, float* SXxy, float* SYxy, float* SYyy,
     float* pmlxVux, float* pmlyVuy, float* pmlxVuy, float* pmlyVux,
     float* Vuxx, float* Vuxy, float* Vuyx, float* Vuyy)
 {
     int k = 0;
     arr[k++] = vux; arr[k++] = vuy; arr[k++] = txx; arr[k++] = tyy; arr[k++] = txy;
     arr[k++] = pmlxSxx; arr[k++] = pmlySxy; arr[k++] = pmlxSxy; arr[k++] = pmlySyy;
     arr[k++] = SXxx; arr[k++] = SXxy; arr[k++] = SYxy; arr[k++] = SYyy;
     arr[k++] = pmlxVux; arr[k++] = pmlyVuy; arr[k++] = pmlxVuy; arr[k++] = pmlyVux;
     arr[k++] = Vuxx; arr[k++] = Vuxy; arr[k++] = Vuyx; arr[k++] = Vuyy;
 }
 
 //--------------------------------------------------------------------------
 int main(void)
 {
     const int NX_ext = NX_EXT;
     const int NY_ext = NY_EXT;
     const int n2 = NX_ext * NY_ext;
     const int BLOCK_X = BLOCK_SIZE_X;
     const int BLOCK_Y = BLOCK_SIZE_Y;
     const int blend_threads = 256;
 
     int sx = NX_ext / 2;
     int sy = 40;
     int sn = IDX2(sx, sy, NX_ext);
 
     int NT = 3500;
     float H = 0.01f;
     float RC = 1.0e-6f;
     float DT = 1.0e-6f;
     float DP = (float)NP * H;
     float F0 = 2e3f;
     float T0 = 1.2f / F0;
     float Vpmax = 4270.f;
     float Vsmin = 1500.f;
 
     dim3 Block1((unsigned)BLOCK_SIZE_X, (unsigned)BLOCK_SIZE_Y, 1u);
     dim3 Grid1((unsigned)divUp(NX_ext, BLOCK_SIZE_X), (unsigned)divUp(NY_ext, BLOCK_SIZE_Y), 1u);
     dim3 BlockSrc((unsigned)BLOCK_SIZE_X, (unsigned)BLOCK_SIZE_Y, 1u);
     dim3 GridSrc((unsigned)divUp(NX_ext, BLOCK_SIZE_X), (unsigned)divUp(NY_ext, BLOCK_SIZE_Y), 1u);
 
     size_t mem2 = (size_t)n2 * sizeof(float);
 
     float* vp = (float*)calloc((size_t)NX * NY, sizeof(float));
     float* vs = (float*)calloc((size_t)NX * NY, sizeof(float));
     float* rho = (float*)calloc((size_t)NX * NY, sizeof(float));
     float* lamda = (float*)calloc((size_t)NX * NY, sizeof(float));
     float* lamda2u = (float*)calloc((size_t)NX * NY, sizeof(float));
     float* mu = (float*)calloc((size_t)NX * NY, sizeof(float));
     create_model_2d(vp, vs, rho, lamda, lamda2u, mu, NY, NX);
 
     float* lam_e = (float*)calloc((size_t)n2, sizeof(float));
     float* lam2u_e = (float*)calloc((size_t)n2, sizeof(float));
     float* mu_e = (float*)calloc((size_t)n2, sizeof(float));
     float* rho_e = (float*)calloc((size_t)n2, sizeof(float));
     extmodel2d_pad(lamda, lam_e, NY, NX, NP);
     extmodel2d_pad(lamda2u, lam2u_e, NY, NX, NP);
     extmodel2d_pad(mu, mu_e, NY, NX, NP);
     extmodel2d_pad(rho, rho_e, NY, NX, NP);
     free(vp); free(vs); free(rho); free(lamda); free(lamda2u); free(mu);
 
     float* rho_tx = (float*)calloc((size_t)n2, sizeof(float));
     float* rho_ty = (float*)calloc((size_t)n2, sizeof(float));
     float* muxy = (float*)calloc((size_t)n2, sizeof(float));
 
     int ix, iy;
     for (iy = 1; iy < NY_ext - 1; iy++)
     {
         for (ix = 1; ix < NX_ext - 1; ix++)
         {
             int o = IDX2(ix, iy, NX_ext);
             rho_tx[o] = 2.f / (rho_e[o] + rho_e[IDX2(ix + 1, iy, NX_ext)]);
             rho_ty[o] = 2.f / (rho_e[o] + rho_e[IDX2(ix, iy + 1, NX_ext)]);
         }
     }
     for (iy = 1; iy < NY_ext - 1; iy++)
     {
         for (ix = 1; ix < NX_ext - 1; ix++)
         {
             int o = IDX2(ix, iy, NX_ext);
             float m0 = mu_e[o], m1 = mu_e[IDX2(ix + 1, iy, NX_ext)], m2 = mu_e[IDX2(ix, iy + 1, NX_EXT)], m3 = mu_e[IDX2(ix + 1, iy + 1, NX_EXT)];
             if (m0 == 0.f || m1 == 0.f || m2 == 0.f || m3 == 0.f)
                 muxy[o] = 0.f;
             else
                 muxy[o] = 1.f / (0.25f * (1.f / m0 + 1.f / m1 + 1.f / m2 + 1.f / m3));
         }
     }
     free(rho_e);
 
     float* dxi = (float*)calloc((size_t)n2, sizeof(float));
     float* dxi2 = (float*)calloc((size_t)n2, sizeof(float));
     float* dyj = (float*)calloc((size_t)n2, sizeof(float));
     float* dyj2 = (float*)calloc((size_t)n2, sizeof(float));
     float* e_dxi = (float*)calloc((size_t)n2, sizeof(float));
     float* e_dxi2 = (float*)calloc((size_t)n2, sizeof(float));
     float* e_dyj = (float*)calloc((size_t)n2, sizeof(float));
     float* e_dyj2 = (float*)calloc((size_t)n2, sizeof(float));
 
     float xoleft = DP, xoright = (NX_ext - 1) * H - DP;
     float yoleft = DP, yoright = (NY_ext - 1) * H - DP;
     float d0, v0 = 4000.f;
     for (iy = 0; iy < NY_ext; iy++)
     {
         for (ix = 0; ix < NX_ext; ix++)
         {
             float x = ix * H, y = iy * H;
             int o = IDX2(ix, iy, NX_EXT);
             e_dxi[o] = e_dxi2[o] = e_dyj[o] = e_dyj2[o] = 1.f;
 
             if (x < xoleft)
             {
                 d0 = 3.f * v0 * logf(1.f / RC) / (2.f * DP);
                 dxi[o] = d0 * powf((xoleft - x) / DP, 2.f);
                 dxi2[o] = d0 * powf((xoleft - x - 0.5f * H) / DP, 2.f);
                 e_dxi[o] = expf(-dxi[o] * DT);
                 e_dxi2[o] = expf(-dxi2[o] * DT);
             }
             if (x > xoright)
             {
                 d0 = 3.f * v0 * logf(1.f / RC) / (2.f * DP);
                 dxi[o] = d0 * powf((x - xoright) / DP, 2.f);
                 dxi2[o] = d0 * powf((x + 0.5f * H - xoright) / DP, 2.f);
                 e_dxi[o] = expf(-dxi[o] * DT);
                 e_dxi2[o] = expf(-dxi2[o] * DT);
             }
             if (y < yoleft)
             {
                 d0 = 3.f * v0 * logf(1.f / RC) / (2.f * DP);
                 dyj[o] = d0 * powf((yoleft - y) / DP, 2.f);
                 dyj2[o] = d0 * powf((yoleft - y - 0.5f * H) / DP, 2.f);
                 e_dyj[o] = expf(-dyj[o] * DT);
                 e_dyj2[o] = expf(-dyj2[o] * DT);
             }
             if (y > yoright)
             {
                 d0 = 3.f * v0 * logf(1.f / RC) / (2.f * DP);
                 dyj[o] = d0 * powf((y - yoright) / DP, 2.f);
                 dyj2[o] = d0 * powf((y + 0.5f * H - yoright) / DP, 2.f);
                 e_dyj[o] = expf(-dyj[o] * DT);
                 e_dyj2[o] = expf(-dyj2[o] * DT);
             }
         }
     }
 
 #define HALLOC(name) float* name = (float*)calloc((size_t)n2, sizeof(float));
     HALLOC(h_vux) HALLOC(h_vuy)
         HALLOC(h_txx) HALLOC(h_tyy) HALLOC(h_txy)
         HALLOC(h_pmlxSxx) HALLOC(h_pmlySxy) HALLOC(h_pmlxSxy) HALLOC(h_pmlySyy)
         HALLOC(h_SXxx) HALLOC(h_SXxy) HALLOC(h_SYxy) HALLOC(h_SYyy)
         HALLOC(h_pmlxVux) HALLOC(h_pmlyVuy) HALLOC(h_pmlxVuy) HALLOC(h_pmlyVux)
         HALLOC(h_Vuxx) HALLOC(h_Vuxy) HALLOC(h_Vuyx) HALLOC(h_Vuyy)
 #undef HALLOC
 
     float* sis = (float*)calloc((size_t)NY_ext * (size_t)NT, sizeof(float));
 
     float* d_lam, * d_lam2, * d_rtx, * d_rty, * d_muxy;
     float* d_vux, * d_vuy, * d_txx, * d_tyy, * d_txy;
     float* d_SXxx, * d_SXxy, * d_SYxy, * d_SYyy;
     float* d_Vuxx, * d_Vuxy, * d_Vuyx, * d_Vuyy;
     float* d_pmlxSxx, * d_pmlySxy, * d_pmlxSxy, * d_pmlySyy;
     float* d_pmlxVux, * d_pmlyVuy, * d_pmlxVuy, * d_pmlyVux;
     float* d_dxi, * d_dxi2, * d_dyj, * d_dyj2, * d_e_dxi, * d_e_dxi2, * d_e_dyj, * d_e_dyj2;
 
     float* d_vux_bak, * d_vuy_bak, * d_txx_bak, * d_tyy_bak, * d_txy_bak;
     float* d_SXxx_bak, * d_SXxy_bak, * d_SYxy_bak, * d_SYyy_bak;
     float* d_Vuxx_bak, * d_Vuxy_bak, * d_Vuyx_bak, * d_Vuyy_bak;
     float* d_pmlxSxx_bak, * d_pmlySxy_bak, * d_pmlxSxy_bak, * d_pmlySyy_bak;
     float* d_pmlxVux_bak, * d_pmlyVuy_bak, * d_pmlxVuy_bak, * d_pmlyVux_bak;
 
     float* d_vux_rs, * d_vuy_rs, * d_txx_rs, * d_tyy_rs, * d_txy_rs;
     float* d_SXxx_rs, * d_SXxy_rs, * d_SYxy_rs, * d_SYyy_rs;
     float* d_Vuxx_rs, * d_Vuxy_rs, * d_Vuyx_rs, * d_Vuyy_rs;
     float* d_pmlxSxx_rs, * d_pmlySxy_rs, * d_pmlxSxy_rs, * d_pmlySyy_rs;
     float* d_pmlxVux_rs, * d_pmlyVuy_rs, * d_pmlxVuy_rs, * d_pmlyVux_rs;
 
     float* d_vux_cur, * d_vuy_cur, * d_txx_cur, * d_tyy_cur, * d_txy_cur;
     float* d_SXxx_cur, * d_SXxy_cur, * d_SYxy_cur, * d_SYyy_cur;
     float* d_Vuxx_cur, * d_Vuxy_cur, * d_Vuyx_cur, * d_Vuyy_cur;
     float* d_pmlxSxx_cur, * d_pmlySxy_cur, * d_pmlxSxy_cur, * d_pmlySyy_cur;
     float* d_pmlxVux_cur, * d_pmlyVuy_cur, * d_pmlxVuy_cur, * d_pmlyVux_cur;
 
     float* d_vux_wrk, * d_vuy_wrk, * d_txx_wrk, * d_tyy_wrk, * d_txy_wrk;
     float* d_SXxx_wrk, * d_SXxy_wrk, * d_SYxy_wrk, * d_SYyy_wrk;
     float* d_Vuxx_wrk, * d_Vuxy_wrk, * d_Vuyx_wrk, * d_Vuyy_wrk;
     float* d_pmlxSxx_wrk, * d_pmlySxy_wrk, * d_pmlxSxy_wrk, * d_pmlySyy_wrk;
     float* d_pmlxVux_wrk, * d_pmlyVuy_wrk, * d_pmlxVuy_wrk, * d_pmlyVux_wrk;
 
 #define CMAL(p) cudaMalloc(&(p), mem2);
     CMAL(d_lam) CMAL(d_lam2) CMAL(d_rtx) CMAL(d_rty) CMAL(d_muxy)
         CMAL(d_vux) CMAL(d_vuy) CMAL(d_txx) CMAL(d_tyy) CMAL(d_txy)
         CMAL(d_pmlxSxx) CMAL(d_pmlySxy) CMAL(d_pmlxSxy) CMAL(d_pmlySyy)
         CMAL(d_SXxx) CMAL(d_SXxy) CMAL(d_SYxy) CMAL(d_SYyy)
         CMAL(d_pmlxVux) CMAL(d_pmlyVuy) CMAL(d_pmlxVuy) CMAL(d_pmlyVux)
         CMAL(d_Vuxx) CMAL(d_Vuxy) CMAL(d_Vuyx) CMAL(d_Vuyy)
         CMAL(d_dxi) CMAL(d_dxi2) CMAL(d_dyj) CMAL(d_dyj2)
         CMAL(d_e_dxi) CMAL(d_e_dxi2) CMAL(d_e_dyj) CMAL(d_e_dyj2)
         CMAL(d_vux_bak) CMAL(d_vuy_bak) CMAL(d_txx_bak) CMAL(d_tyy_bak) CMAL(d_txy_bak)
         CMAL(d_pmlxSxx_bak) CMAL(d_pmlySxy_bak) CMAL(d_pmlxSxy_bak) CMAL(d_pmlySyy_bak)
         CMAL(d_SXxx_bak) CMAL(d_SXxy_bak) CMAL(d_SYxy_bak) CMAL(d_SYyy_bak)
         CMAL(d_pmlxVux_bak) CMAL(d_pmlyVuy_bak) CMAL(d_pmlxVuy_bak) CMAL(d_pmlyVux_bak)
         CMAL(d_Vuxx_bak) CMAL(d_Vuxy_bak) CMAL(d_Vuyx_bak) CMAL(d_Vuyy_bak)
         CMAL(d_vux_rs) CMAL(d_vuy_rs) CMAL(d_txx_rs) CMAL(d_tyy_rs) CMAL(d_txy_rs)
         CMAL(d_pmlxSxx_rs) CMAL(d_pmlySxy_rs) CMAL(d_pmlxSxy_rs) CMAL(d_pmlySyy_rs)
         CMAL(d_SXxx_rs) CMAL(d_SXxy_rs) CMAL(d_SYxy_rs) CMAL(d_SYyy_rs)
         CMAL(d_pmlxVux_rs) CMAL(d_pmlyVuy_rs) CMAL(d_pmlxVuy_rs) CMAL(d_pmlyVux_rs)
         CMAL(d_Vuxx_rs) CMAL(d_Vuxy_rs) CMAL(d_Vuyx_rs) CMAL(d_Vuyy_rs)
         CMAL(d_vux_cur) CMAL(d_vuy_cur) CMAL(d_txx_cur) CMAL(d_tyy_cur) CMAL(d_txy_cur)
         CMAL(d_pmlxSxx_cur) CMAL(d_pmlySxy_cur) CMAL(d_pmlxSxy_cur) CMAL(d_pmlySyy_cur)
         CMAL(d_SXxx_cur) CMAL(d_SXxy_cur) CMAL(d_SYxy_cur) CMAL(d_SYyy_cur)
         CMAL(d_pmlxVux_cur) CMAL(d_pmlyVuy_cur) CMAL(d_pmlxVuy_cur) CMAL(d_pmlyVux_cur)
         CMAL(d_Vuxx_cur) CMAL(d_Vuxy_cur) CMAL(d_Vuyx_cur) CMAL(d_Vuyy_cur)
         CMAL(d_vux_wrk) CMAL(d_vuy_wrk) CMAL(d_txx_wrk) CMAL(d_tyy_wrk) CMAL(d_txy_wrk)
         CMAL(d_pmlxSxx_wrk) CMAL(d_pmlySxy_wrk) CMAL(d_pmlxSxy_wrk) CMAL(d_pmlySyy_wrk)
         CMAL(d_SXxx_wrk) CMAL(d_SXxy_wrk) CMAL(d_SYxy_wrk) CMAL(d_SYyy_wrk)
         CMAL(d_pmlxVux_wrk) CMAL(d_pmlyVuy_wrk) CMAL(d_pmlxVuy_wrk) CMAL(d_pmlyVux_wrk)
         CMAL(d_Vuxx_wrk) CMAL(d_Vuxy_wrk) CMAL(d_Vuyx_wrk) CMAL(d_Vuyy_wrk)
 #undef CMAL
 
     cudaMemcpy(d_lam, lam_e, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_lam2, lam2u_e, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_rtx, rho_tx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_rty, rho_ty, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_muxy, muxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_vux, h_vux, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_vuy, h_vuy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_txx, h_txx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_tyy, h_tyy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_txy, h_txy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlxSxx, h_pmlxSxx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlySxy, h_pmlySxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlxSxy, h_pmlxSxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlySyy, h_pmlySyy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_SXxx, h_SXxx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_SXxy, h_SXxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_SYxy, h_SYxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_SYyy, h_SYyy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlxVux, h_pmlxVux, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlyVuy, h_pmlyVuy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlxVuy, h_pmlxVuy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_pmlyVux, h_pmlyVux, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_Vuxx, h_Vuxx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_Vuxy, h_Vuxy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_Vuyx, h_Vuyx, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_Vuyy, h_Vuyy, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_dxi, dxi, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_dxi2, dxi2, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_dyj, dyj, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_dyj2, dyj2, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_e_dxi, e_dxi, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_e_dxi2, e_dxi2, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_e_dyj, e_dyj, mem2, cudaMemcpyHostToDevice);
     cudaMemcpy(d_e_dyj2, e_dyj2, mem2, cudaMemcpyHostToDevice);
 
     float* mains[N_STATE_FIELDS];
     float* baks[N_STATE_FIELDS];
     float* rss[N_STATE_FIELDS];
     float* curs[N_STATE_FIELDS];
     float* wrks[N_STATE_FIELDS];
     fill_state_ptrs(mains, d_vux, d_vuy, d_txx, d_tyy, d_txy, d_pmlxSxx, d_pmlySxy, d_pmlxSxy, d_pmlySyy,
         d_SXxx, d_SXxy, d_SYxy, d_SYyy, d_pmlxVux, d_pmlyVuy, d_pmlxVuy, d_pmlyVux, d_Vuxx, d_Vuxy, d_Vuyx, d_Vuyy);
     fill_state_ptrs(baks, d_vux_bak, d_vuy_bak, d_txx_bak, d_tyy_bak, d_txy_bak, d_pmlxSxx_bak, d_pmlySxy_bak, d_pmlxSxy_bak, d_pmlySyy_bak,
         d_SXxx_bak, d_SXxy_bak, d_SYxy_bak, d_SYyy_bak, d_pmlxVux_bak, d_pmlyVuy_bak, d_pmlxVuy_bak, d_pmlyVux_bak, d_Vuxx_bak, d_Vuxy_bak, d_Vuyx_bak, d_Vuyy_bak);
     fill_state_ptrs(rss, d_vux_rs, d_vuy_rs, d_txx_rs, d_tyy_rs, d_txy_rs, d_pmlxSxx_rs, d_pmlySxy_rs, d_pmlxSxy_rs, d_pmlySyy_rs,
         d_SXxx_rs, d_SXxy_rs, d_SYxy_rs, d_SYyy_rs, d_pmlxVux_rs, d_pmlyVuy_rs, d_pmlxVuy_rs, d_pmlyVux_rs, d_Vuxx_rs, d_Vuxy_rs, d_Vuyx_rs, d_Vuyy_rs);
     fill_state_ptrs(curs, d_vux_cur, d_vuy_cur, d_txx_cur, d_tyy_cur, d_txy_cur, d_pmlxSxx_cur, d_pmlySxy_cur, d_pmlxSxy_cur, d_pmlySyy_cur,
         d_SXxx_cur, d_SXxy_cur, d_SYxy_cur, d_SYyy_cur, d_pmlxVux_cur, d_pmlyVuy_cur, d_pmlxVuy_cur, d_pmlyVux_cur, d_Vuxx_cur, d_Vuxy_cur, d_Vuyx_cur, d_Vuyy_cur);
     fill_state_ptrs(wrks, d_vux_wrk, d_vuy_wrk, d_txx_wrk, d_tyy_wrk, d_txy_wrk, d_pmlxSxx_wrk, d_pmlySxy_wrk, d_pmlxSxy_wrk, d_pmlySyy_wrk,
         d_SXxx_wrk, d_SXxy_wrk, d_SYxy_wrk, d_SYyy_wrk, d_pmlxVux_wrk, d_pmlyVuy_wrk, d_pmlxVuy_wrk, d_pmlyVux_wrk, d_Vuxx_wrk, d_Vuxy_wrk, d_Vuyx_wrk, d_Vuyy_wrk);
 
     float best_dt = 6.f * H / (7.f * sqrtf(2.f) * Vpmax);
     if (DT >= best_dt)
         printf("时间步长过大，应小于 %f\n", best_dt);
     if (Vsmin / (F0 * H) < 15.f)
         printf("空间步长可能过大，注意频散\n");
     int n_picard = picard_iters_for_source(F0);
     const int env_np = picard_iters_from_env_override();
     if (env_np >= 0)
     {
         n_picard = env_np;
         printf("QN_N_PICARD 环境变量覆盖: n_picard=%d\n", n_picard);
     }
     printf("==== CrankSGFDM 宏 N_PICARD=%d  实际 n_picard=%d  上限 N_PICARD_MAX=%d  F0=%g",
         N_PICARD, n_picard, N_PICARD_MAX, (double)F0);
 #ifdef PICARD_FREQ_BOOST
     printf("  PICARD_FREQ_BOOST=ON");
 #else
     printf("  PICARD_FREQ_BOOST=OFF");
 #endif
     printf(" ====\n");
     printf("(若宏 N_PICARD 与你在 .cu 里改的数字不一致，说明 nvcc/CMake 带了 -DN_PICARD=… 覆盖了源文件)\n");
     if (n_picard == 0)
         printf("n_picard=0：跳过 Picard，每步在加源后仅做一次 R（与显式子步同构）。\n");
     printf("每步 FD_V+FD_T 共 %d 次（1×R(u^s) + Picard 内 %d×R）\n", 1 + n_picard, n_picard);
     const int picard_trace = (getenv("QN_PICARD_TRACE") != NULL && getenv("QN_PICARD_TRACE")[0] != '0');
 
     clock_t t0 = clock();
     for (int it = 0; it < NT; it++)
     {
         float tt = (float)it * DT;
         if (it % 100 == 0)
             printf("it=%d\n", it);
         float I_sou = -(1.f - 2.f * PI * PI * F0 * F0 * (tt - T0) * (tt - T0)) * expf(-PI * PI * F0 * F0 * (tt - T0) * (tt - T0));
         Source2D <<<GridSrc, BlockSrc>>> (d_txx, d_tyy, I_sou, sn);
         cudaDeviceSynchronize();
 
         copy_state_fields(baks, mains, mem2);
         copy_state_fields(wrks, baks, mem2);
         launch_VT<BLOCK_X, BLOCK_Y>(
             d_vux_wrk, d_vuy_wrk, d_rtx, d_rty, d_txx_wrk, d_tyy_wrk, d_txy_wrk,
             d_pmlxSxx_wrk, d_pmlySxy_wrk, d_pmlxSxy_wrk, d_pmlySyy_wrk,
             d_SXxx_wrk, d_SXxy_wrk, d_SYxy_wrk, d_SYyy_wrk,
             d_pmlxVux_wrk, d_muxy, d_pmlyVuy_wrk, d_pmlxVuy_wrk, d_pmlyVux_wrk,
             d_Vuxx_wrk, d_Vuxy_wrk, d_Vuyx_wrk, d_Vuyy_wrk,
             d_lam2, d_lam, d_e_dxi, d_dxi, d_e_dxi2, d_dxi2, d_e_dyj, d_dyj, d_e_dyj2, d_dyj2, DT, H, Grid1, Block1);
         cudaDeviceSynchronize();
         copy_state_fields(rss, wrks, mem2);
         copy_state_fields(mains, rss, mem2);
 
         if (n_picard > 0)
         {
             for (int pic = 0; pic < n_picard; pic++)
             {
                 copy_state_fields(curs, mains, mem2);
                 copy_state_fields(wrks, mains, mem2);
                 launch_VT<BLOCK_X, BLOCK_Y>(
                     d_vux_wrk, d_vuy_wrk, d_rtx, d_rty, d_txx_wrk, d_tyy_wrk, d_txy_wrk,
                     d_pmlxSxx_wrk, d_pmlySxy_wrk, d_pmlxSxy_wrk, d_pmlySyy_wrk,
                     d_SXxx_wrk, d_SXxy_wrk, d_SYxy_wrk, d_SYyy_wrk,
                     d_pmlxVux_wrk, d_muxy, d_pmlyVuy_wrk, d_pmlxVuy_wrk, d_pmlyVux_wrk,
                     d_Vuxx_wrk, d_Vuxy_wrk, d_Vuyx_wrk, d_Vuyy_wrk,
                     d_lam2, d_lam, d_e_dxi, d_dxi, d_e_dxi2, d_dxi2, d_e_dyj, d_dyj, d_e_dyj2, d_dyj2, DT, H, Grid1, Block1);
                 cudaDeviceSynchronize();
                 blend_all_fields(baks, rss, curs, wrks, mains, n2, blend_threads);
                 cudaDeviceSynchronize();
                 if (picard_trace && it == 0)
                 {
                     float txx_sn;
                     cudaMemcpy(&txx_sn, d_txx + sn, sizeof(float), cudaMemcpyDeviceToHost);
                     printf("  [QN_PICARD_TRACE] it=0 pic=%d txx[sn]=%e\n", pic, (double)txx_sn);
                 }
             }
         }
 
         cudaMemcpy(h_txx, d_txx, mem2, cudaMemcpyDeviceToHost);
         for (iy = 0; iy < NY_ext; iy++)
             sis[(size_t)iy * NT + it] = h_txx[IDX2(sx, iy, NX_ext)];
         if (it % 100 == 0)
         {
             char snap_name[64];
             snprintf(snap_name, sizeof(snap_name), "CrankLOPOMtxx_it%04d_2d.dat", it);
             wfile2d(snap_name, h_txx, NY_ext, NX_ext);
         }
     }
     printf("%f s\n", (float)(clock() - t0) / CLOCKS_PER_SEC);
     wfile2d("CrankLOPOMsisx_2d.dat", sis, NY_ext, NT);
     printf("2D grid NY_ext=%d NX_ext=%d NT=%d CN Picard=%d\n", NY_ext, NX_ext, NT, n_picard);
 
 #define CFREE(p) cudaFree(p);
     CFREE(d_lam) CFREE(d_lam2) CFREE(d_rtx) CFREE(d_rty) CFREE(d_muxy)
         CFREE(d_vux) CFREE(d_vuy) CFREE(d_txx) CFREE(d_tyy) CFREE(d_txy)
         CFREE(d_pmlxSxx) CFREE(d_pmlySxy) CFREE(d_pmlxSxy) CFREE(d_pmlySyy)
         CFREE(d_SXxx) CFREE(d_SXxy) CFREE(d_SYxy) CFREE(d_SYyy)
         CFREE(d_pmlxVux) CFREE(d_pmlyVuy) CFREE(d_pmlxVuy) CFREE(d_pmlyVux)
         CFREE(d_Vuxx) CFREE(d_Vuxy) CFREE(d_Vuyx) CFREE(d_Vuyy)
         CFREE(d_dxi) CFREE(d_dxi2) CFREE(d_dyj) CFREE(d_dyj2)
         CFREE(d_e_dxi) CFREE(d_e_dxi2) CFREE(d_e_dyj) CFREE(d_e_dyj2)
         CFREE(d_vux_bak) CFREE(d_vuy_bak) CFREE(d_txx_bak) CFREE(d_tyy_bak) CFREE(d_txy_bak)
         CFREE(d_pmlxSxx_bak) CFREE(d_pmlySxy_bak) CFREE(d_pmlxSxy_bak) CFREE(d_pmlySyy_bak)
         CFREE(d_SXxx_bak) CFREE(d_SXxy_bak) CFREE(d_SYxy_bak) CFREE(d_SYyy_bak)
         CFREE(d_pmlxVux_bak) CFREE(d_pmlyVuy_bak) CFREE(d_pmlxVuy_bak) CFREE(d_pmlyVux_bak)
         CFREE(d_Vuxx_bak) CFREE(d_Vuxy_bak) CFREE(d_Vuyx_bak) CFREE(d_Vuyy_bak)
         CFREE(d_vux_rs) CFREE(d_vuy_rs) CFREE(d_txx_rs) CFREE(d_tyy_rs) CFREE(d_txy_rs)
         CFREE(d_pmlxSxx_rs) CFREE(d_pmlySxy_rs) CFREE(d_pmlxSxy_rs) CFREE(d_pmlySyy_rs)
         CFREE(d_SXxx_rs) CFREE(d_SXxy_rs) CFREE(d_SYxy_rs) CFREE(d_SYyy_rs)
         CFREE(d_pmlxVux_rs) CFREE(d_pmlyVuy_rs) CFREE(d_pmlxVuy_rs) CFREE(d_pmlyVux_rs)
         CFREE(d_Vuxx_rs) CFREE(d_Vuxy_rs) CFREE(d_Vuyx_rs) CFREE(d_Vuyy_rs)
         CFREE(d_vux_cur) CFREE(d_vuy_cur) CFREE(d_txx_cur) CFREE(d_tyy_cur) CFREE(d_txy_cur)
         CFREE(d_pmlxSxx_cur) CFREE(d_pmlySxy_cur) CFREE(d_pmlxSxy_cur) CFREE(d_pmlySyy_cur)
         CFREE(d_SXxx_cur) CFREE(d_SXxy_cur) CFREE(d_SYxy_cur) CFREE(d_SYyy_cur)
         CFREE(d_pmlxVux_cur) CFREE(d_pmlyVuy_cur) CFREE(d_pmlxVuy_cur) CFREE(d_pmlyVux_cur)
         CFREE(d_Vuxx_cur) CFREE(d_Vuxy_cur) CFREE(d_Vuyx_cur) CFREE(d_Vuyy_cur)
         CFREE(d_vux_wrk) CFREE(d_vuy_wrk) CFREE(d_txx_wrk) CFREE(d_tyy_wrk) CFREE(d_txy_wrk)
         CFREE(d_pmlxSxx_wrk) CFREE(d_pmlySxy_wrk) CFREE(d_pmlxSxy_wrk) CFREE(d_pmlySyy_wrk)
         CFREE(d_SXxx_wrk) CFREE(d_SXxy_wrk) CFREE(d_SYxy_wrk) CFREE(d_SYyy_wrk)
         CFREE(d_pmlxVux_wrk) CFREE(d_pmlyVuy_wrk) CFREE(d_pmlxVuy_wrk) CFREE(d_pmlyVux_wrk)
         CFREE(d_Vuxx_wrk) CFREE(d_Vuxy_wrk) CFREE(d_Vuyx_wrk) CFREE(d_Vuyy_wrk)
 #undef CFREE
 
     free(lam_e); free(lam2u_e); free(mu_e);
     free(rho_tx); free(rho_ty); free(muxy);
     free(dxi); free(dxi2); free(dyj); free(dyj2);
     free(e_dxi); free(e_dxi2); free(e_dyj); free(e_dyj2);
     free(h_vux); free(h_vuy); free(h_txx); free(h_tyy); free(h_txy);
     free(h_pmlxSxx); free(h_pmlySxy); free(h_pmlxSxy); free(h_pmlySyy);
     free(h_SXxx); free(h_SXxy); free(h_SYxy); free(h_SYyy);
     free(h_pmlxVux); free(h_pmlyVuy); free(h_pmlxVuy); free(h_pmlyVux);
     free(h_Vuxx); free(h_Vuxy); free(h_Vuyx); free(h_Vuyy);
     free(sis);
     return 0;
 }
 
 void wfile2d(const char* fn, const float* data, int nrows, int ncols)
 {
     FILE* fp = fopen(fn, "wt");
     if (!fp) return;
     for (int i = 0; i < nrows; i++)
     {
         for (int j = 0; j < ncols; j++)
             fprintf(fp, "%e ", data[(size_t)i * (size_t)ncols + (size_t)j]);
         fprintf(fp, "\n");
     }
     fclose(fp);
 }
 
 void create_model_2d(float* vp, float* vs, float* rho, float* lamda, float* lamda2u, float* mu, int ny, int nx)
 {
     for (int iy = 0; iy < ny; iy++)
     {
         for (int ix = 0; ix < nx; ix++)
         {
             int o = IDX2(ix, iy, nx);
             if (ix > (nx / 2) - 10 && ix < (nx / 2) + 10)
             {
                 vp[o] = 1500.f; vs[o] = 0.f; rho[o] = 1000.f;
                 mu[o] = rho[o] * vs[o] * vs[o];
                 lamda2u[o] = rho[o] * vp[o] * vp[o];
                 lamda[o] = lamda2u[o] - 2.f * mu[o];
             }
             else if(iy > 200 && iy < 210)
             {
                 vp[o] = 6000.f; vs[o] = 3000.f; rho[o] = 2700.f;
                 mu[o] = rho[o] * vs[o] * vs[o];
                 lamda2u[o] = rho[o] * vp[o] * vp[o];
                 lamda[o] = lamda2u[o] - 2.f * mu[o];
             }
             else
             {
                 vp[o] = 5000.f; vs[o] = 3000.f; rho[o] = 2450.f;
                 mu[o] = rho[o] * vs[o] * vs[o];
                 lamda2u[o] = rho[o] * vp[o] * vp[o];
                 lamda[o] = lamda2u[o] - 2.f * mu[o];
             }
         }
     }
 }
 
 void extmodel2d_pad(const float* src, float* dst, int ny, int nx, int np)
 {
     int nx2 = nx + 2 * np;
     int ny2 = ny + 2 * np;
     for (int iy = 0; iy < ny2; iy++)
     {
         for (int ix = 0; ix < nx2; ix++)
         {
             int sx = ix - np;
             if (sx < 0) sx = 0;
             if (sx >= nx) sx = nx - 1;
             int sy = iy - np;
             if (sy < 0) sy = 0;
             if (sy >= ny) sy = ny - 1;
             dst[IDX2(ix, iy, nx2)] = src[IDX2(sx, sy, nx)];
         }
     }
 }
 
 