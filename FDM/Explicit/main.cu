#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define PI 3.1415926535f
#define NY 200
#define NX 200
#define NP 30
#define BLOCK_SIZE_X 8
#define BLOCK_SIZE_Y 16

/* 设备端与主机共用：扩边后 x 快变 */
#define NX_EXT (NX + 2 * NP)
#define NY_EXT (NY + 2 * NP)
#define IDX2(ix, iy, nxe) ((ix) + (iy) * (nxe))

int divUp(int a, int b) { return (a - 1) / b + 1; }

void wfile2d(const char* fn, const float* data, int nrows, int ncols);
void create_model_2d(float* vp, float* vs, float* rho, float* lamda, float* lamda2u, float* mu, int ny, int nx);
void extmodel2d_pad(const float* src, float* dst, int ny, int nx, int np);

//--------------------------------------------------------------------------
/* 平面 x–y 内爆炸源：txx, tyy, tzz 同加（平面应变仍保留 σzz 方程） */
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

/* 动量方程：仅 ∂x、∂y；z 向分裂场系数在主机置为 e=1,d=0，与三维公式兼容 */
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

	txx[o] += DT * (lamda_ext[o] * (uyy) + lamda2u_ext[o] * uxx);
	tyy[o] += DT * (lamda_ext[o] * (uxx) + lamda2u_ext[o] * uyy);
	txy[o] += muxy[o] * DT * (uxy + uyx);
}

//--------------------------------------------------------------------------
int main(void)
{
	const int NX_ext = NX_EXT;
	const int NY_ext = NY_EXT;
	const int n2 = NX_ext * NY_ext;
	const int BLOCK_X = BLOCK_SIZE_X;
	const int BLOCK_Y = BLOCK_SIZE_Y;

	int sx = NX_ext / 2;
	int sy = NY_ext / 2 - 40;
	int sn = IDX2(sx, sy, NX_ext);

	int NT = 4500;
	float H = 10.0f;
	float RC = 1.0e-6f;
	float DT = 1.0e-4f;
	float DP = (float)NP * H;
	float F0 = 30.0f;
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
			float m0 = mu_e[o], m1 = mu_e[IDX2(ix + 1, iy, NX_ext)], m2 = mu_e[IDX2(ix, iy + 1, NX_ext)], m3 = mu_e[IDX2(ix + 1, iy + 1, NX_ext)];
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
			int o = IDX2(ix, iy, NX_ext);
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
		HALLOC(h_txx) HALLOC(h_tyy) HALLOC(h_txy) HALLOC(h_tyz)
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

#define CMAL(p) cudaMalloc(&(p), mem2);
	CMAL(d_lam) CMAL(d_lam2) CMAL(d_rtx) CMAL(d_rty) CMAL(d_muxy)
		CMAL(d_vux) CMAL(d_vuy) CMAL(d_txx) CMAL(d_tyy) CMAL(d_txy)
		CMAL(d_pmlxSxx) CMAL(d_pmlySxy) CMAL(d_pmlxSxy) CMAL(d_pmlySyy)
		CMAL(d_SXxx) CMAL(d_SXxy) CMAL(d_SYxy) CMAL(d_SYyy)
		CMAL(d_pmlxVux) CMAL(d_pmlyVuy) CMAL(d_pmlxVuy) CMAL(d_pmlyVux)
		CMAL(d_Vuxx) CMAL(d_Vuxy) CMAL(d_Vuyx) CMAL(d_Vuyy)
		CMAL(d_dxi) CMAL(d_dxi2) CMAL(d_dyj) CMAL(d_dyj2)
		CMAL(d_e_dxi) CMAL(d_e_dxi2) CMAL(d_e_dyj) CMAL(d_e_dyj2)
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

	float best_dt = 6.f * H / (7.f * sqrtf(2.f) * Vpmax);
	if (DT >= best_dt)
		printf("时间步长过大，应小于 %f\n", best_dt);
	if (Vsmin / (F0 * H) < 15.f)
		printf("空间步长可能过大，注意频散\n");

	clock_t t0 = clock();
	for (int it = 0; it < NT; it++)
	{
		float tt = (float)it * DT;
		if (it % 100 == 0)
			printf("it=%d\n", it);
		float I_sou = -(1.f - 2.f * PI * PI * F0 * F0 * (tt - T0) * (tt - T0)) * expf(-PI * PI * F0 * F0 * (tt - T0) * (tt - T0));
		Source2D << <GridSrc, BlockSrc >> > (d_txx, d_tyy, I_sou, sn);
		cudaDeviceSynchronize();

		FD_V_2D<BLOCK_X, BLOCK_Y> << <Grid1, Block1 >> > (
			d_vux, d_vuy, d_rtx, d_rty, d_txx, d_tyy, d_txy,
			d_pmlxSxx, d_pmlySxy, d_pmlxSxy, d_pmlySyy,
			d_SXxx, d_SXxy, d_SYxy, d_SYyy,
			d_e_dxi, d_dxi, d_e_dxi2, d_dxi2, d_e_dyj, d_dyj, d_e_dyj2, d_dyj2, DT, H);
		cudaDeviceSynchronize();

		FD_T_2D<BLOCK_X, BLOCK_Y> << <Grid1, Block1 >> > (
			d_vux, d_vuy, d_txx, d_tyy, d_txy, d_lam2, d_lam,
			d_pmlxVux, d_muxy, d_pmlyVuy, d_pmlxVuy, d_pmlyVux,
			d_Vuxx, d_Vuxy, d_Vuyx, d_Vuyy,
			d_e_dxi, d_dxi, d_e_dxi2, d_dxi2, d_e_dyj2, d_dyj2, d_e_dyj, d_dyj, DT, H);
		cudaDeviceSynchronize();

		cudaMemcpy(h_txx, d_txx, mem2, cudaMemcpyDeviceToHost);
		for (iy = 0; iy < NY_ext; iy++)
			sis[(size_t)iy * NT + it] = h_txx[IDX2(sx, iy, NX_ext)];
		if (it % 100 == 0)
		{
			char snap_name[64];
			snprintf(snap_name, sizeof(snap_name), "LOPOMtxx_it%04d_2d.dat", it);
			wfile2d(snap_name, h_txx, NY_ext, NX_ext);
		}
	}
	printf("%f s\n", (float)(clock() - t0) / CLOCKS_PER_SEC);
	wfile2d("LOPOMsisx_2d.dat", sis, NY_ext, NT);
	printf("2D grid NY_ext=%d NX_ext=%d NT=%d\n", NY_ext, NX_ext, NT);

#define CFREE(p) cudaFree(p);
	CFREE(d_lam) CFREE(d_lam2) CFREE(d_rtx) CFREE(d_rty) CFREE(d_muxy)
		CFREE(d_vux) CFREE(d_vuy) CFREE(d_txx) CFREE(d_tyy) CFREE(d_txy)
		CFREE(d_pmlxSxx) CFREE(d_pmlySxy) CFREE(d_pmlxSxy) CFREE(d_pmlySyy)
		CFREE(d_SXxx) CFREE(d_SXxy) CFREE(d_SYxy) CFREE(d_SYyy)
		CFREE(d_pmlxVux) CFREE(d_pmlyVuy) CFREE(d_pmlxVuy) CFREE(d_pmlyVux)
		CFREE(d_Vuxx) CFREE(d_Vuxy) CFREE(d_Vuyx) CFREE(d_Vuyy)
		CFREE(d_dxi) CFREE(d_dxi2) CFREE(d_dyj) CFREE(d_dyj2)
		CFREE(d_e_dxi) CFREE(d_e_dxi2) CFREE(d_e_dyj) CFREE(d_e_dyj2)
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
			// if (((ix - nx / 2s) * (ix - nx / 2) + (iy - ny / 2) * (iy - ny / 2)) <= 100)
			if(iy > 100 && iy < 120)
             {
                 vp[o] = 6000.f; vs[o] = 3400.f; rho[o] = 2770.f;
                 mu[o] = rho[o] * vs[o] * vs[o];
                 lamda2u[o] = rho[o] * vp[o] * vp[o];
                 lamda[o] = lamda2u[o] - 2.f * mu[o];
             }
			else
			{
				vp[o] = 5000.f; vs[o] = 3300.f; rho[o] = 2450.f;
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
