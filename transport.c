#include <Accelerate/Accelerate.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

int main()
{
	int info1, lda1, *ipiv1, lwork1, n;
	int i, j, k, iu, ju, ib, jb, M1, N1, ip, *ipc;
	int it, flag, num_W, number, counter, w;
	float t, We, e, e0;
	float complex Tr;
	float delta_top, delta_bottom, vv;
	int Nx, Ny, Nz;
	float *xx, *yy, *zz, *xxp, *yyp, *zzp, dx, dy, dz, dx2, dy2, dz2, dd, dd2;
	float complex *gama_G1, *gama_G2, *gama_1, *gama_2, *green, *green_a;
	float complex *d, *dr, *dl, *Dd, *C, *A, *A1, *B, *B1, *E, *F, *E1, *F1, *work1;
	float complex *sigma_leftlead, *sigma_rightlead, *sigmastar_leftlead, *sigmastar_rightlead;
	float complex *T1_leftlead, *T2_leftlead, *T1_rightlead, *T2_rightlead, *T1_dr, *T1_dl;
	float complex *Tx, *Ty, *Tz, *Txd, *Tyd, *Tzd;
	float complex *v_sia, *MM0, *zeeman_top, *zeeman_bottom;

	const float pi = 3.1415926;
	const float eta = 0.000001;
	const int Nm = 4;

	/*========== Parameters ========== */
	const float aa = 1.0;
	const float AA = 0.5;
	const float BB = 0.25;
	const float MM = 0.28;

	Nx = 15;
	Ny = Nx;
	Nz = 3;

	delta_top = 0.15;
	delta_bottom = 0.15;

	const int NSyz = Ny * Nz;
	const int Nl = NSyz * 4;
	const int ND = Nx * Ny * Nz * 4;

	MM0 = malloc(Nm * Nm * sizeof(float complex));
	Tx = malloc(Nm * Nm * sizeof(float complex));
	Ty = malloc(Nm * Nm * sizeof(float complex));
	Tz = malloc(Nm * Nm * sizeof(float complex));
	Txd = malloc(Nm * Nm * sizeof(float complex));
	Tyd = malloc(Nm * Nm * sizeof(float complex));
	Tzd = malloc(Nm * Nm * sizeof(float complex));

	zeeman_top = malloc(Nm * Nm * sizeof(float complex));
	zeeman_bottom = malloc(Nm * Nm * sizeof(float complex));
	v_sia = malloc(Nm * Nm * sizeof(float complex));

#include "create_TT_matrices.c"

//	Print a hopping matrix for test
/*	printf("Txd =\n");
	for (int i = 0; i < 4; i++) {
    		for (int j = 0; j < 4; j++) {
        	float complex z = Txd[i + j*4];   
        	printf("(%7.3f %+7.3fi) ",
                crealf(z), cimagf(z));
    }
    printf("\n");
}
printf("\n");
*/
	
	xx = malloc(ND * sizeof(float));
	yy = malloc(ND * sizeof(float));
	zz = malloc(ND * sizeof(float));

	xxp = malloc(ND * sizeof(float));
	yyp = malloc(ND * sizeof(float));
	zzp = malloc(ND * sizeof(float));

	ipc = malloc(ND * sizeof(int));

	for (i = 0; i < ND; i++)
	{
		iu = i / 4;
		xx[i] = (iu / (Ny * Nz)) * aa;
		yy[i] = (iu % Ny) * aa;
		zz[i] = (iu % (Ny * Nz)) / Ny * aa;
	}

	for (i = 0; i < ND; i++)
	{
		iu = i / 4;
		yyp[i] = (iu / (Nx * Nz)) * aa;
		xxp[i] = (iu % Nx) * aa;
		zzp[i] = (iu % (Nx * Nz)) / Nx * aa;
	}

	for (i = 0; i < ND; i++)
		for (ip = 0; ip < ND; ip++)
		{
			dx = xx[i] - xxp[ip];
			dy = yy[i] - yyp[ip];
			dz = zz[i] - zzp[ip];
			dx2 = dx * dx;
			dy2 = dy * dy;
			dz2 = dz * dz;
			dd2 = dx2 + dy2 + dz2;
			if ((dd2 < 0.1 * aa * aa) && (i % 4 == ip % 4))
				ipc[ip] = i;
		}
	green = malloc(ND * ND * sizeof(float complex));

	gama_1 = malloc(ND * ND * sizeof(float complex));

	gama_2 = malloc(ND * ND * sizeof(float complex));

	gama_G1 = malloc(ND * ND * sizeof(float complex));

	gama_G2 = malloc(ND * ND * sizeof(float complex));

	sigma_leftlead = malloc(Nl * Nl * sizeof(float complex));

	sigma_rightlead = malloc(Nl * Nl * sizeof(float complex));

	sigmastar_leftlead = malloc(Nl * Nl * sizeof(float complex));

	sigmastar_rightlead = malloc(Nl * Nl * sizeof(float complex));

	T1_leftlead = malloc(Nl * Nl * sizeof(float complex));
	T2_leftlead = malloc(Nl * Nl * sizeof(float complex));
	T1_rightlead = malloc(Nl * Nl * sizeof(float complex));
	T2_rightlead = malloc(Nl * Nl * sizeof(float complex));

	for (e = 0.02; e <= 0.080001; e += 1.02)
	{
#include "Rightlead_green.c"
#include "Leftlead_green.c"

		memset(gama_1, 0, ND * ND * sizeof(float complex));
		memset(gama_2, 0, ND * ND * sizeof(float complex));

		memset(sigmastar_rightlead, 0, Nl * Nl * sizeof(float complex));
		memset(sigmastar_leftlead, 0, Nl * Nl * sizeof(float complex));

#include "initialize_green_D.c"

		for (int i = 0; i < Nl; i++)
		{
			for (int j = 0; j < Nl; j++)
			{
				sigmastar_rightlead[j * Nl + i] = conj(sigma_rightlead[i * Nl + j]);
				sigmastar_leftlead[j * Nl + i] = conj(sigma_leftlead[i * Nl + j]);
			}
		}

		for (i = 0; i < Nl; i++)
		{
			for (j = 0; j < Nl; j++)
			{
				int idx_s = j * Nl + i;
				int idx1 = ipc[j] * ND + ipc[i];
				int idx2 = (j + ND - Nl) * ND + (i + ND - Nl);

				gama_1[idx1] = I * (sigma_rightlead[idx_s] - sigmastar_rightlead[idx_s]);

				gama_2[idx2] = I * (sigma_leftlead[idx_s] - sigmastar_leftlead[idx_s]);
			}
		}

		free(dr);
		free(dl);

		memset(gama_G1, 0, ND * ND * sizeof(float complex));
		memset(gama_G2, 0, ND * ND * sizeof(float complex));

		for (i = 0; i < ND; i++)
		{
			for (k = 0; k < ND; k++)
			{
				for (j = 0; j < ND; j++)
				{
					gama_G1[k * ND + i] += gama_1[j * ND + i] * green[k * ND + j];
				}
			}
		}

		green_a = malloc(ND * ND * sizeof(float complex));
		memset(green_a, 0, ND * ND * sizeof(float complex));

		for (i = 0; i < ND; i++)
		{
			for (j = 0; j < ND; j++)
			{
				green_a[j * ND + i] = conj(green[i * ND + j]);
			}
		}

		for (int i = 0; i < ND; i++)
		{
			for (int k = 0; k < ND; k++)
			{
				for (int j = 0; j < ND; j++)
				{
					gama_G2[k * ND + i] += gama_2[j * ND + i] * green_a[k * ND + j];
				}
			}
		}

		free(green_a);

		Tr = 0.0 + 0.0 * I;

		for (int i = 0; i < ND; i++)
		{
			for (int j = 0; j < ND; j++)
			{
				Tr += gama_G1[j * ND + i] * gama_G2[i * ND + j];
			}
		}

		printf("%.3f\t%.6f\n", e, creal(Tr));

	} /* End of loop for energy e */

	free(xx);
	free(yy);
	free(zz);
	free(xxp);
	free(yyp);
	free(zzp);
	free(ipc);

	free(MM0);
	free(Tx);
	free(Ty);
	free(Tz);
	free(Txd);
	free(Tyd);
	free(Tzd);
	free(zeeman_top);
	free(zeeman_bottom);
	free(v_sia);

	free(gama_G1);
	free(gama_G2);
	free(gama_1);
	free(gama_2);
	free(green);
	free(sigma_leftlead);
	free(sigma_rightlead);
	free(sigmastar_leftlead);
	free(sigmastar_rightlead);
	free(T1_leftlead);
	free(T2_leftlead);
	free(T1_rightlead);
	free(T2_rightlead);

	return 0;
}
