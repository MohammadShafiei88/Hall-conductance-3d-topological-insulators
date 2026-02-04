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
	float t, We, e, e0, pi, Tr_up_up_r, Tr_i, Tr_r;
	float *xx, *yy, *zz, *xxp, *yyp, *zzp, dx, dy, dz, dx2, dy2, dz2, dd, dd2;
	float complex *gama_G1, *gama_G2, *gama_1, *gama_2, *green, *green_a;
	int Nx, Ny, Nz;
	float complex *d, *dr, *dl, *Dd, *C, *A, *A1, *B, *B1, *E, *F, *E1, *F1, *work1;
	float complex *sigma_leftlead, *sigma_rightlead, *sigmastar_leftlead, *sigmastar_rightlead;
	float complex *T1_leftlead, *T2_leftlead, *T1_rightlead, *T2_rightlead, *T1_dr, *T1_dl;
	float complex *MM0, *Tx, *Ty, *Tz, *Txd, *Tyd, *Tzd;
	float complex *zeeman_top, *zeeman_bottom, *v_sia;
	float AA, BB, MM, delta_top, delta_bottom, aa;

	extern void cgetrf_(int* m, int* n, float complex* a, int* lda, int* ipiv, int* info);
	extern void cgetri_(int* n, float complex* a, int* lda, int* ipiv, float complex* work, int* lwork, int* info);
	extern void zgetrf_(int* m, int* n, double complex* a, int* lda, int* ipiv, int* info);
	extern void zgetri_(int* n, double complex* a, int* lda, int* ipiv, double complex* work, int* lwork, int* info);


	/* ========== Size of the system ========== */
	Nx = 15;
	Ny = Nx;
	Nz = 3;

	const int NSyz = Ny * Nz;
	const int Nl = NSyz * 4;
	const int ND = Nx * Ny * Nz * 4;
	int Nm = 4;

	const double eta = 0.0000001;

	/* ========== Model parameters ========== */
	AA = 0.50;
	BB = 0.25;
	MM = 0.30;
	aa = 1.0;
	delta_top = 0.1;
	delta_bottom = 0.1;

	MM0 = malloc(Nm * Nm * sizeof(float));
	Tx = malloc(Nm * Nm * sizeof(float complex));
	Ty = malloc(Nm * Nm * sizeof(float complex));
	Tz = malloc(Nm * Nm * sizeof(float complex));
	Txd = malloc(Nm * Nm * sizeof(float complex));
	Tyd = malloc(Nm * Nm * sizeof(float complex));
	Tzd = malloc(Nm * Nm * sizeof(float complex));
	v_sia = malloc(Nm * Nm * sizeof(float));
	zeeman_top = malloc(Nm * Nm * sizeof(float));
	zeeman_bottom = malloc(Nm * Nm * sizeof(float));

	memset(MM0, 0, Nm * Nm * sizeof(float));
	memset(Tx, 0, Nm * Nm * sizeof(float complex));
	memset(Ty, 0, Nm * Nm * sizeof(float complex));
	memset(Tz, 0, Nm * Nm * sizeof(float complex));
	memset(Txd, 0, Nm * Nm * sizeof(float complex));
	memset(Tyd, 0, Nm * Nm * sizeof(float complex));
	memset(Tzd, 0, Nm * Nm * sizeof(float complex));
	memset(zeeman_top, 0, Nm * Nm * sizeof(float));
	memset(zeeman_bottom, 0, Nm * Nm * sizeof(float));
	memset(v_sia, 0, Nm * Nm * sizeof(float));

	/* ========== Creating hopping matrices ========== */
	Tx[0] = BB;
	Tx[3] = -0.5 * AA * I;
	Tx[5] = -BB;
	Tx[6] = -0.5 * AA * I;
	Tx[9] = -0.5 * AA * I;
	Tx[10] = BB;
	Tx[12] = -0.5 * AA * I;
	Tx[15] = -BB;

	Ty[0] = BB;
	Ty[3] = 0.5 * AA;
	Ty[5] = -BB;
	Ty[6] = 0.5 * AA;
	Ty[9] = -0.5 * AA;
	Ty[10] = BB;
	Ty[12] = -0.5 * AA;
	Ty[15] = -BB;

	Tz[0] = BB;
	Tz[1] = -0.5 * AA * I;
	Tz[4] = -0.5 * AA * I;
	Tz[5] = -BB;
	Tz[10] = BB;
	Tz[11] = 0.5 * AA * I;
	Tz[14] = 0.5 * AA * I;
	Tz[15] = -BB;

	Txd[0] = BB;
	Txd[3] = 0.5 * AA * I;
	Txd[5] = -BB;
	Txd[6] = 0.5 * AA * I;
	Txd[9] = 0.5 * AA * I;
	Txd[10] = BB;
	Txd[12] = 0.5 * AA * I;
	Txd[15] = -BB;

	Tyd[0] = BB;
	Tyd[3] = -0.5 * AA;
	Tyd[5] = -BB;
	Tyd[6] = -0.5 * AA;
	Tyd[9] = 0.5 * AA;
	Tyd[10] = BB;
	Tyd[12] = 0.5 * AA;
	Tyd[15] = -BB;

	Tzd[0] = BB;
	Tzd[1] = 0.5 * AA * I;
	Tzd[4] = 0.5 * AA * I;
	Tzd[5] = -BB;
	Tzd[10] = BB;
	Tzd[11] = -0.5 * AA * I;
	Tzd[14] = -0.5 * AA * I;
	Tzd[15] = -BB;

	MM0[0] = MM - 6 * BB;
	MM0[5] = -(MM - 6 * BB);
	MM0[10] = MM - 6 * BB;
	MM0[15] = -(MM - 6 * BB);

	zeeman_top[0] = delta_top;
	zeeman_top[5] = delta_top;
	zeeman_top[10] = -delta_top;
	zeeman_top[15] = -delta_top;

	zeeman_bottom[0] = delta_bottom;
	zeeman_bottom[5] = delta_bottom;
	zeeman_bottom[10] = -delta_bottom;
	zeeman_bottom[15] = -delta_bottom;

	xx = (float *)malloc(ND * sizeof(float));
	yy = (float *)malloc(ND * sizeof(float));
	zz = (float *)malloc(ND * sizeof(float));

	xxp = (float *)malloc(ND * sizeof(float));
	yyp = (float *)malloc(ND * sizeof(float));
	zzp = (float *)malloc(ND * sizeof(float));

	ipc = (int *)malloc(ND * sizeof(int));

	double thresh = 0.1 * aa * aa;

	for (i = 0; i < ND; i++)
	{
		int imod = i & 3;

		double xi = xx[i];
		double yi = yy[i];
		double zi = zz[i];

		for (ip = 0; ip < ND; ip++)
		{

			if ((ip & 3) != imod)
				continue;

			double dx = xi - xxp[ip];
			double dy = yi - yyp[ip];
			double dz = zi - zzp[ip];

			double dd2 = dx * dx + dy * dy + dz * dz;

			if (dd2 < thresh)
				ipc[ip] = i;
		}
	}

	green = malloc(ND * ND * sizeof(float complex));

	gama_1 = malloc((ND) * (ND) * sizeof(float complex));

	gama_2 = malloc((ND) * (ND) * sizeof(float complex));

	gama_G1 = malloc((ND)*ND * sizeof(float complex));

	gama_G2 = malloc((ND)*ND * sizeof(float complex));

	sigma_leftlead = malloc((Nl) * (Nl) * sizeof(float complex));

	sigma_rightlead = malloc((Nl) * (Nl) * sizeof(float complex));

	sigmastar_leftlead = malloc((Nl) * (Nl) * sizeof(float complex));

	sigmastar_rightlead = malloc((Nl) * (Nl) * sizeof(float complex));

	T1_leftlead = malloc((Nl) * (Nl) * sizeof(float complex));

	T2_leftlead = malloc((Nl) * (Nl) * sizeof(float complex));

	T1_rightlead = malloc((Nl) * (Nl) * sizeof(float complex));

	T2_rightlead = malloc((Nl) * (Nl) * sizeof(float complex));

	for (e = 0.02; e <= 0.1000001; e += 1.01)
	{

		/* ========== Right lead ========== */

		d = malloc(Nl * Nl * sizeof(float complex));
		Dd = malloc(Nl * Nl * sizeof(float complex));
		/* C:copy of Dd  */
		C = malloc(Nl * Nl * sizeof(float complex));
		A = malloc(Nl * Nl * sizeof(float complex));
		A1 = malloc(Nl * Nl * sizeof(float complex));
		B = malloc(Nl * Nl * sizeof(float complex));
		B1 = malloc(Nl * Nl * sizeof(float complex));
		/* E=AxDd(-1)xB=E1xB */
		E = malloc(Nl * Nl * sizeof(float complex));
		/* F=BxDd(-1)xA=F1xA */
		F = malloc(Nl * Nl * sizeof(float complex));
		/* E1=AxDd(-1) */
		E1 = malloc(Nl * Nl * sizeof(float complex));
		/* F1=BxDd(-1) */
		F1 = malloc(Nl * Nl * sizeof(float complex));

		/* ========== Initialize right lead ========== */

		memset(d, 0, Nl * Nl * sizeof(float complex));
		memset(Dd, 0, Nl * Nl * sizeof(float complex));

		for (i = 0; i < Nl; i++)
		{
			int iu = i >> 2;
			int ib = i & 3;

			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (int j = 0; j < Nl; j++)
			{
				int ju = j >> 2;
				int jb = j & 3;

				float dx = xi - xx[j];
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd2 = dx2 + dy2 + dz2;
				float dd = sqrtf(dd2);

				int m = IDX(i, j);
				int s = SIDX(ib, jb);

				// ---- onsite term ----
				if (iu == ju)
				{
					d[m] += (e + eta * I) - MM0[s] - zeeman_bottom[s];
				}

				// ---- hopping ----
				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{

					if (j < i)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							d[m] += -Tyd[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							d[m] += -Txd[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							d[m] += -Tzd[s];
					}
					else if (j > i)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							d[m] += -Ty[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							d[m] += -Tx[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							d[m] += -Tz[s];
					}
				}
			}
		}

		memcpy(Dd, d, Nl * Nl * sizeof(float complex));

		memset(A, 0, Nl * Nl * sizeof(float complex));
		memset(B, 0, Nl * Nl * sizeof(float complex));
		memset(C, 0, Nl * Nl * sizeof(float complex));


		

		/************ matrix A **************/
#define IDX(i, j) ((j) * Nl + (i))
#define SIDX(i, j) ((j) * 4 + (i))

		for (int i = 0; i < Nl; i++)
		{
			int iu = i >> 2;
			int ib = i & 3;

			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (int j = 0; j < Nl; j++)
			{
				int ju = j >> 2;
				int jb = j & 3;

				float dx = xi - xx[j] - aa;
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd2 = dx2 + dy2 + dz2;
				float dd = sqrtf(dd2);

				int m = IDX(i, j);
				int s = SIDX(ib, jb);

				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
						A[m] += Tyd[s];
					if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
						A[m] += Txd[s];
					if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
						A[m] += Tzd[s];
				}
			}
		}

		for (int i = 0; i < Nl; i++)
		{
			for (int j = 0; j < Nl; j++)
			{
				int m = IDX(i, j);
				int n = IDX(j, i);
				B[m] = conjf(A[n]);
			}
		}

		/* ========== End of intialize_rightlead ========== */

		it = 1;

		while (it < 30)
		{

			memset(E1, 0, Nl * Nl * sizeof(float complex));
			memset(E, 0, Nl * Nl * sizeof(float complex));
			memset(F1, 0, Nl * Nl * sizeof(float complex));
			memset(F, 0, Nl * Nl * sizeof(float complex));

			memcpy(C, Dd, Nl * Nl * sizeof(float complex));

			int info;
			int *ipiv = malloc(Nl * sizeof(int));
			int lda = Nl;
			int lwork = 33 * Nl;
			float complex *work = malloc(lwork * sizeof(float complex));

			cgetrf_(&Nl, &Nl, Dd, &lda, ipiv, &info);
			cgetri_(&Nl, Dd, &lda, ipiv, work, &lwork, &info);

			free(ipiv);
			free(work);

			// ---- E1 = A x Dd ----
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						E1[IDX(k, i)] += A[IDX(j, i)] * Dd[IDX(k, j)];

			// ---- E = E1 x B ----
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						E[IDX(k, i)] += E1[IDX(j, i)] * B[IDX(k, j)];

			// ---- F1 = B x Dd ----
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						F1[IDX(k, i)] += B[IDX(j, i)] * Dd[IDX(k, j)];

			// ---- F = F1 x A ----
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						F[IDX(k, i)] += F1[IDX(j, i)] * A[IDX(k, j)];

			// ---- d = d - E ----
			for (i = 0; i < Nl * Nl; i++)
				d[i] -= E[i];

			// ---- Dd = C - E - F ----
			for (i = 0; i < Nl * Nl; i++)
				Dd[i] = C[i] - E[i] - F[i];

			// ---- A1 = E1 x A ----
			memset(A1, 0, Nl * Nl * sizeof(float complex));
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						A1[IDX(k, i)] += E1[IDX(j, i)] * A[IDX(k, j)];
			memcpy(A, A1, Nl * Nl * sizeof(float complex));

			// ---- B1 = F1 x B ----
			memset(B1, 0, Nl * Nl * sizeof(float complex));
			for (i = 0; i < Nl; i++)
				for (k = 0; k < Nl; k++)
					for (j = 0; j < Nl; j++)
						B1[IDX(k, i)] += F1[IDX(j, i)] * B[IDX(k, j)];
			memcpy(B, B1, Nl * Nl * sizeof(float complex));

			it++;
		}

		int *ipiv = malloc(Nl * sizeof(int));
		int lda = Nl;
		int lwork = 33 * Nl;
		float complex *work = malloc(lwork * sizeof(float complex));

		cgetrf_(&Nl, &Nl, d, &lda, ipiv, &info1);
		cgetri_(&Nl, d, &lda, ipiv, work, &lwork, &info1);

		free(work);
		free(ipiv);

		memset(dr, 0, Nl * Nl * sizeof(float complex));
		memcpy(dr, d, Nl * Nl * sizeof(float complex));

		free(E);
		free(F);
		free(E1);
		free(F1);
		free(d);
		free(Dd);
		free(A);
		free(A1);
		free(B);
		free(B1);
		free(C);

		/* ========== End of right lead ========== */

		/* ========== Left lead ========== */

		d = malloc(Nl * Nl * sizeof(float complex));
		Dd = malloc(Nl * Nl * sizeof(float complex));
		C = malloc(Nl * Nl * sizeof(float complex));
		A = malloc(Nl * Nl * sizeof(float complex));
		A1 = malloc(Nl * Nl * sizeof(float complex));
		B = malloc(Nl * Nl * sizeof(float complex));
		B1 = malloc(Nl * Nl * sizeof(float complex));
		E = malloc(Nl * Nl * sizeof(float complex));
		F = malloc(Nl * Nl * sizeof(float complex));
		E1 = malloc(Nl * Nl * sizeof(float complex));
		F1 = malloc(Nl * Nl * sizeof(float complex));

		memset(d, 0, Nl * Nl * sizeof(float complex));
		memset(Dd, 0, Nl * Nl * sizeof(float complex));
		memset(A, 0, Nl * Nl * sizeof(float complex));
		memset(B, 0, Nl * Nl * sizeof(float complex));
		memset(C, 0, Nl * Nl * sizeof(float complex));

		/* ========== Initialize left lead ========== */
		for (int i = 0; i < Nl; i++)
		{
			int iu = i / 4;
			int ib = i % 4;
			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (int j = 0; j < Nl; j++)
			{
				int ju = j / 4;
				int jb = j % 4;

				float dx = xi - xx[j];
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd = sqrtf(dx2 + dy2 + dz2);

				int m = IDX(i, j);
				int s = SIDX(ib, jb);

				// ---- onsite term ----
				if (iu == ju)
				{
					d[m] = (e + eta * I) - MM0[s] - zeeman_bottom[s];
				}

				// ---- hopping ----
				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (j < i)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							d[m] -= Txd[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							d[m] -= Tyd[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							d[m] -= Tzd[s];
					}
					else if (j > i)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							d[m] -= Tx[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							d[m] -= Ty[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							d[m] -= Tz[s];
					}
				}
			}
		}

		/* copy d -> Dd */
		memcpy(Dd, d, Nl * Nl * sizeof(float complex));

		/* ========== Initialize matrix A for left lead ========== */
		for (int i = 0; i < Nl; i++)
		{
			int iu = i / 4;
			int ib = i % 4;
			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (int j = 0; j < Nl; j++)
			{
				int ju = j / 4;
				int jb = j % 4;

				float dx = xi - xx[j] + aa; // left lead shift
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd = sqrtf(dx2 + dy2 + dz2);

				int m = IDX(i, j);
				int s = SIDX(ib, jb);

				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
						A[m] += Tx[s];
					if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
						A[m] += Ty[s];
					if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
						A[m] += Tz[s];
				}
			}
		}

		/* B = Hermitian conjugate of A */
		for (int i = 0; i < Nl; i++)
			for (int j = 0; j < Nl; j++)
				B[IDX(i, j)] = conjf(A[IDX(j, i)]);

		/* ========== End of initialize left lead ========== */
		it = 1;
		while (it < 30)
		{
			memset(E1, 0, Nl * Nl * sizeof(float complex));
			memset(E, 0, Nl * Nl * sizeof(float complex));
			memset(F1, 0, Nl * Nl * sizeof(float complex));
			memset(F, 0, Nl * Nl * sizeof(float complex));
			memcpy(C, Dd, Nl * Nl * sizeof(float complex));

			int info;
			int *ipiv = malloc(Nl * sizeof(int));
			int lda = Nl;
			int lwork = 33 * Nl;
			float complex *work = malloc(lwork * sizeof(float complex));

			cgetrf_(&Nl, &Nl, Dd, &lda, ipiv, &info);
			cgetri_(&Nl, Dd, &lda, ipiv, work, &lwork, &info);

			free(ipiv);
			free(work);

			/* E1 = A x Dd */
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						E1[IDX(k, i)] += A[IDX(j, i)] * Dd[IDX(k, j)];

			/* E = E1 x B */
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						E[IDX(k, i)] += E1[IDX(j, i)] * B[IDX(k, j)];

			/* F1 = B x Dd */
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						F1[IDX(k, i)] += B[IDX(j, i)] * Dd[IDX(k, j)];

			/* F = F1 x A */
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						F[IDX(k, i)] += F1[IDX(j, i)] * A[IDX(k, j)];

			/* d = d - E */
			for (int i = 0; i < Nl * Nl; i++)
				d[i] -= E[i];

			/* Dd = C - E - F */
			for (int i = 0; i < Nl * Nl; i++)
				Dd[i] = C[i] - E[i] - F[i];

			/* A1 = E1 x A */
			memset(A1, 0, Nl * Nl * sizeof(float complex));
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						A1[IDX(k, i)] += E1[IDX(j, i)] * A[IDX(k, j)];
			memcpy(A, A1, Nl * Nl * sizeof(float complex));

			/* B1 = F1 x B */
			memset(B1, 0, Nl * Nl * sizeof(float complex));
			for (int i = 0; i < Nl; i++)
				for (int k = 0; k < Nl; k++)
					for (int j = 0; j < Nl; j++)
						B1[IDX(k, i)] += F1[IDX(j, i)] * B[IDX(k, j)];
			memcpy(B, B1, Nl * Nl * sizeof(float complex));

			it++;
		}

		/* Final inverse on d */
		int *ipiv = malloc(Nl * sizeof(int));
		int lda = Nl;
		int lwork = 33 * Nl;
		float complex *work = malloc(lwork * sizeof(float complex));

		cgetrf_(&Nl, &Nl, d, &lda, ipiv, &info1);
		cgetri_(&Nl, d, &lda, ipiv, work, &lwork, &info1);

		free(ipiv);
		free(work);

		/* Copy d -> dl */
		dl = malloc(Nl * Nl * sizeof(float complex));
		memcpy(dl, d, Nl * Nl * sizeof(float complex));

		/* Free memory */
		free(E);
		free(F);
		free(E1);
		free(F1);
		free(d);
		free(Dd);
		free(A);
		free(A1);
		free(B);
		free(B1);
		free(C);

		/* ========== End of left lead ========== */

		// ---------- Initialize gamma matrices ----------
		memset(gama_1, 0, ND * ND * sizeof(float complex));
		memset(gama_2, 0, ND * ND * sizeof(float complex));

		// ---------- Initialize sigma star matrices ----------
		memset(sigmastar_rightlead, 0, Nl * Nl * sizeof(float complex));
		memset(sigmastar_leftlead, 0, Nl * Nl * sizeof(float complex));

		// ---------- Initialize Green matrix ----------
		memset(green, 0, ND * ND * sizeof(float complex));

		for (i = 0; i < ND; i++)
		{
			int iu = i / 4;
			int ib = i % 4;

			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (j = 0; j < ND; j++)
			{
				int ju = j / 4;
				int jb = j % 4;

				float dx = xi - xx[j];
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd2 = dx2 + dy2 + dz2;
				float dd = sqrtf(dd2);

				int m = j * ND + i;	 // linear index
				int s = jb * 4 + ib; // subblock index

				// ---- onsite term ----
				if (iu == ju)
				{
					green[m] += (e + eta * I) - MM0[s] - zeeman_bottom[s];
				}

				// ---- nearest neighbor hopping ----
				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (i > j)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							green[m] -= Txd[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							green[m] -= Tyd[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							green[m] -= Tzd[s];
					}
					else if (j > i)
					{
						if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
							green[m] -= Tx[s];
						if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
							green[m] -= Ty[s];
						if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
							green[m] -= Tz[s];
					}
				}
			}
		}

		/* ========== End of matrix green ==========  */

		// add_sigma_to_green

		memset(T1_rightlead, 0, Nl * Nl * sizeof(float complex));
		memset(T2_rightlead, 0, Nl * Nl * sizeof(float complex));
		memset(sigma_rightlead, 0, Nl * Nl * sizeof(float complex));

		for (i = 0; i < Nl; i++)
		{
			int iu = i / 4;
			int ib = i % 4;
			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (j = 0; j < Nl; j++)
			{
				int ju = j / 4;
				int jb = j % 4;

				float dx = xi - xx[j] - aa; // shift for right lead
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd = sqrtf(dx2 + dy2 + dz2);

				int s = jb * 4 + ib; // subblock index

				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
						T2_rightlead[j * Nl + i] += Ty[s];
					if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
						T2_rightlead[j * Nl + i] += Tx[s];
					if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
						T2_rightlead[j * Nl + i] += Tz[s];
				}
			}
		}

		// ---------- Compute T1_rightlead as conjugate transpose ----------
		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				T1_rightlead[j * Nl + i] = conjf(T2_rightlead[i * Nl + j]);

		// ---------- Compute T1_dr = T1_rightlead x dr ----------
		memset(T1_dr, 0, Nl * Nl * sizeof(float complex));

		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				for (int k = 0; k < Nl; k++)
					T1_dr[j * Nl + i] += T1_rightlead[k * Nl + i] * dr[j * Nl + k];

		// ---------- Compute sigma_rightlead = T1_dr x T2_rightlead ----------

		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				for (int k = 0; k < Nl; k++)
					sigma_rightlead[j * Nl + i] += T1_dr[k * Nl + i] * T2_rightlead[j * Nl + k];

		memset(T1_leftlead, 0, Nl * Nl * sizeof(float complex));
		memset(T2_leftlead, 0, Nl * Nl * sizeof(float complex));
		memset(sigma_leftlead, 0, Nl * Nl * sizeof(float complex));

		// ---------- Initialize T2_leftlead ----------

		for (i = 0; i < Nl; i++)
		{
			int iu = i / 4;
			int ib = i % 4;
			float xi = xx[i];
			float yi = yy[i];
			float zi = zz[i];

			for (j = 0; j < Nl; j++)
			{
				int ju = j / 4;
				int jb = j % 4;

				float dx = xi - xx[j] + aa; // shift for left lead
				float dy = yi - yy[j];
				float dz = zi - zz[j];

				float dx2 = dx * dx;
				float dy2 = dy * dy;
				float dz2 = dz * dz;
				float dd = sqrtf(dx2 + dy2 + dz2);

				int s = jb * 4 + ib; // subblock index

				if (dd > 0.99f * aa && dd < 1.01f * aa)
				{
					if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
						T2_leftlead[j * Nl + i] += Txd[s];
					if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
						T2_leftlead[j * Nl + i] += Tyd[s];
					if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
						T2_leftlead[j * Nl + i] += Tzd[s];
				}
			}
		}

		// ---------- Compute T1_leftlead as conjugate transpose ----------
		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				T1_leftlead[j * Nl + i] = conjf(T2_leftlead[i * Nl + j]);

		// ---------- Compute T1_dl = T1_leftlead x dl ----------
		memset(T1_dl, 0, Nl * Nl * sizeof(float complex));

		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				for (int k = 0; k < Nl; k++)
					T1_dl[j * Nl + i] += T1_leftlead[k * Nl + i] * dl[j * Nl + k];

		// ---------- Compute sigma_leftlead = T1_dl x T2_leftlead ----------

		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
				for (int k = 0; k < Nl; k++)
					sigma_leftlead[j * Nl + i] += T1_dl[k * Nl + i] * T2_leftlead[j * Nl + k];

		// ---------- Add lead self-energies into Green's function ----------
		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
			{
				green[ipc[j] * ND + ipc[i]] += -sigma_rightlead[j * Nl + i];
				green[(j + ND - Nl) * ND + (i + ND - Nl)] += -sigma_leftlead[j * Nl + i];
			}

		// end of add_sigma_to_green

		lda1 = ND;
		N1 = ND;
		M1 = ND;
			lwork1 = 33 * ND;
		work1 = (double complex *)malloc(lwork1 * sizeof(double complex));

		ipiv1 = (int *)malloc(ND * sizeof(int));
		zgetrf_(&M1, &N1, green, &lda1, ipiv1, &info1);
		zgetri_(&N1, green, &lda1, ipiv1, work1, &lwork1, &info1);
		free(work1);
		free(ipiv1);

		// end of initialize_green_D

		// ---------- Compute sigmastar matrices as conjugate transpose ----------
		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
			{
				sigmastar_rightlead[j * Nl + i] = conjf(sigma_rightlead[i * Nl + j]);
				sigmastar_leftlead[j * Nl + i] = conjf(sigma_leftlead[i * Nl + j]);
			}

		// ---------- Compute gamma matrices ----------
		for (i = 0; i < Nl; i++)
			for (j = 0; j < Nl; j++)
			{
				int gi = ipc[j];
				int gj = ipc[i];
				gama_1[gi * ND + gj] = -sigma_rightlead[j * Nl + i] + sigmastar_rightlead[j * Nl + i];
				gama_2[(j + ND - Nl) * ND + (i + ND - Nl)] = -sigma_leftlead[j * Nl + i] + sigmastar_leftlead[j * Nl + i];
			}

		// Free lead Green's function copies
		free(dr);
		free(dl);

		// ---------- Compute gama_G1 = gama_1 * green ----------
		for (i = 0; i < ND; i++)
			for (k = 0; k < ND; k++)
			{
				gama_G1[k * ND + i] = 0.0 + 0.0 * I;
				for (j = 0; j < ND; j++)
					gama_G1[k * ND + i] += gama_1[j * ND + i] * green[k * ND + j];
			}

		// ---------- Compute gama_G2 = gama_2 * green_a ----------
		float complex *green_a_c = malloc(ND * ND * sizeof(float complex));

		// Compute green_a as conjugate transpose of green
		for (i = 0; i < ND; i++)
			for (j = 0; j < ND; j++)
				green_a_c[j * ND + i] = conjf(green[i * ND + j]);

		for (i = 0; i < ND; i++)
			for (k = 0; k < ND; k++)
			{
				gama_G2[k * ND + i] = 0.0 + 0.0 * I;
				for (j = 0; j < ND; j++)
					gama_G2[k * ND + i] += gama_2[j * ND + i] * green_a_c[k * ND + j];
			}

		free(green_a_c);

		// ---------- Compute transmission Tr = Tr(gamma_G1 * gamma_G2) ----------
		float complex Tr = 0.0 + 0.0 * I;

		for (i = 0; i < ND; i++)
			for (j = 0; j < ND; j++)
				Tr += gama_G1[j * ND + i] * gama_G2[i * ND + j];

		printf("%.3f\t%.6f\n", e, crealf(Tr));

	} /* End of for e energy */

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
// end of transport.c
