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

#include "initialize_leftlead.c"
it = 1;

while (it < 30)
{
	memset(E, 0, Nl * Nl * sizeof(float complex));
	memset(F, 0, Nl * Nl * sizeof(float complex));
	memset(E1, 0, Nl * Nl * sizeof(float complex));
	memset(F1, 0, Nl * Nl * sizeof(float complex));
	/*initialize Green*/

	/*initialize C:copy of Dd*/
	for (i = 0; i < Nl; i++)
	{
		for (j = 0; j < Nl; j++)
		{
			C[j * Nl + i] = Dd[j * Nl + i];
		}
	}

	M1 = Nl;
	N1 = Nl;
	lda1 = Nl;
	lwork1 = 33 * Nl;
	work1 = malloc(lwork1 * sizeof(float complex));
	ipiv1 = malloc(Nl * sizeof(int));
	cgetrf_(&M1, &N1, Dd, &lda1, ipiv1, &info1);
	cgetri_(&N1, Dd, &lda1, ipiv1, work1, &lwork1, &info1);
	free(work1);
	free(ipiv1);

	/*initialize E=AxDd(-1)xB=E1xB*/
	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				E1[k * Nl + i] += A[j * Nl + i] * Dd[k * Nl + j];
			}
		}
	}

	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				E[k * Nl + i] += E1[j * Nl + i] * B[k * Nl + j];
			}
		}
	}

	/*initialize F=BxDd(-1)xA=F1xA*/
	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				F1[k * Nl + i] += B[j * Nl + i] * Dd[k * Nl + j];
			}
		}
	}

	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				F[k * Nl + i] += F1[j * Nl + i] * A[k * Nl + j];
			}
		}
	}

	/*initialize d=d-E*/
	for (int i = 0; i < Nl; i++)
	{
		for (int j = 0; j < Nl; j++)
		{
			d[j * Nl + i] -= E[j * Nl + i];
		}
	}

	/*initialize Dd=C-E-F*/
	for (int i = 0; i < Nl; i++)
	{
		for (int j = 0; j < Nl; j++)
		{
			Dd[j * Nl + i] = C[j * Nl + i] - E[j * Nl + i] - F[j * Nl + i];
		}
	}

	memset(A1, 0, Nl * Nl * sizeof(float complex));

	/*initialize A=E1xA*/
	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				A1[k * Nl + i] += E1[j * Nl + i] * A[k * Nl + j];
			}
		}
	}

	for (i = 0; i < Nl; i++)
	{
		for (j = 0; j < Nl; j++)
		{
			A[j * Nl + i] = A1[j * Nl + i];
		}
	}

	memset(B1, 0, Nl * Nl * sizeof(float complex));

	/*initialize B=F1xB*/
	for (int i = 0; i < Nl; i++)
	{
		for (int k = 0; k < Nl; k++)
		{
			for (int j = 0; j < Nl; j++)
			{
				B1[k * Nl + i] += F1[j * Nl + i] * B[k * Nl + j];
			}
		}
	}

	for (i = 0; i < Nl; i++)
	{
		for (j = 0; j < Nl; j++)
		{
			B[j * Nl + i] = B1[j * Nl + i];
		}
	}

	it += 1;
}

M1 = Nl;
N1 = Nl;
lda1 = Nl;
lwork1 = 33 * Nl;
work1 = malloc(lwork1 * sizeof(float complex));
ipiv1 = malloc(N1 * sizeof(int));

cgetrf_(&M1, &N1, d, &lda1, ipiv1, &info1);
cgetri_(&N1, d, &lda1, ipiv1, work1, &lwork1, &info1);
free(work1);
free(ipiv1);

dl = malloc((Nl) * (Nl) * sizeof(float complex));

memset(dl, 0, Nl * Nl * sizeof(float complex));

for (i = 0; i < (Nl); i++)
{
	for (j = 0; j < (Nl); j++)
	{
		dl[j * Nl + i] = d[j * Nl + i];
	}
}

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
