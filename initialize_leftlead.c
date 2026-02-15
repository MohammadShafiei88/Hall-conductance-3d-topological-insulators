memset(d, 0, Nl * Nl * sizeof(float complex));
memset(Dd, 0, Nl * Nl * sizeof(float complex));

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        int iu = i / 4;
        int ju = j / 4;

        int ib = i % 4;
        int jb = j % 4;

        float dx = xx[i] - xx[j];
        float dy = yy[i] - yy[j];
        float dz = zz[i] - zz[j];

        float dx2 = dx * dx;
        float dy2 = dy * dy;
        float dz2 = dz * dz;

        float dd2 = dx2 + dy2 + dz2;
        float dd = sqrtf(dd2);

        // diagonal block
        if (iu == ju)
        {
            d[j * Nl + i] += (e - MM0[jb * 4 + ib]) + I * eta;

            // if (xx[i] < (Nx/2)*aa)
            d[j * Nl + i] += -zeeman_bottom[jb * 4 + ib];
            // if (xx[i] > 0.999*(Nx/2)*aa)
            //   d[j*Nl + i] += zeeman_top[jb*4 + ib];
        }

        // nearest-neighbor coupling
        if (dd < 1.01f * aa && dd > 0.99f * aa)
        {
            if (j < i)
            {
                if (dx2 < 1.01f * aa * aa && dx2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Txd[jb * 4 + ib];
                if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Tyd[jb * 4 + ib];
                if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Tzd[jb * 4 + ib];
            }
            else if (j > i)
            {
                if (dx2 < 1.01f * aa * aa && dx2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Tx[jb * 4 + ib];
                if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Ty[jb * 4 + ib];
                if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                    d[j * Nl + i] += -Tz[jb * 4 + ib];
            }
        }
    }
}
/* end of matrix d  */

for (i = 0; i < Nl; i++)
{
    for (j = 0; j < Nl; j++)
    {
        Dd[j * Nl + i] = d[j * Nl + i];
    }
}

memset(A, 0, Nl * Nl * sizeof(float complex));
memset(B, 0, Nl * Nl * sizeof(float complex));
memset(C, 0, Nl * Nl * sizeof(float complex));

/************ matrix A **************/
for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        int iu = i / 4;
        int ju = j / 4;

        int ib = i % 4;
        int jb = j % 4;

        float dx = xx[i] - xx[j] + aa; // for left lead
        float dy = yy[i] - yy[j];
        float dz = zz[i] - zz[j];

        float dx2 = dx * dx;
        float dy2 = dy * dy;
        float dz2 = dz * dz;

        float dd2 = dx2 + dy2 + dz2;
        float dd = sqrtf(dd2);

        // nearest-neighbor coupling
        if (dd < 1.01f * aa && dd > 0.99f * aa)
        {
            if (dx2 < 1.01f * aa * aa && dx2 > 0.99f * aa * aa)
                A[j * Nl + i] += Tx[jb * 4 + ib];

            if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                A[j * Nl + i] += Ty[jb * 4 + ib];

            if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                A[j * Nl + i] += Tz[jb * 4 + ib];
        }
    }
}

/************************/

for (i = 0; i < (Nl); i++)
{
    for (j = 0; j < (Nl); j++)
    {
        B[j * Nl + i] = conjf(A[i * Nl + j]);
    }
}
