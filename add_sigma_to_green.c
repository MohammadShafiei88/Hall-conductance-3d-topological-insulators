memset(T1_rightlead, 0, Nl * Nl * sizeof(float complex));
memset(T2_rightlead, 0, Nl * Nl * sizeof(float complex));
memset(sigma_rightlead, 0, Nl * Nl * sizeof(float complex));

/************ matrix A of right lead using for T2_rightlead **************/
for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        int iu = i / 4;
        int ju = j / 4;

        int ib = i % 4;
        int jb = j % 4;

        float dx = xx[i] - xx[j] - aa; // for right lead
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
                T2_rightlead[j * Nl + i] += Ty[jb * 4 + ib];

            if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                T2_rightlead[j * Nl + i] += Tx[jb * 4 + ib];

            if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                T2_rightlead[j * Nl + i] += Tz[jb * 4 + ib];
        }
    }
} /* end of matrix A using for T2_rightlead  */

/************************/

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        T1_rightlead[j * Nl + i] = conjf(T2_rightlead[i * Nl + j]);
    }
}

T1_dr = malloc(Nl * Nl * sizeof(float complex));
T1_dl = malloc(Nl * Nl * sizeof(float complex));

memset(T1_dr, 0, Nl * Nl * sizeof(float complex));

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        for (int k = 0; k < Nl; k++)
        {
            T1_dr[j * Nl + i] += T1_rightlead[k * Nl + i] * dr[j * Nl + k];
        }
    }
}

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        for (int k = 0; k < Nl; k++)
        {
            sigma_rightlead[j * Nl + i] += T1_dr[k * Nl + i] * T2_rightlead[j * Nl + k];
        }
    }
}

memset(T1_leftlead, 0, Nl * Nl * sizeof(float complex));
memset(T2_leftlead, 0, Nl * Nl * sizeof(float complex));
memset(sigma_leftlead, 0, Nl * Nl * sizeof(float complex));

/************ matrix A left using as matrix T2_leftlead **************/
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

        float dd = sqrt(dx2 + dy2 + dz2);

        if (dd > 0.99f * aa && dd < 1.01f * aa)
        {
            if (dx2 > 0.99f * aa * aa && dx2 < 1.01f * aa * aa)
                T2_leftlead[j * Nl + i] += Txd[jb * 4 + ib];

            if (dy2 > 0.99f * aa * aa && dy2 < 1.01f * aa * aa)
                T2_leftlead[j * Nl + i] += Tyd[jb * 4 + ib];

            if (dz2 > 0.99f * aa * aa && dz2 < 1.01f * aa * aa)
                T2_leftlead[j * Nl + i] += Tzd[jb * 4 + ib];
        }
    }
}
/* end of matrix A right using as T2_leftlead */
/************************/

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        // transpose and take complex conjugate
        T1_leftlead[j * Nl + i] = conjf(T2_leftlead[i * Nl + j]);
    }
}

memset(T1_dl, 0, Nl * Nl * sizeof(float complex));

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        for (int k = 0; k < Nl; k++)
        {
            T1_dl[j * Nl + i] += T1_leftlead[k * Nl + i] * dl[j * Nl + k];
        }
    }
}

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        for (int k = 0; k < Nl; k++)
        {
            sigma_leftlead[j * Nl + i] += T1_dl[k * Nl + i] * T2_leftlead[j * Nl + k];
        }
    }
}

for (int i = 0; i < Nl; i++)
{
    for (int j = 0; j < Nl; j++)
    {
        green[ipc[j] * ND + ipc[i]] += -sigma_rightlead[j * Nl + i];
    }
}

for (i = 0; i < (Nl); i++)
{
    for (j = 0; j < (Nl); j++)
    {
        green[(j + ND - Nl) * ND + (i + ND - Nl)] += (-sigma_leftlead[j * Nl + i]);
    }
}

free(T1_dr);
free(T1_dl);
