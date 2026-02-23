memset(green, 0, ND * ND * sizeof(float complex));

for (int i = 0; i < ND; i++)
{
    for (int j = 0; j < ND; j++)
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

        if (iu == ju)
        {
            green[j * ND + i] += (e - MM0[jb * 4 + ib] + I * eta);

            // if (xx[i] < (Nx/2)*aa)
            green[j * ND + i] += -zeeman_bottom[jb * 4 + ib];
            // if (xx[i] > 0.999f*(Nx/2)*aa)
            //   green[j*ND + i] += zeeman_top[jb*4 + ib];
        }

        // nearest-neighbor coupling
        if (dd < 1.01f * aa && dd > 0.99f * aa)
        {
            if (i > j)
            {
                if (dx2 < 1.01f * aa * aa && dx2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Txd[jb * 4 + ib];
                if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Tyd[jb * 4 + ib];
                if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Tzd[jb * 4 + ib];
            }
            else if (j > i)
            {
                if (dx2 < 1.01f * aa * aa && dx2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Tx[jb * 4 + ib];
                if (dy2 < 1.01f * aa * aa && dy2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Ty[jb * 4 + ib];
                if (dz2 < 1.01f * aa * aa && dz2 > 0.99f * aa * aa)
                    green[j * ND + i] += -Tz[jb * 4 + ib];
            }
        }
    }
}

#include "add_sigma_to_green.c"

lda1 = ND;
N1 = ND;
M1 = ND;
lwork1 = 33 * ND;
work1 = malloc(lwork1 * sizeof(float complex));
ipiv1 = malloc(ND * sizeof(int));
cgetrf_(&M1, &N1, green, &lda1, ipiv1, &info1);
cgetri_(&N1, green, &lda1, ipiv1, work1, &lwork1, &info1);
free(work1);
free(ipiv1);
