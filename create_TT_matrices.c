memset(Tx, 0, Nm * Nm * sizeof(float complex));
memset(Ty, 0, Nm * Nm * sizeof(float complex));
memset(Tz, 0, Nm * Nm * sizeof(float complex));
memset(Txd, 0, Nm * Nm * sizeof(float complex));
memset(Tyd, 0, Nm * Nm * sizeof(float complex));
memset(Tzd, 0, Nm * Nm * sizeof(float complex));

memset(MM0, 0, Nm * Nm * sizeof(float complex));

memset(zeeman_top, 0, Nm * Nm * sizeof(float complex));
memset(zeeman_bottom, 0, Nm * Nm * sizeof(float complex));
memset(v_sia, 0, Nm * Nm * sizeof(float complex));

Tx[0] = BB;
Tx[12] = -I * 0.5 * AA;
Tx[5] = -BB;
Tx[9] = -I * 0.5 * AA;
Tx[6] = -I * 0.5 * AA;
Tx[10] = BB;
Tx[3] = -I * 0.5 * AA;
Tx[15] = -BB;

Ty[0] = BB;
Ty[12] = -0.5 * AA;
Ty[5] = -BB;
Ty[9] = -0.5 * AA;
Ty[6] = 0.5 * AA;
Ty[10] = BB;
Ty[3] = 0.5 * AA;
Ty[15] = -BB;

Tz[0] = BB;
Tz[4] = -I * 0.5 * AA;
Tz[1] = -I * 0.5 * AA;
Tz[5] = -BB;
Tz[10] = BB;
Tz[14] = I * 0.5 * AA;
Tz[11] = I * 0.5 * AA;
Tz[15] = -BB;

for (int i = 0; i < Nm; i++)
{
	for (int j = 0; j < Nm; j++)
	{
		Txd[j * Nm + i] = conj(Tx[i * Nm + j]);
		Tyd[j * Nm + i] = conj(Ty[i * Nm + j]);
		Tzd[j * Nm + i] = conj(Tz[i * Nm + j]);
	}
}

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

v_sia[0] = -vv;
v_sia[5] = vv;
v_sia[10] = vv;
v_sia[15] = -vv;
