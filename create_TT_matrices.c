
for(i=0;i<Nm;i++)
{
	for(j=0;j<Nm;j++)
	{
		Tx[2*(Nm*j+i)]    =0.0;
		Tx[2*(Nm*j+i)+1]  =0.0;
		Ty[2*(Nm*j+i)]    =0.0;
		Ty[2*(Nm*j+i)+1]  =0.0;
		Tz[2*(Nm*j+i)]    =0.0;
		Tz[2*(Nm*j+i)+1]  =0.0;	
		Tzg[2*(Nm*j+i)]    =0.0;
		Tzg[2*(Nm*j+i)+1]  =0.0;	

		zeeman_top[2*(Nm*j+i)]    =0.0;
		zeeman_top[2*(Nm*j+i)+1]  =0.0;
	        zeeman_bottom[2*(Nm*j+i)]    =0.0;
		zeeman_bottom[2*(Nm*j+i)+1]  =0.0;
		zeeman_inplane[2*(Nm*j+i)]    =0.0;
		zeeman_inplane[2*(Nm*j+i)+1]  =0.0;

	        Txd[2*(Nm*j+i)]    =0.0;
		Txd[2*(Nm*j+i)+1]  =0.0;
		Tyd[2*(Nm*j+i)]    =0.0;
		Tyd[2*(Nm*j+i)+1]  =0.0;
		Tzd[2*(Nm*j+i)]    =0.0;
		Tzd[2*(Nm*j+i)+1]  =0.0;
	
		Tzdg[2*(Nm*j+i)]    =0.0;
		Tzdg[2*(Nm*j+i)+1]  =0.0;
		
		MM0[2*(Nm*j+i)]    =0.0;
		MM0[2*(Nm*j+i)+1]  =0.0;

		v_sia[2*(Nm*j+i)]    =0.0;
		v_sia[2*(Nm*j+i)+1]  =0.0;
	}
}
Tx[0]  = BB;
Tx[25]  = -0.5*AA;
Tx[10] = -BB;
Tx[19] = -0.5*AA;
Tx[13] = -0.5*AA;
Tx[20] = BB;
Tx[7] = -0.5*AA;
Tx[30] = -BB;
Tx[2] = lambda*A0*cos(theta)/(0.255);
Tx[8] = lambda*A0*cos(theta)/(0.255);
Tx[22] =-lambda*A0*cos(theta)/(0.255);
Tx[28] = -lambda*A0*cos(theta)/(0.255);

Ty[0]  = BB;
Ty[24]  = -0.5*AA;
Ty[10] = -BB;
Ty[18] = -0.5*AA;
Ty[12] = 0.5*AA;
Ty[20] = BB;
Ty[6] = 0.5*AA;
Ty[30] = -BB;
Ty[3] = -lambda*A0*sin(theta)/(0.255);
Ty[9] = lambda*A0*sin(theta)/(0.255);
Ty[23] =-lambda*A0*sin(theta)/(0.255);
Ty[29] =lambda*A0*sin(theta)/(0.255);

Tz[0]  = BB;
Tz[9]  =-0.5*AAz;
Tz[3]  = -0.5*AAz;
Tz[10] = -BB;
Tz[20] = BB;
Tz[29] = 0.5*AAz;
Tz[23] = 0.5*AAz;
Tz[30] = -BB;



Tzg[0]  = BB;
Tzg[9]  = -zeta*0.5*AAz;
Tzg[3]  = -zeta*0.5*AAz;
Tzg[10] = -BB;
Tzg[20] = BB;
Tzg[29] = zeta*0.5*AAz;
Tzg[23] = zeta*0.5*AAz;
Tzg[30] = -BB;





Txd[0]  = BB;
Txd[25]  = 0.5*AA;
Txd[10] = -BB;
Txd[19] = 0.5*AA;
Txd[13] = 0.5*AA;
Txd[20] = BB;
Txd[7] = 0.5*AA;
Txd[30] = -BB;
Tx[2] = lambda*A0*cos(theta)/(0.255);
Tx[8] = lambda*A0*cos(theta)/(0.255);
Tx[22] =-lambda*A0*cos(theta)/(0.255);
Tx[28] = -lambda*A0*cos(theta)/(0.255);

Tyd[0]  = BB;
Tyd[24]  = 0.5*AA;
Tyd[10] = -BB;
Tyd[18] = 0.5*AA;
Tyd[12] = -0.5*AA;
Tyd[20] = BB;
Tyd[6] = -0.5*AA;
Tyd[30] = -BB;
Tyd[3] = -lambda*A0*sin(theta)/(0.255);
Tyd[9] = lambda*A0*sin(theta)/(0.255);
Tyd[23] = -lambda*A0*sin(theta)/(0.255);
Tyd[29] = lambda*A0*sin(theta)/(0.255);

Tzd[0]  = BB;
Tzd[9]  = 0.5*AAz;
Tzd[3]  = 0.5*AAz;
Tzd[10] = -BB;
Tzd[20] = BB;
Tzd[29] = -0.5*AAz;
Tzd[23] = -0.5*AAz;
Tzd[30] = -BB;


Tzdg[0]  = BB;
Tzdg[9]  = zeta*0.5*AAz;
Tzdg[3]  = zeta*0.5*AAz;
Tzdg[10] = -BB;
Tzdg[20] = BB;
Tzdg[29] = -zeta*0.5*AAz;
Tzdg[23] = -zeta*0.5*AAz;
Tzdg[30] = -BB;




MM0[0]  = MM-6*BB;
MM0[10] =-(MM-6*BB);
MM0[20] = MM-6*BB;
MM0[30] =-(MM-6*BB);


zeeman_top[0]  = delta_top;
zeeman_top[10] = delta_top;
zeeman_top[20] = -delta_top;
zeeman_top[30] = -delta_top;


zeeman_bottom[0]  = delta_bottom;
zeeman_bottom[10] = delta_bottom;
zeeman_bottom[20] = -delta_bottom;
zeeman_bottom[30] = -delta_bottom;

zeeman_inplane[4] =  Mx;
zeeman_inplane[5] =  My;
zeeman_inplane[14]=  Mx;
zeeman_inplane[15]=  My;
zeeman_inplane[16]=  Mx;
zeeman_inplane[17]= -My;
zeeman_inplane[26]=  Mx;
zeeman_inplane[27]= -My;


v_sia[0]= -vv;
v_sia[10]= vv;
v_sia[20]= vv;
v_sia[30]= -vv;
