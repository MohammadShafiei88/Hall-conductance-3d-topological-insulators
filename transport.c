//#include <Accelerate/Accelerate.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
int main()
{
								 /*argumans of ZGETRF & ZGETRI*/
	int info1,lda1,*ipiv1,lwork1,n;
	int i,j,k,iu,ju,ib,jb,M1,N1,ip,*ipc; /*i,j:for site l:number of sites&N1:width of lead*/
								 /*it:itration*/
	int it,flag,num_W,number,counter,w;
	double t,We,e,e0,pi,Tr_up_up_r,Tr_i,Tr_r;
	double *xx,*yy,*zz,*xxp,*yyp,*zzp,dx,dy,dz,dx2,dy2,dz2,dd,dd2;
	double *gama_G1,*gama_G2,*gama_1,*gama_2,*green,*green_a;
	int Nz;							 /*E1=AxDd(-1)=> E=E1xB & F1=BxDd(-1)=> F=F1xA*/
	double *d,*dr,*dl,*Dd,*C,*A,*A1,*B,*B1,*E,*F,*E1,*F1,*work1;
	double *sigma_leftlead,*sigma_rightlead,*sigmastar_leftlead,*sigmastar_rightlead;
	double *T1_leftlead,*T2_leftlead,*T1_rightlead,*T2_rightlead,*T1_dr,*T1_dl;
	double *v_sia,*MM0,*Tx,*Ty,*Tz,*Tzg,*zeeman_top,*zeeman_bottom,*zeeman_inplane,*Txd,*Tyd,*Tzd,*Tzdg;
	double Mx,My,AAz,zeta,A0,lambda,phi,theta;

	//   int random;
	//   long random_seed;
	//   double    *Er;
	//   FILE * fp;
	//   char file[40];
	//   double n_phi, y_n, lambda, B_field, plank_electron;
pi = 3.141592653589793238462643383279502884197;

	const double aa = 1.0;
	const double AA = 0.5;
	const double BB = 0.25;
	const double MM = 0.28;
//	const double delta_top = 0.15;
//   	const double delta_bottom = 0.15;
	//const double vv = 0.0;
	const double DD = 0.0;	
	int Nx,Ny;
	double delta_top,delta_bottom,vv;
	Nx= 20;
	Ny = Nx;
	AAz=0.5;

	lambda = 0.128/1000;
	//phi=pi/2;
	//theta=pi/6;

	zeta=1.0;
	for(phi=0.0; phi<=1.001*pi/2; phi+=10000*pi/60)
	{
	
	A0 = 0.0*sin(phi);
	delta_top = 0.1*cos(phi);
	delta_bottom = delta_top;
	vv=0.0;
	My = Mx = 0.0;
	
	for(theta=0.0;theta<=1.001*pi/2;theta+=1000*pi/6)
	{
	
	//const int Nx = 15;  /* Nx and Ny should be the same, otherwise Nl for right_lead should be changed */
	//const int Ny = Nx;
	Nz = 3;
	const int NSyz = Ny*Nz;
	const int Nl = NSyz*4;
	const int ND = Nx*Ny*Nz*4;

	const int Nm = 4;

	//    v=ND-ly;

	//    plank_electron=4.1356*100000;  /*Tesla * angestrom^2 */
	//    printf("%.5f\t",plank_electron);   /* ************in ra bepoers******* */
	//    B_field=0.0;  /*Tesla*/
	//    lambda=(-1/(plank_electron))*(2*pi)*(a*a*sqrt(3.0)/4.0)*B_field;
	//    printf("%.5f\t",lambda);   /*************in ra bepoers********/

	//    Lv=0.0;    //Electric field
	const double eta=0.0000001;
	//    num_W=40;   //the number of disorders
	//    number=35+1;   //the number of realization + 1 coloumn disorder
	//    R=1.0/35;   //for avarage
	//    printf("%.5f\n",R);

	// num_W=4;   //the number of disorders
	// number=2+1;   //the number of realization + 1 coloumn disorder
	// R=1.0;   //for avarage
	//We=0.0;     //disorder

	/*
		for(i=0;i<ND;i++)
		{
		iu = i/4;
		xx[i] = iu/(Ny*Nz)*sqrt(3.0)*0.5*aa;
		yy[i] = iu%Ny*aa + (iu/(Ny*Nz)+1)%2*0.5*aa;
		zz[i] = (iu%(Ny*Nz))/Ny*cc;
		}
	*/

	MM0=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tx=(double *)malloc(Nm*Nm*2*sizeof(double));
	Ty=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tz=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tzg=(double *)malloc(Nm*Nm*2*sizeof(double));
	Txd=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tyd=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tzd=(double *)malloc(Nm*Nm*2*sizeof(double));
	Tzdg=(double *)malloc(Nm*Nm*2*sizeof(double));
 	zeeman_inplane=(double *)malloc(Nm*Nm*2*sizeof(double));
 	zeeman_top=(double *)malloc(Nm*Nm*2*sizeof(double));
 	zeeman_bottom=(double *)malloc(Nm*Nm*2*sizeof(double));
	v_sia=(double *)malloc(Nm*Nm*2*sizeof(double));
	//Tzd=(double *)malloc(Nm*Nm*2*sizeof(double));
	


	#include "create_TT_matrices.c"



	xx=(double *)malloc(ND*sizeof(double));
	yy=(double *)malloc(ND*sizeof(double));
	zz=(double *)malloc(ND*sizeof(double));

	xxp=(double *)malloc(ND*sizeof(double));
	yyp=(double *)malloc(ND*sizeof(double));
	zzp=(double *)malloc(ND*sizeof(double));

	ipc=(int *)malloc(ND*sizeof(int));

	for(i=0;i<ND;i++)
	{
		iu = i/4;
		xx[i] = (iu/(Ny*Nz))*aa;
		yy[i] = (iu%Ny)*aa;
		zz[i] = (iu%(Ny*Nz))/Ny*aa;
	}

	for(i=0;i<ND;i++)
	{
		iu = i/4;
		yyp[i] = (iu/(Nx*Nz))*aa;
		xxp[i] = (iu%Nx)*aa;
		zzp[i] = (iu%(Nx*Nz))/Nx*aa;
	}

	for(i=0;i<ND;i++)
	for(ip=0;ip<ND;ip++)
	{
		dx = xx[i] - xxp[ip];
		dy = yy[i] - yyp[ip];
		dz = zz[i] - zzp[ip];
		dx2 = dx*dx; dy2 = dy*dy; dz2 = dz*dz;
		dd2 = dx2 + dy2 + dz2;
		if((dd2 < 0.1*aa*aa)&&(i%4==ip%4))   ipc[ip]=i; 
//		if((dd2 < 0.1*aa*aa)&&(i%4==ip%4))  { ipc[ip]=i; printf("%u\t%u\n",ip,i); }
// ipc[i]=i;
	}

	green=(double *)malloc(ND*ND*2*sizeof(double));
//	gama_1=(double *)malloc((Nl)*(Nl)*2*sizeof(double));
	gama_1=(double *)malloc((ND)*(ND)*2*sizeof(double));
//	gama_2=(double *)malloc((Nl)*(Nl)*2*sizeof(double));
	gama_2=(double *)malloc((ND)*(ND)*2*sizeof(double));
//	gama_G1=(double *)malloc((Nl)*ND*2*sizeof(double));
	gama_G1=(double *)malloc((ND)*ND*2*sizeof(double));
//	gama_G2=(double *)malloc((Nl)*ND*2*sizeof(double));
	gama_G2=(double *)malloc((ND)*ND*2*sizeof(double));
	sigma_leftlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	sigma_rightlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	sigmastar_leftlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	sigmastar_rightlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	T1_leftlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	T2_leftlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	T1_rightlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	T2_rightlead=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
	
//	e=0.02;
	for(e=0.02;e<=0.080001;e+=1.02)
{
	#include "Rightlead_green.c"
	#include "Leftlead_green.c"
        
	
	/*
	printf("dr_up_up matrix:\n");
		for(i=0;i<(Nl+2);i++)
		  {for(k=0;k<(Nl+2);k++)
			 { printf("%.4f\t%.4f\n",dr_up_up[2*(k*(Nl+2)+i)],dr_up_up[2*(k*(Nl+2)+i)+1]);
			   printf("%.4f\t%.4f\n",dl_up_up[2*(k*(Nl+2)+i)],dl_up_up[2*(k*(Nl+2)+i)+1]);
			  }
		   }

	*/

	for(i=0;i<(ND);i++)
	{
		for(j=0;j<(ND);j++)
		{
			gama_1[2*(j*(ND)+i)]=0.0;
			gama_2[2*(j*(ND)+i)]=0.0;
			gama_1[2*(j*(ND)+i)+1]=0.0;
			gama_2[2*(j*(ND)+i)+1]=0.0;
		}
	}

	for(i=0;i<(Nl);i++)


	{
		for(j=0;j<(Nl);j++)
		{
			sigmastar_rightlead[2*(j*(Nl)+i)]=0.0;
			sigmastar_rightlead[2*(j*(Nl)+i)+1]=0.0;
			sigmastar_leftlead[2*(j*(Nl)+i)]=0.0;
			sigmastar_leftlead[2*(j*(Nl)+i)+1]=0.0;
		}
	}

	#include "initialize_green_D.c"

//	printf("It Works correctly!!\n");
	for(i=0;i<(Nl);i++)
	{
		for(j=0;j<(Nl);j++)
		{
			sigmastar_rightlead[2*(j*(Nl)+i)]=sigma_rightlead[2*(i*(Nl)+j)];
			sigmastar_rightlead[2*(j*(Nl)+i)+1]=-sigma_rightlead[2*(i*(Nl)+j)+1];
			sigmastar_leftlead[2*(j*(Nl)+i)]=sigma_leftlead[2*(i*(Nl)+j)];
			sigmastar_leftlead[2*(j*(Nl)+i)+1]=-sigma_leftlead[2*(i*(Nl)+j)+1];
		}
	}

	/*
	 printf("sigma_up_up_right&leftlead:\n");
	   for(i=0;i<ND;i++)
		 {for(j=0;j<ND;j++)
			{    printf("%.4f   %u %u  %.5f\t%.5f\t  %.5f\t%.5f\n",e,i,j,sigma_up_up_rightlead[2*(j*ND_T+i)],sigma_up_up_rightlead[2*(j*ND_T+i)+1],sigmastar_up_up_rightlead[2*(j*ND_T+i)],sigmastar_up_up_rightlead[2*(j*ND_T+i)+1]);
			}
		 }

	*/

	for(i=0;i<(Nl);i++)
	{
		for(j=0;j<(Nl);j++)
		{
			gama_1[2*(ipc[j]*(ND)+ipc[i])]=-sigma_rightlead[2*(j*(Nl)+i)+1]+sigmastar_rightlead[2*(j*(Nl)+i)+1];
			gama_1[2*(ipc[j]*(ND)+ipc[i])+1]=sigma_rightlead[2*(j*(Nl)+i)]-sigmastar_rightlead[2*(j*(Nl)+i)];
		}
	}

	for(i=0;i<(Nl);i++)
	{
		for(j=0;j<(Nl);j++)
		{
			gama_2[2*((j+ND-Nl)*(ND)+(i+ND-Nl))]=-sigma_leftlead[2*(j*(Nl)+i)+1]+sigmastar_leftlead[2*(j*(Nl)+i)+1];
			gama_2[2*((j+ND-Nl)*(ND)+(i+ND-Nl))+1]=sigma_leftlead[2*(j*(Nl)+i)]-sigmastar_leftlead[2*(j*(Nl)+i)];
		}
	}

	/*
	   printf("gama_up_up_1,_2:\n");
	   for(i=0;i<ND_T;i++)
		 {for(j=0;j<ND;j++)
			{    printf("%.4f   %u %u  %.5f\t%.5f\t  %.5f\t%.5f\n",e,i,j,gama_up_up_1[2*(j*ND+i)],gama_up_up_1[2*(j*ND+i)+1],gama_up_up_2[2*(j*ND+i)],gama_up_up_2[2*(j*ND+i)+1]);
			}
		 }
	  */

	free(dr);
	free(dl);

	/*
	printf("green_HD matrix:\n");
		for(i=0;i<ND;i++)
		  {for(k=0;k<ND;k++)
			 { printf("%u %u  %.4f\t%.4f\n",i,k,green[2*(k*ND+i)],green[2*(k*ND+i)+1]);
			  }
		   }

	printf("\n\n");
	*/

	for(i=0;i<(ND);i++)
	{
		for(j=0;j<ND;j++)
		{
			gama_G1[2*(j*(ND)+i)]=0.0;
			gama_G1[2*(j*(ND)+i)+1]=0.0;
			gama_G2[2*(j*(ND)+i)]=0.0;
			gama_G2[2*(j*(ND)+i)+1]=0.0;
		}
	}

	for(i=0;i<(ND);i++)
	{
		for(k=0;k<ND;k++)
		{
			for(j=0;j<(ND);j++)
			{
				gama_G1[2*(k*(ND)+i)]+=gama_1[2*(j*(ND)+i)]*green[2*(k*ND+j)]-gama_1[2*(j*(ND)+i)+1]*green[2*(k*ND+j)+1];

				gama_G1[2*(k*(ND)+i)+1]+=gama_1[2*(j*(ND)+i)]*green[2*(k*ND+j)+1]+gama_1[2*(j*(ND)+i)+1]*green[2*(k*ND+j)];

			}
		}
	}

	green_a=(double *)malloc((ND)*ND*2*sizeof(double));

	for(i=0;i<(ND);i++)
	{
		for(j=0;j<ND;j++)
		{
			green_a[2*(j*(ND)+i)]=0.0;
			green_a[2*(j*(ND)+i)+1]=0.0;
		}
	}

	for(i=0;i<ND;i++)
	{
		for(j=0;j<ND;j++)
		{
			green_a[2*(j*ND+i)]=green[2*(i*ND+j)];
			green_a[2*(j*ND+i)+1]=-green[2*(i*ND+j)+1];
		}
	}
	for(i=0;i<(ND);i++)
	{
		for(k=0;k<ND;k++)
		{
			for(j=0;j<(ND);j++)
			{
				gama_G2[2*(k*(ND)+i)]+=gama_2[2*(j*(ND)+i)]*green_a[2*(k*(ND)+j)]-gama_2[2*(j*(ND)+i)+1]*green_a[2*(k*(ND)+j)+1];

				gama_G2[2*(k*(ND)+i)+1]+=gama_2[2*(j*(ND)+i)]*green_a[2*(k*(ND)+j)+1]+gama_2[2*(j*(ND)+i)+1]*green_a[2*(k*(ND)+j)];

			}
		}
	}

	free(green_a);

	Tr_r=Tr_i=0.0;
	for(i=0;i<(ND);i++)
		for(j=0;j<(ND);j++)
	{
		Tr_r+=gama_G1[2*(j*(ND)+i)]*gama_G2[2*(i*(ND)+j)]-gama_G1[2*(j*(ND)+i)+1]*gama_G2[2*(i*(ND)+j)+1];
		Tr_i+=gama_G1[2*(j*(ND)+i)]*gama_G2[2*(i*(ND)+j)+1]+gama_G1[2*(j*(ND)+i)+1]*gama_G2[2*(i*(ND)+j)];
	}


printf("%.6f\n",Tr_r);

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
	free(Tzg);
	free(Txd);
	free(Tyd);
	free(Tzd);
	free(Tzdg);
	free(zeeman_top);
	free(zeeman_bottom);
	free(zeeman_inplane);
	free(v_sia);
	//free(Tzd);

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
}//end of Mx loop
}//end of delta
	return 0;
}
