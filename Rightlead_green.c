
d =(double *)malloc(Nl*Nl*2*sizeof(double));
Dd=(double *)malloc(Nl*Nl*2*sizeof(double));
								 /*C:copy of Dd*/
C =(double *)malloc(Nl*Nl*2*sizeof(double));
A =(double *)malloc(Nl*Nl*2*sizeof(double));
A1=(double *)malloc(Nl*Nl*2*sizeof(double));
B =(double *)malloc(Nl*Nl*2*sizeof(double));
B1=(double *)malloc(Nl*Nl*2*sizeof(double));
								 /*E=AxDd(-1)xB=E1xB*/
E =(double *)malloc(Nl*Nl*2*sizeof(double));
								 /*F=BxDd(-1)xA=F1xA*/
F =(double *)malloc(Nl*Nl*2*sizeof(double));
								 /*E1=AxDd(-1)*/
E1=(double *)malloc(Nl*Nl*2*sizeof(double));
								 /*F1=BxDd(-1)*/
F1=(double *)malloc(Nl*Nl*2*sizeof(double));

#include "initialize_rightlead.c"
it = 1;




//while(it<40)
while(it<30)
{

	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			E1[2*(j*Nl+i)]=0.0;
			E1[2*(j*Nl+i)+1]=0.0;
			E[2*(j*Nl+i)]=0.0;
			E[2*(j*Nl+i)+1]=0.0;
			F1[2*(j*Nl+i)]=0.0;
			F1[2*(j*Nl+i)+1]=0.0;
			F[2*(j*Nl+i)]=0.0;
			F[2*(j*Nl+i)+1]=0.0;
		}
	}
	/*initialize Green*/

	/*initialize C:copy of Dd*/
	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			C[2*(j*Nl+i)]=Dd[2*(j*Nl+i)];
			C[2*(j*Nl+i)+1]=Dd[2*(j*Nl+i)+1];
		}
	}
/*
if(it==20)
for(i=0;i<2*Nl*Nl;i++)
printf("%d\t%.6f\n",i,Dd[i]);
*/

	M1=Nl;
	N1=Nl;
	lda1=Nl;
	lwork1=33*Nl;
	work1=(double *)malloc(lwork1*2*sizeof(double));
	ipiv1=(int *)malloc(Nl*sizeof(int));
//	ZGETRF(&M1,&N1,Dd,&lda1,ipiv1,&info1);
//	ZGETRI(&N1,Dd,&lda1,ipiv1,work1,&lwork1,&info1);
	zgetrf_(&M1,&N1,Dd,&lda1,ipiv1,&info1);
	zgetri_(&N1,Dd,&lda1,ipiv1,work1,&lwork1,&info1);
	free(work1);
	free(ipiv1);



	/*initialize E=AxDd(-1)xB=E1xB*/
	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				E1[2*(k*Nl+i)]+=A[2*(j*Nl+i)]*Dd[2*(k*Nl+j)]-A[2*(j*Nl+i)+1]*Dd[2*(k*Nl+j)+1];
				E1[2*(k*Nl+i)+1]+=A[2*(j*Nl+i)]*Dd[2*(k*Nl+j)+1]+A[2*(j*Nl+i)+1]*Dd[2*(k*Nl+j)];
			}
		}
	}


	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				E[2*(k*Nl+i)]+=E1[2*(j*Nl+i)]*B[2*(k*Nl+j)]-E1[2*(j*Nl+i)+1]*B[2*(k*Nl+j)+1];
				E[2*(k*Nl+i)+1]+=E1[2*(j*Nl+i)]*B[2*(k*Nl+j)+1]+E1[2*(j*Nl+i)+1]*B[2*(k*Nl+j)];
			}
		}
	}

	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				F1[2*(k*Nl+i)]+=B[2*(j*Nl+i)]*Dd[2*(k*Nl+j)]-B[2*(j*Nl+i)+1]*Dd[2*(k*Nl+j)+1];
				F1[2*(k*Nl+i)+1]+=B[2*(j*Nl+i)]*Dd[2*(k*Nl+j)+1]+B[2*(j*Nl+i)+1]*Dd[2*(k*Nl+j)];
			}
		}
	}


	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				F[2*(k*Nl+i)]+=F1[2*(j*Nl+i)]*A[2*(k*Nl+j)]-F1[2*(j*Nl+i)+1]*A[2*(k*Nl+j)+1];
				F[2*(k*Nl+i)+1]+=F1[2*(j*Nl+i)]*A[2*(k*Nl+j)+1]+F1[2*(j*Nl+i)+1]*A[2*(k*Nl+j)];
			}
		}
	}

	/*initialize d=d-E*/
	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			d[2*(j*Nl+i)]=d[2*(j*Nl+i)]-E[2*(j*Nl+i)];
			d[2*(j*Nl+i)+1]=d[2*(j*Nl+i)+1]-E[2*(j*Nl+i)+1];
		}
	}

	/*initialize Dd=C-E-F*/
	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			Dd[2*(j*Nl+i)]=C[2*(j*Nl+i)]-E[2*(j*Nl+i)]-F[2*(j*Nl+i)];
			Dd[2*(j*Nl+i)+1]=C[2*(j*Nl+i)+1]-E[2*(j*Nl+i)+1]-F[2*(j*Nl+i)+1];
		}
	}
/*
if(it==38||it==39)
for(i=0;i<2*Nl*Nl;i++)
printf("%d\t%.6f\n",i,Dd[i]);


if(it==38||it==39)
for(i=0;i<2*Nl*Nl;i++)
printf("%d\t%.6f\n",i,A[i]);


if(it==38||it==39)
for(i=0;i<2*Nl*Nl;i++)
printf("%d\t%.6f\n",i,Dd[i]);
*/
	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			A1[2*(j*Nl+i)]=0.0;
			A1[2*(j*Nl+i)+1]=0.0;
		}
	}

	/*initialize A=E1xA*/
	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				A1[2*(k*Nl+i)]+=E1[2*(j*Nl+i)]*A[2*(k*Nl+j)]-E1[2*(j*Nl+i)+1]*A[2*(k*Nl+j)+1];
				A1[2*(k*Nl+i)+1]+=E1[2*(j*Nl+i)]*A[2*(k*Nl+j)+1]+E1[2*(j*Nl+i)+1]*A[2*(k*Nl+j)];
			}
		}
	}

	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			A[2*(j*Nl+i)]=A1[2*(j*Nl+i)];
			A[2*(j*Nl+i)+1]=A1[2*(j*Nl+i)+1];
		}
	}

	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			B1[2*(j*Nl+i)]=0.0;
			B1[2*(j*Nl+i)+1]=0.0;
		}
	}

	/*initialize B=F1xB*/
	for(i=0;i<Nl;i++)
	{
		for(k=0;k<Nl;k++)
		{
			for(j=0;j<Nl;j++)
			{
				B1[2*(k*Nl+i)]+=F1[2*(j*Nl+i)]*B[2*(k*Nl+j)]-F1[2*(j*Nl+i)+1]*B[2*(k*Nl+j)+1];
				B1[2*(k*Nl+i)+1]+=F1[2*(j*Nl+i)]*B[2*(k*Nl+j)+1]+F1[2*(j*Nl+i)+1]*B[2*(k*Nl+j)];
			}
		}
	}

	for(i=0;i<Nl;i++)
	{
		for(j=0;j<Nl;j++)
		{
			B[2*(j*Nl+i)]=B1[2*(j*Nl+i)];
			B[2*(j*Nl+i)+1]=B1[2*(j*Nl+i)+1];
		}
	}

	// printf("info1=%d\n\n",info1);
	// printf("lwork1=%u\n",lwork1);

	it+=1;
}

// return 0;	


M1=Nl;
N1=Nl;
lda1=Nl;
lwork1=33*Nl;
work1=(double *)malloc(lwork1*2*sizeof(double));
ipiv1=(int *)malloc(N1*sizeof(int));
//ZGETRF(&M1,&N1,d,&lda1,ipiv1,&info1);
//ZGETRI(&N1,d,&lda1,ipiv1,work1,&lwork1,&info1);
zgetrf_(&M1,&N1,d,&lda1,ipiv1,&info1);
zgetri_(&N1,d,&lda1,ipiv1,work1,&lwork1,&info1);
free(work1);
free(ipiv1);
dr=(double *)malloc((Nl)*(Nl)*2*sizeof(double));

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		dr[2*(j*(Nl)+i)]=0.0;
		dr[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		dr[2*(j*(Nl)+i)]=d[2*(j*Nl+i)];
		dr[2*(j*(Nl)+i)+1]=d[2*(j*Nl+i)+1];
	}
}


/*
	for(i=0;i<Nl;i++)
	{for(j=0;j<Nl;j++)
	printf("i=%d  , j=%d , dr_up_up[2*(%d*Nl+%d)]=%g ,   dr_up_up[2*(%d*Nl+%d)+1]=%g\n",i,j,j,i,dr_up_up[2*(j*Nl+i)],j,i,dr_up_up[2*(j*Nl+i)+1]);
	}
*/
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
