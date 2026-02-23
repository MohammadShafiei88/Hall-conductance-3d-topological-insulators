
for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		d[2*(j*(Nl)+i)]=0.0;
		d[2*(j*(Nl)+i)+1]=0.0;
		Dd[2*(j*(Nl)+i)]=0.0;
		Dd[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<Nl;i++)
{

	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j];
	dy = yy[i] - yy[j];
	dz = zz[i] - zz[j];

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	dd2 = dx2+dy2+dz2;
	dd = sqrt(dd2);

	
	if(iu==ju)
		{  d[2*(j*(Nl)+i)] += e-MM0[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += eta-MM0[2*(jb*(4)+ib)+1];  
		if(xx[i]<(Nx/2)*aa)
		{   d[2*(j*(Nl)+i)] += -zeeman_bottom[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -zeeman_bottom[2*(jb*(4)+ib)+1]; }  
		if(xx[i]>0.999*(Nx/2)*aa)
		{   d[2*(j*(Nl)+i)] += zeeman_top[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += zeeman_top[2*(jb*(4)+ib)+1]; }  
		}

	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	    if(j<i)
	    {
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Txd[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Tyd[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Tzd[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tzd[2*(jb*(4)+ib)+1];  }
		 }
		 if(j>i)   
		    {
		    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Tx[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Ty[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  d[2*(j*(Nl)+i)] += -Tz[2*(jb*(4)+ib)];    d[2*(j*(Nl)+i)+1] += -Tz[2*(jb*(4)+ib)+1];  }
		}
		}
	}								 
}   /* end of matrix d  */




for(i=0;i<Nl;i++)
{
	for(j=0;j<Nl;j++)
	{
		Dd[2*(j*Nl+i)]=d[2*(j*Nl+i)];
		Dd[2*(j*Nl+i)+1]=d[2*(j*Nl+i)+1];
	}
}


for(i=0;i<Nl;i++)
{
	for(j=0;j<Nl;j++)
	{
		A[2*(j*Nl+i)]=0.0;
		A[2*(j*Nl+i)+1]=0.0;
		B[2*(j*Nl+i)]=0.0;
		B[2*(j*Nl+i)+1]=0.0;
		C[2*(j*Nl+i)]=0.0;
		C[2*(j*Nl+i)+1]=0.0;
	}
}


/************ matrix A **************/
for(i=0;i<Nl;i++)
{
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j] + aa;	 /*  for left lead  */
	dy = yy[i] - yy[j];
	dz = zz[i] - zz[j];

	dx2 = dx*dx;
	dy2 = dy*dy;
	dz2 = dz*dz;

	dd2 = dx2+dy2+dz2;
	dd = sqrt(dd2);
	
	
if((dd<1.01*aa)&&(dd>0.99*aa))
	{
	
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  A[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1] += Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  A[2*(j*(Nl)+i)] += Ty[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1] += Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  A[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    A[2*(j*(Nl)+i)+1] += Tz[2*(jb*(4)+ib)+1]; } 
	}
		}
		}						 /* end of matrix A  */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)

	{
		B[2*(j*(Nl)+i)]=A[2*(i*(Nl)+j)];
		B[2*(j*(Nl)+i)+1]=-A[2*(i*(Nl)+j)+1];
	}
}


/*

   for(i=0;i<(Nl+2);i++)
	{for(j=0;j<(Nl+2);j++)
	  printf("i=%d  , j=%d , dr[2*(j*(Nl+2)+i)]=%.6f ,   dr[2*(j*(Nl+2)+i)+1]=%.6f ,  A[2*(j*(Nl+2)+i)]=%.6f  ,  B[2*(j*(Nl+2)+i)]=%.6f\n",i,j,d[2*(j*(Nl+2)+i)],d[2*(j*(Nl+2)+i)+1],A[2*(j*(Nl+2)+i)],B[2*(j*(Nl+2)+i)]);
	 }

  */

/*
printf("left lead\n\n");
   for(i=0;i<(Nl);i++)
        {for(j=0;j<(Nl);j++)
          printf("i=%d  , j=%d , dr_r=%.6f , dr_i=%.6f , A_r=%.6f , A_i=%.6f , B_r=%.6f , B_i=%.6f\n",i,j,d[2*(j*(Nl)+i)],d[2*(j*(Nl)+i)+1],A[2*(j*(Nl)+i)],A[2*(j*(Nl)+i)+1],B[2*(j*(Nl)+i)],B[2*(j*(Nl)+i)+1]);
         }

printf("\n\n");
*/
