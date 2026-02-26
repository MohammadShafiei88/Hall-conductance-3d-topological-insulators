

for(i=0;i<ND;i++)
{
	for(j=0;j<ND;j++)
	{
		green[2*(j*ND+i)]=0.0;
		green[2*(j*ND+i)+1]=0.0;
	}
}


for(i=0;i<ND;i++)
{
	for(j=0;j<ND;j++)
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
		{  green[2*(j*(ND)+i)] += e-MM0[2*(jb*(4)+ib)]-zeeman_inplane[2*(jb*(4)+ib)]; 
   		   green[2*(j*(ND)+i)+1] += eta-MM0[2*(jb*(4)+ib)+1]-zeeman_inplane[2*(jb*(4)+ib)+1];
		
		if(xx[i]<(Nx/2)*aa)
		{   green[2*(j*(ND)+i)] += -zeeman_bottom[2*(jb*(4)+ib)];    
		    green[2*(j*(ND)+i)+1] += -zeeman_bottom[2*(jb*(4)+ib)+1]; }  
		if(xx[i]>0.999*(Nx/2)*aa)
		{   green[2*(j*(ND)+i)] += zeeman_top[2*(jb*(4)+ib)];    
		    green[2*(j*(ND)+i)+1] += zeeman_top[2*(jb*(4)+ib)+1]; }  
		}
	if((dd<1.01*aa)&&(dd>0.99*aa))
	{
    if(i>j)
    	{
	    if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Txd[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Tyd[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Tzdg[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tzdg[2*(jb*(4)+ib)+1];  }
		}
	
		if(j>i)
		{
	    
		if((dx2<1.01*aa*aa)&&(dx2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Tx[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Ty[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  green[2*(j*(ND)+i)] += -Tzg[2*(jb*(4)+ib)];    green[2*(j*(ND)+i)+1] += -Tzg[2*(jb*(4)+ib)+1];  }
		}
		}
	}	
}		
/* end of matrix green  */

/********************disordrer?*******************/

/* ALREADY COMMENTED OUT
	  for(i=0;i<ND;i++)
		 { green[2*((3*i)*ND+(3*i))]+=Er[i];
		  green[2*((3*i+1)*ND+(3*i+1))]+=Er[i];
		  green[2*((3*i+2)*ND+(3*i+2))]+=Er[i]; }
*/
/******************************** end **********************/

#include"add_sigma_to_green.c"

/* ALREADY COMMENTED OUT 
  printf("HDG matrix:\n");
	for(i=0;i<ND_T;i++)
	  {for(k=0;k<ND_T;k++)
		 { printf("%.4f  %u %u   %.5f\t%.5f\t   %.5f\t%.5f\n",e,i,k,green_up_up[2*(k*ND_T+i)],green_up_up[2*(k*ND_T+i)+1],green_down_down[2*(k*ND_T+i)],green_down_down[2*(k*ND_T+i)+1]);
		  }
	   }

printf("\n\n");
*/

lda1=ND;
N1=ND;
M1=ND;
lwork1=33*ND;
work1=(double *)malloc(lwork1*2*sizeof(double));
ipiv1=(int *)malloc(ND*sizeof(int));
//ZGETRF(&M1,&N1,green,&lda1,ipiv1,&info1);
//ZGETRI(&N1,green,&lda1,ipiv1,work1,&lwork1,&info1);
zgetrf_(&M1,&N1,green,&lda1,ipiv1,&info1);
zgetri_(&N1,green,&lda1,ipiv1,work1,&lwork1,&info1);
free(work1);
free(ipiv1);


