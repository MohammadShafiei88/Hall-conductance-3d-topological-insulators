
for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_rightlead[2*(j*(Nl)+i)]=0.0;
		T1_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T2_rightlead[2*(j*(Nl)+i)]=0.0;
		T2_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		sigma_rightlead[2*(j*(Nl)+i)]=0.0;
		sigma_rightlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


/************ matrix A of right lead using for T2_rightlead **************/
for(i=0;i<Nl;i++)
	for(j=0;j<Nl;j++)
{
	iu = i/4;
	ju = j/4;

	ib = i%4;
	jb = j%4;

	dx = xx[i] - xx[j] - aa;	 /*  for right lead  */
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
		    {  T2_rightlead[2*(j*(Nl)+i)] += Ty[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1] += Ty[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  T2_rightlead[2*(j*(Nl)+i)] += Tx[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1] += Tx[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  T2_rightlead[2*(j*(Nl)+i)] += Tz[2*(jb*(4)+ib)];    T2_rightlead[2*(j*(Nl)+i)+1] += Tz[2*(jb*(4)+ib)+1];  }
	
}	

	
		

}							 /* end of matrix A using for T2_rightlead  */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_rightlead[2*(j*(Nl)+i)]=T2_rightlead[2*(i*(Nl)+j)];
		T1_rightlead[2*(j*(Nl)+i)+1]=-T2_rightlead[2*(i*(Nl)+j)+1];
	}
}




T1_dr=(double*)malloc((Nl)*(Nl)*2*sizeof(double));
T1_dl=(double*)malloc((Nl)*(Nl)*2*sizeof(double));

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_dr[2*(j*(Nl)+i)]=0.0;
		T1_dr[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			T1_dr[2*(j*(Nl)+i)]+=T1_rightlead[2*(k*(Nl)+i)]*dr[2*(j*(Nl)+k)]-T1_rightlead[2*(k*(Nl)+i)+1]*dr[2*(j*(Nl)+k)+1];
			T1_dr[2*(j*(Nl)+i)+1]+=T1_rightlead[2*(k*(Nl)+i)]*dr[2*(j*(Nl)+k)+1]+T1_rightlead[2*(k*(Nl)+i)+1]*dr[2*(j*(Nl)+k)];
		}
	}
}




for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			sigma_rightlead[2*(j*(Nl)+i)]+=T1_dr[2*(k*(Nl)+i)]*T2_rightlead[2*(j*(Nl)+k)]-T1_dr[2*(k*(Nl)+i)+1]*T2_rightlead[2*(j*(Nl)+k)+1];
			sigma_rightlead[2*(j*(Nl)+i)+1]+=T1_dr[2*(k*(Nl)+i)]*T2_rightlead[2*(j*(Nl)+k)+1]+T1_dr[2*(k*(Nl)+i)+1]*T2_rightlead[2*(j*(Nl)+k)];
		}
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_leftlead[2*(j*(Nl)+i)]=0.0;
		T1_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T2_leftlead[2*(j*(Nl)+i)]=0.0;
		T2_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		sigma_leftlead[2*(j*(Nl)+i)]=0.0;
		sigma_leftlead[2*(j*(Nl)+i)+1]=0.0;
	}

}


/************ matrix A left using as matrix T2_leftlead **************/
for(i=0;i<Nl;i++)
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
		    { T2_leftlead[2*(j*(Nl)+i)] += Txd[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1] += Txd[2*(jb*(4)+ib)+1];  }
		
	    if((dy2<1.01*aa*aa)&&(dy2>0.99*aa*aa))
		    {  T2_leftlead[2*(j*(Nl)+i)] += Tyd[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1] += Tyd[2*(jb*(4)+ib)+1];  }
		
	    if((dz2<1.01*aa*aa)&&(dz2>0.99*aa*aa))
		    {  T2_leftlead[2*(j*(Nl)+i)] += Tzd[2*(jb*(4)+ib)];    T2_leftlead[2*(j*(Nl)+i)+1] += Tzd[2*(jb*(4)+ib)+1];  }
	}	

}
	
								 /* end of matrix A right using as T2_leftlead */


/************************/

for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_leftlead[2*(j*(Nl)+i)]=T2_leftlead[2*(i*(Nl)+j)];
		T1_leftlead[2*(j*(Nl)+i)+1]=-T2_leftlead[2*(i*(Nl)+j)+1];
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		T1_dl[2*(j*(Nl)+i)]=0.0;
		T1_dl[2*(j*(Nl)+i)+1]=0.0;
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			T1_dl[2*(j*(Nl)+i)]+=T1_leftlead[2*(k*(Nl)+i)]*dl[2*(j*(Nl)+k)]-T1_leftlead[2*(k*(Nl)+i)+1]*dl[2*(j*(Nl)+k)+1];
			T1_dl[2*(j*(Nl)+i)+1]+=T1_leftlead[2*(k*(Nl)+i)]*dl[2*(j*(Nl)+k)+1]+T1_leftlead[2*(k*(Nl)+i)+1]*dl[2*(j*(Nl)+k)];
		}

	}
}




for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		for(k=0;k<(Nl);k++)
		{
			sigma_leftlead[2*(j*(Nl)+i)]+=T1_dl[2*(k*(Nl)+i)]*T2_leftlead[2*(j*(Nl)+k)]-T1_dl[2*(k*(Nl)+i)+1]*T2_leftlead[2*(j*(Nl)+k)+1];
			sigma_leftlead[2*(j*(Nl)+i)+1]+=T1_dl[2*(k*(Nl)+i)]*T2_leftlead[2*(j*(Nl)+k)+1]+T1_dl[2*(k*(Nl)+i)+1]*T2_leftlead[2*(j*(Nl)+k)];
		}
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		green[2*(ipc[j]*ND+ipc[i])]+=(-sigma_rightlead[2*(j*(Nl)+i)]);
		green[2*(ipc[j]*ND+ipc[i])+1]+=(-sigma_rightlead[2*(j*(Nl)+i)+1]);
	}
}


for(i=0;i<(Nl);i++)
{
	for(j=0;j<(Nl);j++)
	{
		green[2*(((j+ND-Nl))*ND+((i+ND-Nl)))]+=(-sigma_leftlead[2*(j*(Nl)+i)]);
		green[2*(((j+ND-Nl))*ND+((i+ND-Nl)))+1]+=(-sigma_leftlead[2*(j*(Nl)+i)+1]);
	}
}


free(T1_dr);
free(T1_dl);


