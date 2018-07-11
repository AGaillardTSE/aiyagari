/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2018                     **/
/*************************************************/




// GET RESIDUALS OF THE HISTOGRAM //
void weighthist2(double x, double xgrid, double *residout, int *ixgrid, double vector[], int dimweight)
{
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        //printf("getresid: (ixgrid>(dimweight-1))");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        //printf("getresid: (ixgrid<0)");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    *residout=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
}




void SIMULATION(double *save, double *startdist, double *enddist, double *capitalout, double taxK, double *taxoutL)
{

    // INITIALIZATION //
	int i,y,k,*isave;
	double *distold,*distin,critdist,*saveres, totrevK, govexp, grossGDP, tempxgrid, distverif, *distasset;
    distold = (double *) calloc((ifulldim), sizeof(double));
	distin = (double *) calloc((ifulldim), sizeof(double));
    saveres = (double *) calloc((ifulldim), sizeof(double));
    isave = (int *) calloc((ifulldim), sizeof(int));
    distasset = (double *) calloc((maxigrid), sizeof(double));


    // Compute new distribution //
    for(y=0;y<maxygrid;y++)
    {
        for(i=0;i<maxigrid;i++)
        {
            distin[inx(i,y)]=startdist[inx(i,y)];
            if(save[inx(i,y)]>Gridmax && save[inx(i,y)]<Gridmax+0.00000001){save[inx(i,y)]=Gridmax-0.00000001;}
            if(save[inx(i,y)]<Gridmin && save[inx(i,y)]>Gridmin-0.00000001){save[inx(i,y)]=Gridmin;}
            tempxgrid = phiinv(save[inx(i,y)]);
            weighthist2(save[inx(i,y)],tempxgrid,&saveres[inx(i,y)],&isave[inx(i,y)],K,maxigrid);
        }
    }

	
    // STATIONARY DISTRIBUTION //
	critdist=1.0;

	while(critdist>epsdist)
	{


        // copy distin in distold //
        bascule(distin,distold,ifulldim);
        
        for(i=0;i<ifulldim;i++)
        {
            distin[i]=0.0;
        }
		
        // Compute new distribution //
		for(y=0;y<maxygrid;y++)
		{
			for(i=0;i<maxigrid;i++)
			{
				for(k=0;k<maxygrid;k++)
				{
					distin[inx(isave[inx(i,y)],k)] += ytrans[y][k]*(1.0-saveres[inx(i,y)])*(distold[inx(i,y)]);
					distin[inx((min((isave[inx(i,y)]+1),(maxigrid-1))),k)] += ytrans[y][k]*saveres[inx(i,y)]*(distold[inx(i,y)]);
				}
			}
		}

        // convergence criterion //
		critdist=0.0;
        distverif = 0.0;
        for(i=0;i<maxigrid;i++)
        {
            for(y=0;y<maxygrid;y++)
            {
                critdist=(max(critdist,fabs(distin[inx(i,y)]-(distold[inx(i,y)]))));
            }
		}
        
	}


    // Aggregate capital + Government revenue //
	*capitalout=0.0;
    totrevK=0.0;
	
	for(y=0;y<maxygrid;y++)
	{
		for(i=0;i<maxigrid;i++)
		{
			*capitalout += distin[inx(i,y)]*K[i];
            enddist[inx(i,y)] = distin[inx(i,y)];
            totrevK += distin[inx(i,y)]*rstar*K[i]*taxK;
		}
	}
    
    grossGDP = pow(*capitalout,alphapar)*pow(Labor, (1.0-alphapar)) + deltapar*(*capitalout);
    govexp = fracG*grossGDP;
    *taxoutL = (govexp - totrevK)/(wstar*Labor);
    
	
    free(distold);
    free(distin);
}
