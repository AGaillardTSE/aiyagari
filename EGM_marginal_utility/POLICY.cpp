/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2018                     **/
/*************************************************/



// POLICY ITERATION //
void POLICY_EGM(double *MU, double *save, double *Kendo, double *Kendo_min, double *cons, double taxK, double taxL, int iter2)
{


// INTEGER //
int i,ii,y,ynext,iter,threshold_ii, Icase;


// INITIALIZATION //
double *MUnew, *eMU, critere, weight, slope;
MUnew = (double *) calloc((ifulldim), sizeof(double));          // Value function on the next time grid, next iteration
eMU = (double *) calloc((ifulldim), sizeof(double));            // expected Marginal Utility



// START LOOP OVER DECISION RULES //
critere=1.0;
iter=0;

#if PM == 0
while (critere > epsilon && iter < 1000)
{
#endif


    #if OMP == 1
    omp_set_dynamic(0);
    omp_set_num_threads(maxygrid);
   
    #pragma omp parallel
    {

    #pragma omp for private(i,y,threshold_ii,ii,Icase,slope,weight)
    #endif


    // MAIN LOOP //
    /** 1. INTERPOLATE THE VALUE FUNCTION AND POLICY FUNCTION **/
    for(y = 0; y < maxygrid; y++){
    
        threshold_ii = 0;
        
        for(i = 0; i < maxigrid; i++){
            
            
            if(K[i] <= Kendo_min[y]){
                save[inx(i,y)] = K[0];
            }
            
            
            if(K[i] > Kendo_min[y]){
            
                ii = max(threshold_ii,0);
                Icase = 0;
                
                if((K[i]<Kendo[inx(ii,y)]+0.0000001) && ii == 0){Icase = 1;}
                
                while((K[i]>Kendo[inx(ii,y)]) && (ii < maxigrid)){
                    if(ii == (maxigrid-1)){Icase = 2;break;}else{ii++;}  // solution is after Yendo, so extrapolate forward
                }
                
                if(Icase == 2){ // case where you extrapolate.
                    slope = (K[(maxigrid-1)] - K[(maxigrid-2)])/(Kendo[inx((maxigrid-1),y)] - Kendo[inx((maxigrid-2),y)]);
                    save[inx(i,y)] = (K[i] - Kendo[inx((maxigrid-1),y)])*slope + K[(maxigrid-1)];
                }
                
                
                if(Icase == 0){ // normal case
                    weight = (K[i] - Kendo[inx((ii-1),y)])/(Kendo[inx(ii,y)] - Kendo[inx((ii-1),y)]);
                    save[inx(i,y)] = inter1d(weight,K[(ii-1)],K[ii]);
                    
                }

                if(Icase == 1){ // normal case
                    save[inx(i,y)] = K[0];
                }
                
                // save localisation of K[i] for next grid point //
                threshold_ii = ii; // for next iteration, then set to the previous solution.
            
            }
            
            cons[inx(i,y)] = K[i]*(1+rstar*(1-taxK)) + wstar*prod[y]*(1 - taxL) - save[inx(i,y)];
            MUnew[inx(i,y)] = MUc(cons[inx(i,y)]);

        } // end igrid
    } // end ygrid
 

    
    
    /** 2. COMPUTE CRITERION + EXPECTATION + ENDOGENOUS GRID **/
    critere=0.0;
    
    #if OMP == 1
    #pragma omp for private(i,y,ynext)
    #endif
    for(y = 0; y < maxygrid; y++){
        for(i = 0; i < maxigrid; i++){
            
            /** 2.1. COMPUTE THE CRITERION **/
            #if PM == 0
            critere = max(critere,fabs(MU[inx(i,y)] - MUnew[inx(i,y)]));
            #endif
            
            MU[inx(i,y)] = MUnew[inx(i,y)];
            
            /** 2.2. COMPUTE EXPECTATION **/
            eMU[inx(i,y)] = 0.0;
            for(ynext = 0; ynext < maxygrid; ynext++){
                eMU[inx(i,y)] += (1+rstar*(1.0-taxK))*betapar*ytrans[y][ynext]*MU[inx(i,ynext)];
            }
 
            /** 2.3. COMPUTE ENDOGENOUS GRID **/
            Kendo[inx(i,y)] = (K[i] + inv_MU(eMU[inx(i,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));

        }
        
        /** 2.4. COMPUTE ENDOGENOUS GRID THRESHOLD **/
        Kendo_min[y] = (K[0] + inv_MU(eMU[inx(0,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));
    }
    
    #if OMP == 1
    }
    #endif
    
    //printf("CONVERGENCE: %d, %20.15f\n", iter, critere);getchar();

    iter++;
   
#if PM == 0
}//end of while loop
#endif



free(eMU);
free(MUnew);

	
}


