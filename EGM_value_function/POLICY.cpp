/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2018                     **/
/*************************************************/


// POLICY ITERATION //
void POLICY_EGM(double *VF, double *save, double *cons, double taxK, double taxL, int iter2)
{


// INTEGER //
int i,ii,y,ynext,iter,threshold_ii, Icase;


// INITIALIZATION //
double *Kendo, *VFnew, *Kendo_min, *eVF, *deVF, critere, weight, slope1, slope2, tempvf;
VFnew = (double *) calloc((ifulldim), sizeof(double));          // Value function on the next time grid, next iteration
Kendo = (double *) calloc((ifulldim), sizeof(double));       // endogenous grid values
Kendo_min = (double *) calloc((maxygrid), sizeof(double));       // endogenous grid values
eVF = (double *) calloc((ifulldim), sizeof(double));     // expected value function
deVF = (double *) calloc((ifulldim), sizeof(double));    // derivative of the expected value function

    
// START LOOP OVER DECISION RULES //
critere=1.0;
iter=0;


bascule(VF,eVF,ifulldim);


// FIRST TIME COMPUTATION OF EXPECTATION + THRESHOLD + ENDO GRID //
for(y = 0; y < maxygrid; y++){

    /** 1. ENDOGENOUS GRID POINT where the borrowing constraint is binding **/
    deVF[inx(0,y)] = (eVF[inx(1,y)]-eVF[inx(0,y)])/(K[1] - K[0]);
    Kendo_min[y] = (K[0] + inv_MU(deVF[inx(0,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK)); // this is the implied first asset level for which we achieve the borrowing constraint, so all grid point below should achieve the gridpoint (monotonicity of the value function)

 
    for(i = 0; i < maxigrid; i++){

        /** 2. COMPUTE THE IMPLIED CONSUMPTION LEVEL **/
        if(i == (maxigrid-1)){deVF[inx(i,y)] = (eVF[inx((maxigrid-1),y)]-eVF[inx((maxigrid-2),y)])/(K[(maxigrid-1)] - K[(maxigrid-2)]);}
        if(i > 0 && i < (maxigrid-1)){deVF[inx(i,y)] = deriv(eVF[inx((i-1),y)],eVF[inx(i,y)],eVF[inx((i+1),y)],K[(i-1)],K[i],K[(i+1)]);}
        Kendo[inx(i,y)] = (K[i] + inv_MU(deVF[inx(i,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));

    } // end igrid
} // end ygrid



#if PM == 0
while (critere > epsilon && iter < 1000)
{
#endif


    // MAIN LOOP //
    critere = 0.0;

    #if OMP == 1
    omp_set_dynamic(0);
    omp_set_num_threads(maxygrid);
   
    #pragma omp parallel
    {

    #pragma omp for private(y,i,ynext,ii,threshold_ii, Icase, weight, slope1, slope2, tempvf)
    #endif
    /** 3. INTERPOLATE THE VALUE FUNCTION AND COMPUTE CRITERION **/
    for(y = 0; y < maxygrid; y++){
    
        threshold_ii = 0;
        
        for(i = 0; i < maxigrid; i++){
            
            
            if(K[i] < Kendo_min[y]){
                save[inx(i,y)] = K[0];
                tempvf = eVF[inx(0,y)];
            }
            
            
            if(K[i] >= Kendo_min[y]){
            
                ii = max(threshold_ii,0);
                Icase = 0;
                
                while((K[i]>Kendo[inx(ii,y)]) && (ii < maxigrid)){
                    if(ii == (maxigrid-1)){Icase = 2;break;}else{ii++;}
                }

                if(Icase == 2){ // case where you extrapolate.
                    slope1 = (eVF[inx((maxigrid-1),y)] - eVF[inx((maxigrid-2),y)])/(Kendo[inx((maxigrid-1),y)] - Kendo[inx((maxigrid-2),y)]);
                    slope2 = (K[(maxigrid-1)] - K[(maxigrid-2)])/(Kendo[inx((maxigrid-1),y)] - Kendo[inx((maxigrid-2),y)]);
 
                    save[inx(i,y)] = (K[i] - Kendo[inx((maxigrid-1),y)])*slope2 + K[(maxigrid-1)];
                    tempvf = (K[i] - Kendo[inx((maxigrid-1),y)])*slope1 + eVF[inx((maxigrid-1),y)];
                }
                
                
                if(Icase == 0){ // normal case
                    weight = (K[i] - Kendo[inx((ii-1),y)])/(Kendo[inx(ii,y)] - Kendo[inx((ii-1),y)]);
                    
                    save[inx(i,y)] = inter1d(weight,K[(ii-1)],K[ii]);
                    tempvf = inter1d(weight,eVF[inx((ii-1),y)],eVF[inx(ii,y)]);
                }
                
                
                // save localisation of K[i] for next grid point //
                threshold_ii = ii; // for next iteration, then set to the previous solution.
            
            }
            
            cons[inx(i,y)] = prod[y]*wstar*(1-taxL) + K[i]*(1+rstar*(1.0-taxK)) - save[inx(i,y)];
            VFnew[inx(i,y)] = U(cons[inx(i,y)]) + tempvf;
            
            // COMPUTE CRITERION //
            critere = max(critere,fabs(VF[inx(i,y)] - VFnew[inx(i,y)]));
            
            // SAVE THE NEW VFI //
            VF[inx(i,y)] = VFnew[inx(i,y)];
            
        } // end igrid

    } // end ygrid
    
    
    
    #if OMP == 1
    #pragma omp for private (y,i,ynext)
    #endif
    for(y = 0; y < maxygrid; y++){
        for(i = 0; i < maxigrid; i++){
        
            eVF[inx(i,y)] = 0.0;
            for(ynext = 0; ynext < maxygrid; ynext++){
                eVF[inx(i,y)] += betapar*ytrans[y][ynext]*VF[inx(i,ynext)];
            }
            
            if(i >= 2){
                deVF[inx((i-1),y)] = deriv(eVF[inx((i-2),y)],eVF[inx((i-1),y)],eVF[inx((i),y)],K[(i-2)],K[(i-1)],K[(i)]);
                //deVF[inx((i-1),y)] = (eVF[inx((i),y)] - eVF[inx((i-2),y)])/(K[(i)]-K[(i-2)]);
                Kendo[inx((i-1),y)] = (K[(i-1)] + inv_MU(deVF[inx((i-1),y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));
            }
            
        }

        deVF[inx(0,y)] = (eVF[inx(1,y)]-eVF[inx(0,y)])/(K[1] - K[0]);
        Kendo_min[y] = (K[0] + inv_MU(deVF[inx(0,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));
        Kendo[inx(0,y)] = (K[0] + inv_MU(deVF[inx(0,y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));
        
        deVF[inx((maxigrid-1),y)] = (eVF[inx((maxigrid-1),y)]-eVF[inx((maxigrid-2),y)])/(K[(maxigrid-1)] - K[(maxigrid-2)]);
        Kendo[inx((maxigrid-1),y)] = (K[(maxigrid-1)] + inv_MU(deVF[inx((maxigrid-1),y)]) - wstar*prod[y]*(1-taxL))/(1+rstar*(1.0-taxK));
        
    }

 
    #if OMP == 1
    } // OMP
    #endif
    

    
    //printf("CONVERGENCE: %d, %20.15f\n", iter, critere);

    iter++;
   
#if PM == 0
}//end of while loop
#endif




free(Kendo_min);
free(eVF);
free(deVF);
free(VFnew);
free(Kendo);
	
}


