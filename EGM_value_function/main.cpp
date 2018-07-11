/*************************************************/
/**          Alexandre GAILLARD - 2017          **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2018                     **/
/*************************************************/


// DECISION RULE: endogenous grid method [CAROLL, C. (2006)] //
// SIMULATION: use of histogram method //




#include "header.hpp"
#include "useful.cpp"
#include "tauchen.cpp"
#include "POLICY.cpp"
#include "SIMULATION.cpp"







// MAIN //
int main(int argc, char* argv[])
{

	
// INITIALIZATION //
int	i,y,k, iter;


timeval t1, t2;
double elapsedTime;



// TRANSITION + STATE Z OF ENTREPRENEURS //
tauchenfun(p_e, m_e, 0.0, std_e, prod, ytrans);
inv_distri(yinv, ytrans);

Labor = 0.0;
printf("PRODUCTIVITY SHOCKS\n");
for(int y = 0; y < maxygrid; y++){
    prod[y] = exp(prod[y]);
    Labor += yinv[y]*prod[y];
    printf("%f\t", prod[y]);
}
printf("\n");
printf("\n");
printf("TRANSITION MATRIX FOR PRODUCTIVITY SHOCKS\n");
for(int y = 0; y < maxygrid; y++){
    for(int k = 0; k < maxygrid; k++){
        printf("%f\t", ytrans[y][k]);
    }
    printf("\n");
}
printf("\n");




// GRID FOR ASSET //
for(i=0;i<maxigrid;i++)
{
	K[i]=phi(i);
}



// MARGINAL UTILITY, VALUES FUNCTION AND POLICIES //
double *VF, *save, *cons;          // for decision rules
double *distin,*distout;                        // for simulation
double capital1,capital0,PIB,critprice,taxL,welfare,taxoutL;    // for equilibrium


// Note for users :: please, always use pointers and save your computer's memory ;) == banish all arrays //
VF = (double *) calloc((ifulldim), sizeof(double));      // value function
save = (double *) calloc((ifulldim), sizeof(double));
cons = (double *) calloc((ifulldim), sizeof(double));

distin = (double *) calloc((ifulldim), sizeof(double));
distout = (double *) calloc((ifulldim), sizeof(double));



// GUESS PRICE & INITIAL DISTRIBUTION //
rstar = 0.041;          // interest rate
wstar = (1.0-alphapar)*pow(alphapar/(rstar+deltapar),alphapar/(1.0-alphapar));      // wage rate
distin[0] = 1.0;        // initial distribution
taxL = 0.3;             // initial taxes



// START BY GUESSING VF //
// I choose as an initial guess such that current asset level is two times the next asset level //
for(i = 0; i < maxigrid; i++){
    for(y = 0; y < maxygrid; y++){
        VF[inx(i,y)] = U(prod[y]*wstar + (1+rstar)*K[i])/(1+betapar);                          // REQUIERE TO BE INCREASING IN K (the case here)
    }
}



/** START EQUILIBRIUM FIXED POINT **/
// start timer
gettimeofday(&t1, NULL);


///// START BIG LOOP OVER INTEREST RATE /////
critprice=1.0;          // convergence criterion
iter = 0;               // iteration 0

while(critprice>epsprice)
{

    
    wstar = (1.0-alphapar)*pow(alphapar/(rstar+deltapar),alphapar/(1.0-alphapar));
    PIB = wstar*Labor/(1.0-alphapar);
    capital0 = Labor*pow(alphapar/(rstar+deltapar),1.0/(1.0-alphapar));
    
    //printf("Cnvg=%f, R=%f, K=%f, GDP=%f, L=%f, K/Y=%f, TaxK=%f\n", critprice, rstar, capital0, PIB, Labor, capital0/PIB, taxK);


    if(betapar*(1.0+rstar*(1.0-taxK))>1){printf("beta condition %20.15f\t%20.15f\n",rstar,betapar*(1.0+rstar*(1.0-taxK)));getchar();}
    

    // SOLVE POLICY FUNCTION //
    POLICY_EGM(VF,save,cons,taxK,taxL,iter);
    
    // SIMULATION //
    SIMULATION(save,distin,distout,&capital1,taxK,&taxoutL);
    bascule(distout,distin,ifulldim);


    // update prices //
    rstar=alphapar*pow(Labor,1.0-alphapar)/pow((relaxsK*capital1+(1.0-relaxsK)*capital0),1.0-alphapar)-deltapar;
    taxL=relaxsT*taxoutL+(1.0-relaxsT)*taxL;
    critprice=max(fabs((capital1-capital0)/capital0),fabs((taxoutL-taxL)/taxL));
    
    iter++;

} // end equilibrium loop.




gettimeofday(&t2, NULL);
elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;
elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;


// SAVE IN FILES //
FILE *pfile, *dfile, *vfile;


// POLICY FUNCTIONS //
pfile=fopen(policyfile, "w"); setbuf (pfile, NULL );
for(y=0;y<maxygrid;y++)
{
    for(i=0;i<maxigrid;i++)
    {
        fprintf(pfile,"%5d\t%20.15f\t%20.15f\n",i,phi(i),save[inx(i,y)]);
    }
    fprintf(pfile,"\n");
}
fclose(pfile);


// DISTRIBUTION //
dfile=fopen(distfile, "w"); setbuf (dfile, NULL );
for(i=0;i<maxigrid;i++)
{
    for(y=0;y<maxygrid;y++)
    {
        fprintf(dfile,"%5d\t%20.15f\t%20.15f\n",i,phi(i),distout[inx(i,y)]);
    }
    fprintf(dfile,"\n");
}
fclose(dfile);


// VALUE FUNCTION //
vfile=fopen(valuefile, "w"); setbuf (vfile, NULL );
for(i=0;i<maxigrid;i++)
{
    for(y=0;y<maxygrid;y++)
    {
        fprintf(vfile,"%5d\t%20.15f\t%20.15f\n",i,phi(i),VF[inx(i,y)]);
    }
    fprintf(vfile,"\n");
}
fclose(vfile);



printf("#########################################################################\n");
printf("TaxK %f (percent), TaxL: %20.15f, welfare: %f \n",taxK*100.0, taxoutL, welfare);
printf("Return to capital (percent) %20.10f, Net (percent) %20.10f\n",rstar*100.0,rstar*(1.0-taxK)*100.0);
printf("Aggregate saving rate (percent) %20.10f\n",(deltapar*alphapar)/(rstar+deltapar)*100.0);
printf("CRRA parameter %20.10f\n",rhopar);
printf("#########################################################################\n");





printf("===================Program ended in %f seconds\n",elapsedTime/1000);

return 0;

}


