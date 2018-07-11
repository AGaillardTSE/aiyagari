/*************************************************/
/**             Alexandre GAILLARD              **/
/**         AIYAGARI MODEL - EGM METHOD         **/
/**                    2018                     **/
/*************************************************/


#define OMP 1

/****************/
//    INCLUDE   //
/****************/
#include <iostream>
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <cstring>
#include <stdio.h>
#include <fstream>
#if OMP == 1
#include <omp.h>
#endif




// INDEX //
#define maxigrid 500 // define the grid of saving (next period wealth)
#define maxygrid 7

#define ifulldim (maxigrid*maxygrid)
#define inx(igridindex,jclassindex) (((jclassindex)*(maxigrid))+(igridindex))
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n))+(x0))

// GRID for assets //
const double Gridmin=0.0;
const double Gridmax=500.0;

const double Echelle1=1.6;
const double grmin=(Gridmin/Echelle1)-1.0;
const double Exponen=log((Gridmax/Echelle1)-grmin)/(maxigrid-1);
#define phi(x) ( Echelle1*(exp(Exponen*(x))+grmin) )
#define phiinv(x) (log((x)/Echelle1-grmin)/Exponen)

double K[maxigrid];




/***************************/
//  CALIBRATION DEFINITION //
/***************************/

// prices //
double wstar;
double rstar;

// Government steady-state //
const double fracG = 0.2; // 20% of gross market G/[f(K,N)+delta*K] = 0.2;
const double taxK = 0.35;

// Production parameters //
const double alphapar = 0.36; //production function parameter
const double deltapar = 0.08; //capital depreciation rate

// Preference parameters //
const double rhopar = 3.0; //CRRA parameter
const double betapar = 0.96; //discount factor

// Marginal utilities //
#define MUc(x) (pow((x),-rhopar))
#define inv_MU(u) (pow((u),(-(1/rhopar))))
#define U(x) (pow((x),(1.0-rhopar))/(1.0-rhopar))

// Income process //
const double p_e = 0.6;
const double std_e = 0.4;
const double m_e = 3;
double prod[maxygrid], ytrans[maxygrid][maxygrid], yinv[maxygrid], Labor;



//Convergence criterion
const double epsilon=	0.000001; //Convergence criterion on labor supply (and other stuff)
const double epsdist=	0.0000001; //Convergence criterion on stationary distribution
const double epsprice=	0.00001; //Convergence criterion on capital stock (interest rate)

//Relaxation parameters
#define relaxsK 0.5
#define relaxsT 0.5


//Output files
const char policyfile[]="policy.out";
const char distfile[]="dist.out";
const char valuefile[]="value.out";
