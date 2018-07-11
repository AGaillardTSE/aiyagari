/*********************/
/** Useful function **/
/*********************/



#define max(a,b) (((a)>(b))?(a):(b))
#define min(a,b) (((a)<(b))?(a):(b))
#define interpol(x,y,z) (y+(x-floor(x))*(z-y))
#define inter1d(x1,y1,y2) ((1.0-(x1))*(y1)+(x1)*(y2))
#define inter2d(x1,x2,y11,y21,y12,y22) ((1.0-(x2))*((1.0-(x1))*(y11)+(x1)*(y21))+(x2)*((1.0-(x1))*(y12)+(x1)*(y22)))


// TO DEFINE GRID SPACING //
#define linspace(x0,xmax,n,i) ((i)*(((xmax)-(x0))/(n))+(x0)) // transform a grid into a value
#define invlinspace(x0,xmax,n,x) ((((x)-(x0))/((xmax)-(x0)))*(n)) // transform a value into a grid
#define expspace(i,xmin,xmax,echelle,n) ( echelle*(exp((log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1))*(i))+((xmin/echelle)-1.0)) ) // transform a grid into a value
#define invexpspace(x,xmin,xmax,echelle,n) ( log((x)/echelle - ((xmin/echelle)-1.0))/(log((xmax/echelle)-((xmin/echelle)-1.0))/(n-1)) ) // transform a value into a grid
// this spacing mix an exponential and then an equispaced grid //
#define fungrid(i,xmin1,xmax1,xmin2,xmax2,echelle,n1,n2) ((i)<(n1))?(echelle*(exp((log((xmax1/echelle)-((xmin1/echelle)-1.0))/(n1-1))*(i))+((xmin1/echelle)-1.0))):((i-n1+1)*(((xmax2)-(xmin2))/(n2))+(xmin2))
#define invfungrid(x,xmin1,xmax1,xmin2,xmax2,echelle,n1,n2) ((x)<=(xmax1))?(log((x)/echelle - ((xmin1/echelle)-1.0))/(log((xmax1/echelle)-((xmin1/echelle)-1.0))/(n1-1))):((((x)-(xmin2))/((xmax2)-(xmin2)))*(n2)+n1-1)




// QUADRATIC INTERPOLATION //
double interQuad1d(const double dx, const double h0, const double h1, const double h2) // dx lies in the interval [i-0.5; i+0.5] (distance is equal to 1), we fake f(i-0.5);
{
double h05, h15, a, b, c, fapprox;

h05 = (h1+h0)/2.0;
h15 = (h2+h1)/2.0;

c = h05;
b = 2.0*(h1-h05);
a = h15-2.0*h1+h05;

fapprox = a*dx*dx+b*dx+c;

return fapprox;

}



// TRILINEAR INTERPOLATION // (source: wikipedia)
double inter3d(const double dx, const double dy, const double dz, const double c000, const double c001, const double c010, const double c011, const double c100, const double c101, const double c110, const double c111)
{
    double c00, c01, c10, c11, c0, c1, c;
    
    c00 = c000*(1.0-dx)+c100*dx;
    c01 = c001*(1.0-dx)+c101*dx;
    c10 = c010*(1.0-dx)+c110*dx;
    c11 = c011*(1.0-dx)+c111*dx;
    
    c0 = c00*(1.0-dy)+c10*dy;
    c1 = c01*(1.0-dy)+c11*dy;
    
    c = c0*(1.0-dz)+c1*dz;
    
    return c;
}



/*// FIND THE WEIGHT FOR QUADRATIC INTERPOLATION //
double weightQuad(double x, double xgrid, int *ixgrid, int *ixgridm1, int *ixgridp1, double vector[], int dimweight)
{
    double dxgrid, vector05, vector15;
    int caseW;
    
    caseW = 0;
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    
    if((xgrid-*ixgrid) >= 0.5) { // also get rid of the case where xgrid > 0.5 when xgrid > dimweight.
        *ixgrid = *ixgrid + 1;
    }
    
    *ixgrid=max(0,*ixgrid);
    
    //printf("%f %f %d \n", x, xgrid, *ixgrid);getchar();
    
    if (*ixgrid>=(dimweight-1))
    {
        xgrid=(dimweight-1)-0.000000001;
        *ixgrid=dimweight-2;
        caseW = 1;
    }
    
    if (*ixgrid<=0)
    {
        xgrid=0.000000001;
        *ixgrid=0;
        caseW = 2;
    }
    
    if(x < 0.0){
        x = 0.0;
    }
    
    //printf("%f %f %d \n", x, xgrid, *ixgrid);getchar();
    
    // interior case //
    if(caseW == 0) {
        vector05 = phi(*ixgrid-0.5);
        vector15 = phi(*ixgrid+0.5);
        *ixgridm1 = *ixgrid-1;
        *ixgridp1 = *ixgrid+1;
        dxgrid=(x-vector05)/(vector15-vector05);
    }
    
    // max case //
    if(caseW == 1) { // do linear interpolation in that case. (sufficient to contrat that ixgrid = dimweight-2).
        *ixgridm1 = *ixgrid;
        *ixgridp1 = *ixgrid+1;
        dxgrid=(x-vector[*ixgrid])/(vector[*ixgrid+1]-vector[*ixgrid]);
    }
    
    // min case //
    if(caseW == 2) {  // do linear interpolation in that case (sufficient to control that ixgrid == 0)
        *ixgridm1 = *ixgrid;
        *ixgridp1 = *ixgrid+1;
        dxgrid=(x-vector[*ixgrid])/(vector[*ixgrid+1]-vector[*ixgrid]);
    }
    
    
    if(dxgrid < 0.0) {
        printf("%f %f %d %d %d %f %f \n", vector15, vector05, *ixgrid, *ixgridm1, *ixgridp1, x, dxgrid);getchar();
    }
    // if i is in [0.5;dimweight-2 + 0.5], this is ok
    // if i is in [-0.5; 0.5], then replace -0.5 by 0, and xgrid by 0.00001
    // if i is in [dimweight-2 + 0.5; dimweight-1+0.5], then replace dimweight-1+0.5 by dimweight and xgrid by (dimweight-1)-0.000001;
    
    
    return dxgrid;
}
*/



// FIND THE WEIGHT FOR INTERPOLATION //
double weightinter(double x, double xgrid, int *ixgrid, double vector[], int dimweight)
{
    double dxgrid;
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        xgrid=(dimweight-1)-0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        xgrid=0.000000001;
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    dxgrid=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
    
    return dxgrid;
    
}




// FOR QSORT
int comparefun2(const void* a, const void* b) 
{ 
	double* da = (double*)a; 
	double* db = (double*)b; 
	int diff1 = (da[0] > db[0]) - (da[0] < db[0]); 
	if (diff1 != 0) return diff1; 
	return (da[1] > db[1]) - (da[1] < db[1]); 
}




// SHIFTER //

namespace {
	inline void shft2(double &a, double &b, const double c)
	{
		a=b;
		b=c;
	}
    
	inline void shft3(double &a, double &b, double &c, const double d)
	{
		a=b;
		b=c;
		c=d;
	}
}

using namespace std;


void bascule(double *VectorIN, double *VectorOUT, int dim)
{
    int i;
    for(i=0; i<dim; i++){VectorOUT[i]=VectorIN[i];}
}



// GET RESIDUALS OF THE HISTOGRAM //
void weighthist(double x, double xgrid, double *residout, int *ixgrid, double vector[], int dimweight)
{
    
    // PUT THESE VALUE ON A GRID //
    *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
    *ixgrid=max(0,*ixgrid);
    
    if (*ixgrid>=(dimweight-1))
    {
        printf("getresid: (ixgrid>(dimweight-1))");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    if (*ixgrid<=0)
    {
        printf("getresid: (ixgrid<0)");getchar();
        *ixgrid=min((dimweight-1),(int)(floor(xgrid)));
        *ixgrid=max(0,*ixgrid);
    }
    
    *residout=(x-vector[*ixgrid])/(vector[(*ixgrid+1)]-vector[*ixgrid]);
}




// INVARIANT DISTRIBUTION //
template<size_t dim, size_t dim2>
void inv_distri(double (&invdist)[dim], double (&prob)[dim][dim2])
{
    double tempdist[dim], critdist, sumdist;
    int i, j;
    
    invdist[1] = 1.0;
    
    critdist = 1.0;
    while(critdist > 0.00000001) {
        
        for(i=0;i<dim;i++){
            tempdist[i] = invdist[i];
        }
        
        // compute the invdist //
        for(i = 0; i<dim; i++){
            invdist[i] = 0.0;
        }
        
        for(i = 0; i<dim; i++){
            for(j = 0; j<dim2; j++) {
                invdist[i] += tempdist[j]*prob[j][i];
            }
        }
        
        critdist = 0.0;
        for(i=0;i<dim;i++) {
            critdist = max(abs(invdist[i] - tempdist[i]), critdist);
        }
        
       //printf("%f %f %f %f %f %f %f", invdist[0], invdist[1], invdist[2], tempdist[0], tempdist[1], tempdist[2], critdist); getchar();
        
    }
    
    // renormalize invdist //
    sumdist = 0.0;
    for(i = 0; i<dim; i++) {
        sumdist += invdist[i];
    }
    for(i=0;i<dim;i++) {
        invdist[i] = invdist[i]/sumdist;
    }
    
}
