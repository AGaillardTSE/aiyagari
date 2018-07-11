
//This program implements the Tauchen (1986) procedure to approximate an AR(1) process with a N-states (node) first order process
//Assume AR(1) is yt=(1-rho)*amu+rho*yt-1+et
//Assume et-->is N(0,sigma*sigma)
//Following Tauchen (1986) the largest node will be placed at m*STDDEV(Yt)


// BY SUMUDU KANKANAMGE (2017) //


//--------Functions definition:

double CDFSTDNormal(double x)
{
	//Function CDFSTDNormal: computes the standard normal CDF using Abramowiz and Stegun (1964) approximation


	// constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
	
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
	
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
	
    return 0.5*(1.0 + sign*y);

	
}




template<size_t Nodes>
void tauchenfun(double rho, double m, double amu, double sigma, double grid[Nodes], double Ptransi[Nodes][Nodes])

{

	
	
	int i,j,k;
	
	
	
	//First lets compute unconditiontal variance of yt
	
	double varyt= (sigma*sigma)/(1-rho*rho);
	
	//Compute stddev of yt
	
	double stdyt= sqrt(varyt);
	
	//Define maximum and minimum grid point
	
	double ynodes[Nodes];
	ynodes[Nodes-1]=m*stdyt;
	ynodes[0]=-ynodes[Nodes-1];
	
	//Define interior nodes
	
	double ynodesinterval=(ynodes[Nodes-1]-ynodes[0])/((Nodes-1)*1.0);
	
	for (i=1; i<(Nodes-1); i++) {
		ynodes[i]=ynodes[i-1]+ynodesinterval;
	}
	
	for (i=0; i<Nodes; i++)
	{
		ynodes[i]=ynodes[i]+amu;
	}
	

	for (i=0; i<Nodes; i++) {
		grid[i]=ynodes[i];
	}
	
	
	//Computing transition probability matrix
	
	double transitionMat[Nodes][Nodes];
	
	
	for (j=0; j<Nodes; j++) 
	{
		for (k=1; k<(Nodes-1); k++) 
		{
			
			transitionMat[j][k]=CDFSTDNormal((ynodes[k]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma)-CDFSTDNormal((ynodes[k]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma);
			
		}
		
		transitionMat[j][0]=CDFSTDNormal((ynodes[0]-(1-rho)*amu-rho*ynodes[j]+ynodesinterval/2.0)/sigma);
		transitionMat[j][Nodes-1]=1.0-CDFSTDNormal((ynodes[Nodes-1]-(1-rho)*amu-rho*ynodes[j]-ynodesinterval/2.0)/sigma);
		
	}
	
	
	for (j=0; j<Nodes; j++) 
	{
		for (k=0; k<(Nodes); k++) 
		{
			Ptransi[j][k] = transitionMat[j][k];
		}
	}
    

}

