#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <math.h>

#include "Opacity.h"

double Fluxlimiter(double GasDensity, double GasTemperature, double gradt);
int main(){

double rho = 5.67E-08;
double tem = 5000.0;
double gradt = 1.0;

double lambda = Fluxlimiter(rho, tem, gradt);

printf("fluxlimiter = %lf \n", lambda);

//FinalizeGasOpacities();

 return 0;
}


double Fluxlimiter(double GasDensity, double GasTemperature, double gradt)
{
  double k_op, sigma, lambda, R;
  k_op = RosselandMeanGasOpacity(GasDensity,GasTemperature);
  sigma = 0.0;
  R = 4*gradt /(GasDensity*(k_op + sigma)*GasTemperature);
  if(R >=0 && R<=2)
    {
      lambda = 2/(3 + pow(9+10*R,0.5));
    }
  else if(R>2)
    {
      lambda = 10/(10*R + 9 + pow(180*R + 81, 0.5)) ;
     }
      else
	{ 
     printf("error in valuse of R");
	 exit(1);
	}
      return lambda;
    
}
