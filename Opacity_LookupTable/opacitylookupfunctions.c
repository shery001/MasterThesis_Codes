#include <stdio.h>
#include <string.h>
#include <stdlib.h> 
#include <math.h>


/*......................declaration.................*/

double **NewRosselandMeanGasOpacityArray;
double *NewGasOpacityArrayDensityValues;
double *NewGasOpacityArrayTemperatureValues;
double GasOpacity_Malygin2012(double *, int, double *, int, double **, double, double);
double RosselandMeanGasOpacity(double , double );
int hunt(double *, int , double );
double LinearInterpolation1D( double , double , double , double , double );
double LinearInterpolation2D(double , double , double , double , double , double , double , double , double , double );
void InitializeGasOpacities();
void FinalizeGasOpacities();
int nRHO  = 99;
int nTEM  = 92; // initialise in InitializeOpacity() because it depends on the table (see flag UseTableDensityTo1em4)
int nTRAD = 1000;

int MalyginLinExtrapol=0;
int DustOpacityFlag = 0;
int UseTableDensityTo1em4 = 1;


/*.........*/

int main(){

double rho = 5.67E-08;
double tem = 5000.0;

double opacity = RosselandMeanGasOpacity(rho, tem);

printf("opacity = %lf \n", opacity);

FinalizeGasOpacities();

 return 0;
}

double RosselandMeanGasOpacity(double GasDensity, double GasTemperature){
  /*......................declaration.................*/
  //double **NewRosselandMeanGasOpacityArray;
  //double *NewGasOpacityArrayDensityValues;
  //double *NewGasOpacityArrayTemperatureValues;
  //int nRHO  = 99;
  //int nTEM  = 42; // initialise in InitializeOpacity() because it depends on the table (see flag UseTableDensityTo1em4)
  //int nTRAD = 10;

  /*.........*/

  InitializeGasOpacities();

  double RosselandMeanOpacity;

  RosselandMeanOpacity =
    GasOpacity_Malygin2012(NewGasOpacityArrayDensityValues, nRHO,
			   NewGasOpacityArrayTemperatureValues, nTEM,
			   NewRosselandMeanGasOpacityArray,
			   GasDensity, GasTemperature);

  //
  // Convert the opacity of the power exponent to the linear scale:
  //
  //RosselandMeanOpacity = pow(10.0, RosselandMeanOpacity);

  return RosselandMeanOpacity;
  
}


double GasOpacity_Malygin2012(double *GasDensities,    int NGasDensities,
			      double *GasTemperatures, int NGasTemperatures,
			      double **D,
			      double GasDensity, double GasTemperature){
  static unsigned long int Nrm=0, Nrp=0, Ntm=0, Ntp=0, Nrlog=20000000;  // Max: 4'294'967'295. 2e4 war zu oft für die kleine Tabelle
  //
  // D - is a 2D array with opacities (Rosseland or Planck)
  //
  int nlrho, nltem;
  int nRHO = NGasDensities;
  int nTEM = NGasTemperatures;
  //
  // Search the nearest lower grid point 'nltem' to the 'GasTemperature':
  //
  nltem = hunt(GasTemperatures, NGasTemperatures, GasTemperature);
  //
  // Search the nearest lower grid point 'nlrho' to the 'GasDensity':
  //
  nlrho = hunt(GasDensities, NGasDensities, GasDensity);
  //printf("nlrho = %d, nltem = %d \n", nlrho, nltem);
  //
  // Check, if GasTemperature and GasDensity are inside the given grid range:
  //
  if(GasDensity < GasDensities[nRHO-1]){
    if (Nrm%Nrlog == 0)
      printf("[Fuer %.0e Male stellvertretend] WARNING: KAPPA GAS OPACITY: Gas opacity for density %0.4E < %0.4E g cm^-3 are not available yet.\n",
	     (double) Nrlog-1., GasDensity, GasDensities[nRHO-1]);
    nlrho = nRHO-2;  // this was =0 before 17.06.2014

    // use the opacity of the lowest gas density or not:
    //
    if(MalyginLinExtrapol==0) {
      if (Nrm%Nrlog == 0) printf("         Den Wert bei rho_min = %e benutzen\n",GasDensities[nRHO-1]);
      GasDensity = GasDensities[nRHO-1];
    } else
      if (Nrm%Nrlog == 0) printf("         Es wird extrapoliert (denn MalyginLinExtrapol = %d)\n",MalyginLinExtrapol);
    Nrm++;
  }
  if(GasDensity > GasDensities[0]){
    if (Nrp%Nrlog == 0)
      printf("[Fuer %.0e Male stellvertretend] WARNING: KAPPA GAS OPACITY: Gas opacity for density %0.4E > %0.4E g cm^-3 are not available yet.\n",
	     (double) Nrlog-1., GasDensity, GasDensities[0]);
    nlrho = 0;       // this was =nRHO-2 before 17.06.2014

    // use the opacity of the highest gas density or not:
    //
    if(MalyginLinExtrapol==0) {
      if (Nrp%Nrlog == 0)
	printf("         Den Wert bei rho_max = %e benutzen\n",GasDensities[0]);
      GasDensity = GasDensities[0];
    } else
      if (Nrp%Nrlog == 0)
	printf("         Es wird extrapoliert (denn MalyginLinExtrapol = %d)\n",MalyginLinExtrapol);
    Nrp++;
  }
  if(GasTemperature < GasTemperatures[0]){
    //
    // normalerweise ist die Staubopazität höher als die Gasopazität -> kein Problem wenn extrapolieren muß
    //
    if (DustOpacityFlag == 0)
      if (Ntm%Nrlog == 0)
	printf("[Fuer %.0e Male stellvertretend] WARNING: KAPPA GAS OPACITY: Gas opacity for temperature %f < %f cgs are not available yet.\n",
	       (double) Nrlog-1., GasTemperature, GasTemperatures[0]);
    nltem = 0;

    // use the opacity of the lowest gas temperature or not:
    if(MalyginLinExtrapol==0) {
      if (DustOpacityFlag == 0)
	if (Ntm%Nrlog == 0)
	  printf("         Den Wert bei T_min = %e benutzen\n",GasTemperatures[0]);
      GasTemperature = GasTemperatures[0];
    } else
      if (DustOpacityFlag == 0)
	if (Ntm%Nrlog == 0)
	  printf("         Es wird extrapoliert (denn MalyginLinExtrapol = %d)\n",MalyginLinExtrapol);
    Ntm++;
  }
  if(GasTemperature > GasTemperatures[nTEM-1]){
    if (Ntp%Nrlog == 0)
      printf("[Fuer %.0e Male stellvertretend] WARNING: KAPPA GAS OPACITY: Gas opacity for temperature %f > %f K is not available yet.\n",
	     (double) Nrlog-1., GasTemperature, GasTemperatures[nTEM-1]);
    nltem = nTEM-2;

    // use the opacity of the highest gas temperature or not:
    if(MalyginLinExtrapol==0) {
      if (Ntp%Nrlog == 0)
	printf("         Den Wert bei T_max = %e benutzen\n",GasTemperatures[nTEM-1]);
      GasTemperature = GasTemperatures[nTEM-1];
    } else
      if (Ntp%Nrlog == 0)
	printf("         Es wird extrapoliert (denn MalyginLinExtrapol = %d)\n",MalyginLinExtrapol);
    Ntp++;
  }
  //
  // Interpolate:
  //
  return LinearInterpolation2D(GasDensities[nlrho+1] ,  GasDensities[nlrho]    ,
			       GasTemperatures[nltem], GasTemperatures[nltem+1],
			       D[nltem][nlrho+1], D[nltem][nlrho], D[nltem+1][nlrho+1], D[nltem+1][nlrho],
			       GasDensity, GasTemperature);
}


int hunt(double * xx, int n, double x){
  
  int jlo = 1;
  
  int jm, jhi, inc;
  int ascnd;
  
  //int n = xx.size();
  ascnd = (xx[n-1] >= xx[0]);
  if (jlo < 0 || jlo > n-1) {
    jlo = -1;
    jhi = n;
  }
  else {
    inc = 1;
    if ((x >= xx[jlo]) == ascnd) {
      if (jlo == n-1) return jlo;
      jhi = jlo+1;
      while ((x >= xx[jhi]) == ascnd) {
	jlo = jhi;
	inc += inc;
	jhi = jlo+inc;
	if (jhi > n-1) {
	  jhi = n;
	  break;
	}
      }
    }
    else {
      if (jlo == 0) {
	jlo = -1;
	return jlo;
      }
      jhi = jlo--;
      while ((x < xx[jlo]) == ascnd) {
	jhi = jlo;
	inc <<= 1;
	if (inc >= jhi) {
	  jlo = -1;
	  break;
	}
	else
	  jlo = jhi-inc;
      }
    }
  }
  while (jhi-jlo != 1) {
    jm = (jhi+jlo) >> 1;
    if ((x >= xx[jm]) == ascnd)
      jlo = jm;
    else
      jhi = jm;
  }
  if (x == xx[n-1])
    jlo = n-2;
  if (x == xx[0])
    jlo = 0;
  
  
  return jlo;
}


double LinearInterpolation1D(
			     double x0, double x1, 
			     double value_at_x0, double value_at_x1, 
			      double x
			     ){
  
  double gradient, value_at_x;
  
  //
  // Check, if x1 is greater than x0:
  //
  if(x0 >= x1){
    printf("ERROR: x0 >= x1\n");
    printf("x0 = %e\n", x0);
    printf("x1 = %e\n", x1);
    printf("x  = %e\n", x);
    exit(1);
  }
  
  //
  // Check, if x is inside the given grid range:
  //
  if(x < x0){
    printf("ERROR: x < x0\n");
    printf("x0 = %e\n", x0);
    printf("x1 = %e\n", x1);
    printf("x  = %e\n", x);
    exit(1);
  }
  if(x > x1){
    printf("ERROR: x > x1\n");
    printf("x0 = %e\n", x0);
    printf("x1 = %e\n", x1);
    printf("x  = %e\n", x);
    exit(1);
  }
  
  gradient = (value_at_x1 - value_at_x0) / (x1 - x0);
  value_at_x = gradient * (x - x0) + value_at_x0;
  
  return value_at_x;
}


double LinearInterpolation2D(
			     double x0, double x1, 
			     double y0, double y1, 
			     double value_at_x0y0, double value_at_x1y0, 
			     double value_at_x0y1, double value_at_x1y1, 
			     double x, double y
			     ){
  
  double value_at_xy;
  double value_at_xy0, value_at_xy1;
  
  //
  // 1D linear interpolation in x:
  //
  value_at_xy0 = LinearInterpolation1D(x0, x1, value_at_x0y0, value_at_x1y0, x);
  value_at_xy1 = LinearInterpolation1D(x0, x1, value_at_x0y1, value_at_x1y1, x);
  
  //
  // 1D linear interpolation in y:
  //
  value_at_xy = LinearInterpolation1D(y0, y1, value_at_xy0, value_at_xy1, y);
  
  return value_at_xy;
}

void InitializeGasOpacities(){

  FILE *GasOpacityFile;
  char path2file[255];
  int i, iRHO, iTEM;

  //
  // Allocate memory for the gas opacity density and temperature grid values:
  //
  if(nRHO<0 || nTEM <0){
    printf("[NewInitializeGasOpacities()] Problem! nRHO = %d oder nTEM = %d <0\n",nRHO,nTEM);
    exit(1);
  }
  NewGasOpacityArrayDensityValues     = malloc(nRHO * sizeof(double));
  NewGasOpacityArrayTemperatureValues = malloc(nTEM * sizeof(double));

  //
  // Read the gas opacity density values:
  //
  //strcpy (path2file, STRINGIZE(BELT_DIR));
  //strcat (path2file,"/src/Makemake2.1/Opacities/");
  strcpy (path2file, "/home/bhavya/Documents/AA_Master_thesis/Opacity/");
  //if (UseTableDensityTo1em4==1)
    //strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.rhovalues.dat");
  //else
    strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.rhovalues.dat");

  GasOpacityFile = fopen(path2file, "r");

  if (GasOpacityFile == NULL) {
    printf("[InitializeGasOpacities()] Problem with `%s'\n", path2file);
    exit(1);
  }

  // Read the 1D array:
  for(iRHO = 0; iRHO < nRHO; ++iRHO){
    fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayDensityValues[iRHO]);
    printf("%0.22lf \n", NewGasOpacityArrayDensityValues[iRHO]);
  }
  fclose(GasOpacityFile);

  //
  // Read the gas opacity temperature values:
  //
  //strcpy (path2file, STRINGIZE(BELT_DIR));
  //strcat (path2file,"/src/Makemake2.1/Opacities/");
  strcpy (path2file, "/home/bhavya/Documents/AA_Master_thesis/Opacity/");
  //if (UseTableDensityTo1em4==1)
    //strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.Tvalues.dat");
  //else
    strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.Tvalues.dat");

  GasOpacityFile = fopen(path2file, "r");

  if (GasOpacityFile == NULL) {
    printf("[InitializeGasOpacities()] Problem with `%s'\n", path2file);
    exit(1);
  }

  // Read the 1D array:
  for(iTEM = 0; iTEM < nTEM; ++iTEM){
    fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayTemperatureValues[iTEM]);
    //    printf("%lf \n", NewGasOpacityArrayTemperatureValues[iTEM]);
  }
  fclose(GasOpacityFile);


  //
  // Allocate memory for the gas Rosseland opacity array:
  //
  NewRosselandMeanGasOpacityArray = malloc (nTEM*sizeof(double *)) ;
  for (i = 0; i < nTEM; i++) {
    NewRosselandMeanGasOpacityArray[i] = malloc (nRHO*sizeof(double)) ;
  }

  //
  // Read the gas Rosseland opacity file into the gas opacity array:
  //
  //strcpy (path2file, STRINGIZE(BELT_DIR));
  //strcat (path2file,"/src/Makemake2.1/Opacities/");
  strcpy (path2file, "/home/bhavya/Documents/AA_Master_thesis/Opacity/");
  //if (UseTableDensityTo1em4==1)
  //strcat(path2file, "Malyginetal2014/kR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.dat");
  //else
    strcat(path2file, "Malyginetal2014/kR_92tx99rho_2013.dat");

  GasOpacityFile = fopen(path2file, "r");

  if (GasOpacityFile == NULL) {
    printf("[InitializeGasOpacities()] Problem with `%s'\n", path2file);
    exit(1);
  }

  // Read the 2D array:
  for(iTEM = 0; iTEM < nTEM; ++iTEM){
    for(iRHO = 0; iRHO < nRHO; ++iRHO){
      fscanf(GasOpacityFile, "%lf", &NewRosselandMeanGasOpacityArray[iTEM][iRHO]);
      //printf("%d %d: %e\n", iTEM, iRHO, NewRosselandMeanGasOpacityArray[iTEM][iRHO]);
    }
  }
  fclose(GasOpacityFile);

}

void FinalizeGasOpacities(){
	free(NewRosselandMeanGasOpacityArray);
	free(NewGasOpacityArrayDensityValues);
	free(NewGasOpacityArrayTemperatureValues);
}
