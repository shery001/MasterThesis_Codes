// AMDG
#include "pluto.h"
#include "Opacity.h"
#include "OpacitySemenov.h"
#include "Radiation.h"
#include "FLD.h"
#include "PhysicalConstantsCGS.h"
#include "DomainDecomposition.h"
#include "ReferenceSystem.h"
#include "MakemakeTools.h" // LinearInterpolation

#if STELLAREVOLUTION == YES
#include "StellarEvolution.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//
// For Bell & Lin (1994):
//    Taking only dust part since gas very wrong near 2000 K?
//
#define  BL94_ONLY_DUST  YES


// ************************
// ** External Parameter **
// ************************
//
int    DustOpacityFlag     = -100;
int    GasOpacityFlag      = -100;
int    OpacityExponent     = -100;
double ConstantDustOpacity = -100.0;
double ConstantGasOpacity  = -100.0;
int NSemenovDustModelTopShapeZK = 7;
char SemenovDustModelTopShape[7];

extern int IrradiationFlag;


// ************************
// ** Internal Parameter **
// ************************
//
int     freq_nr;
double *freq_nu;
double *opac_kabs;

double *meanopac_temp;
int     meanopac_ntemp;
double *meanopac_kapplanck;
double *meanopac_kapross;



// ***********************
// ** Internal Routines **
// ***********************
//
double RosselandMeanOpacity_Routine(double, double, int, int, int);

//
// To read tables in RADMC format and compute the Mean Opacities:
//
int ReadRADMCTablesAndComputeOpacities(FILE *, FILE *);
double bplanck(double, double);
double bplanckdt(double, double);
double integrate(int, double*, double*);

//
// Bell & Lin (1994) fitting routine:
//
double FittedOpacities(double, double);

//
// Gas opacities by Helling et al. (2000):
//
double **RosselandMeanGasOpacityArray;
double **PlanckMeanGasOpacityArray;
double *GasOpacityArrayDensityValues;
double *GasOpacityArrayTemperatureValues;
int nx = 71;
int ny = 71;
//int hunt(double*, int, double);
void InitializeGasOpacities(void);
void FinalizeGasOpacities(void);
//double RosselandMeanGasOpacity(double, double);

//
// Gas opacities by Malygin et al. (2012/2014):
//
double **NewRosselandMeanGasOpacityArray;
double **NewPlanckMeanGasOpacityArray;
double ***NewTwoTemPlanckMeanGasOpacityArray;
double *NewGasOpacityArrayDensityValues;
double *NewGasOpacityArrayTemperatureValues;
double *NewGasOpacityArrayRadiationTemperature;
int nRHO  = 99;
int nTEM  = -42; // initialise in InitializeOpacity() because it depends on the table (see flag UseTableDensityTo1em4)
int nTRAD = 10;
void NewInitializeGasOpacities();
void NewInitializeGasOpacitiesWithTTPL();
void NewFinalizeGasOpacities();
void NewFinalizeGasOpacitiesWithTTPL();
double NewPlanckMeanGasOpacity(double GasDensity, double RadiationTemperature);
double NewTwoTemPlanckMeanGasOpacity(double GasTemperature, double GasDensity, double RadiationTemperature);
double NewRosselandMeanGasOpacity(double GasDensity, double GasTemperature);
double GasOpacity_Malygin2012(double *, int, double *, int, double **, double, double);
double TwoTemPlanckGasOpacity_Malygin2013(double *, int, double *, int, double *, int, double ***, double, double, double);

// Use a table for the Malygin opacities that goes to higher rho? If not, "default" is used (version September 2013)
//
int UseTableDensityTo1em4 = 1;

// Allow linear extrapolation based on first or last two grid values? If not, take last value
//
int MalyginLinExtrapol = 0;


//
// macro operators to convert input into string:
//
#define STRINGIZE2(x) #x
#define STRINGIZE(x) STRINGIZE2(x)



// *************************
// ** InitializeOpacity() **
// *************************
//
int InitializeOpacity(){

	FILE *FrequencyFilePointer, *DustOpacityFilePointer;

#if LOG_OUTPUT != 0
	if(prank == 0) fprintf(LogFile, "###       ... Initialize Opacities ...                                       ###\n");
#else
	if(prank == 0) PetscPrintf(PETSC_COMM_WORLD, "###       ... Initialize Opacities ...                                       ###\n");
#endif

	char path2file[255];

#if LOG_OUTPUT != 0
	if(prank == 0) {
		fprintf(LogFile, "###           ... Using DustOpacityFlag = %-2d                                 ###\n", DustOpacityFlag);
		if(DustOpacityFlag==10)
			fprintf(LogFile, "###               ... Dust model (Semenov): %s                          ###\n", SemenovDustModelTopShape);
	}
#else
	if(prank == 0) {
		PetscPrintf(PETSC_COMM_WORLD, "###           ... Using DustOpacityFlag = %-2d                                 ###\n", DustOpacityFlag);
		if(DustOpacityFlag==10)
			PetscPrintf(PETSC_COMM_WORLD, "###               ... Dust model (Semenov): %s                          ###\n", SemenovDustModelTopShape);
	}
#endif
	switch(DustOpacityFlag){
		case 0:
		case 1:
			break;
		case 2:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"OssenkopfHenning1994/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"OssenkopfHenning1994/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 3:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"DraineLee1984/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"DraineLee1984/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 4:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"LaorDraine1993/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"LaorDraine1993/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 5:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"WeingartnerDraine2001/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"WeingartnerDraine2001/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 6:
#if IRRADIATION == YES
			// TODO Testen!
			printf(" Opacity.c: IrradiationFlag = %d\n",IrradiationFlag);
			if(IrradiationFlag != 0){
				if(prank == 0){
					printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
					printf("ERROR\n");
					printf("ERROR Bell & Lin (1994) opacities (DustOpacityFlag = 6) cannot be used\n");
					printf("ERROR for Irradiation (use IrradiationFlag = 0).\n");
					printf("ERROR\n");
					printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
				}
				QUIT_PLUTO(EXIT_FAILURE);
			}
#endif
			break;
		case 7:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"WeingartnerDraine2001_MC3D/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"WeingartnerDraine2001_MC3D/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 8:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"OssenkopfHenning1994_FUV_EUV/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"OssenkopfHenning1994_FUV_EUV/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		case 9:
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"LaorDraine1993_FUV_EUV/frequency.inp");
			FrequencyFilePointer =   fopen(path2file, "r");
			strcpy (path2file, STRINGIZE(BELT_DIR));
			strcat (path2file,"/src/Makemake2.1/Opacities/");
			strcat (path2file,"LaorDraine1993_FUV_EUV/dustopac_1.inp");
			DustOpacityFilePointer = fopen(path2file, "r");
			ReadRADMCTablesAndComputeOpacities(FrequencyFilePointer, DustOpacityFilePointer);
			break;
		//case 10:
			// Zeichenkette nach Fortran: siehe z.B. http://www.yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
			//
			//semenov03_staubig_init_(SemenovDustModelTopShape,&NSemenovDustModelTopShapeZK);
		//break;
		case 11:
			break;
		default:
			if(prank == 0){
				printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
				printf("ERROR\n");
				printf("ERROR DustOpacityFlag = %d not available, see 'Makemake.conf'.\n", DustOpacityFlag);
				printf("ERROR\n");
				printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
			}
			MPI_Barrier(MPI_COMM_WORLD);
			QUIT_PLUTO(EXIT_FAILURE);
	}


#if LOG_OUTPUT != 0
	if(prank == 0) fprintf(LogFile, "###           ... Using GasOpacityFlag = %-2d                                  ###\n", GasOpacityFlag);
#else
	if(prank == 0) PetscPrintf(PETSC_COMM_WORLD, "###           ... Using GasOpacityFlag = %-2d                                  ###\n", GasOpacityFlag);
#endif

	if (GasOpacityFlag==3 || GasOpacityFlag==4){
#if LOG_OUTPUT != 0
			if(prank == 0) fprintf(LogFile, "###               ... Linear Extrapolation from last two Values: %d           ###\n", MalyginLinExtrapol);
			if(prank == 0) fprintf(LogFile, "###               ... Opacity Table: using to rho=1e-4? %d                    ###\n", UseTableDensityTo1em4);
#else
			if(prank == 0) PetscPrintf(PETSC_COMM_WORLD, "###               ... Linear Extrapolation from last two Values: %d           ###\n", MalyginLinExtrapol);
			if(prank == 0) PetscPrintf(PETSC_COMM_WORLD, "###               ... Opacity Table: using to rho=1e-4? %d                    ###\n", UseTableDensityTo1em4);
#endif
	}

	switch(GasOpacityFlag){
		case 0:
		case 1:
			break;
		case 2:
			InitializeGasOpacities();
			break;
		case 3:
			nTEM = (UseTableDensityTo1em4==1)? 70 :92;
			NewInitializeGasOpacities();
			break;
		case 4:
			if(RadiationFlag==1){
				if (prank==0){
					printf("PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE\n");
					printf("PUZZLE \n");
					printf("PUZZLE Two-temperature opacities but with 1-T radiation transport...\n");
					printf("PUZZLE Strange but legal. Just sharing puzzlement before going on\n");
					printf("PUZZLE \n");
					printf("PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE PUZZLE\n");
				}
			}
			nTEM = (UseTableDensityTo1em4==1)? 70 :92;
			NewInitializeGasOpacitiesWithTTPL();
			break;
			case 5:
				break;
		default:
			if(prank == 0){
				printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
				printf("ERROR\n");
				printf("ERROR GasOpacityFlag = %d not available, see 'Makemake.conf'.\n", GasOpacityFlag);
				printf("ERROR\n");
				printf("ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR ERROR\n");
			}
			MPI_Barrier(MPI_COMM_WORLD);
			QUIT_PLUTO(EXIT_FAILURE);
	}

	return 0;
}




// ***********************
// ** FinalizeOpacity() **
// ***********************
//
// Frees the memory and destroys the objects allocated during InitializeRadiation().
//
int FinalizeOpacity(){

	switch(DustOpacityFlag){
		case 0:
		case 1:
			break;
		case 2:
		case 3:
		case 4:
		case 5:
		case 7:
		case 8:
		case 9:
			free(freq_nu);
			free(opac_kabs);
			free(meanopac_temp);
			free(meanopac_kapplanck);
			free(meanopac_kapross);
			break;
		case 10:
		case 11:
			break;

	}

	switch(GasOpacityFlag){
		case 0:
		case 1:
			break;
		case 2:
			FinalizeGasOpacities();
			break;
		case 3:
			NewFinalizeGasOpacities();
			break;
		case 4:
			NewFinalizeGasOpacitiesWithTTPL();
			break;
		case 5:
			break;
	}

	return 0;
}



// ******************************************
// ** ReadRADMCTablesAndComputeOpacities() **
// ******************************************
//
int ReadRADMCTablesAndComputeOpacities(FILE *FrequencyFilePointer, FILE *DustOpacityFilePointer){

	int i, itemp, inu;
	int nfdum;
	//	PetscErrorCode petscerrorcode;

	//
	// Read the frequency array from file:
	//
	double MaximumFrequency = 0.0;
	if(!(FrequencyFilePointer == NULL)){
		fscanf(FrequencyFilePointer, "%d", &freq_nr);

		//
		// Allocate memory for frequencies:
		//
		//		petscerrorcode = PetscNew(double[freq_nr], &freq_nu);
		freq_nu = malloc(freq_nr * sizeof(double));
		//		if(petscerrorcode){
		//			printf("ERROR: Unable to allocate memory in Opacity.c.\n");
		//			QUIT_PLUTO(EXIT_FAILURE);
		//		}

		//
		// Read frequencies:
		//
		for (i = 0; i < freq_nr; i++) {
			fscanf(FrequencyFilePointer, "%lE", &freq_nu[i]);
			MaximumFrequency = MAX(MaximumFrequency, freq_nu[i]);
		}
		//
		// Close file:
		//
		fclose (FrequencyFilePointer);
	}
	else{
		printf("ERROR: Unable to open the given frequency-input file.\n");
		QUIT_PLUTO(EXIT_FAILURE);
	}

	//
	// Read the dust opacities from file:
	//
	FILE *OpacityTableOut;
	char filename[256];
	int nomeaning;

	if(!(DustOpacityFilePointer == NULL)) {
		fscanf(DustOpacityFilePointer, "%d %d", &nfdum, &nomeaning);
		if (nfdum != freq_nr){
			printf("ERROR: frequency.inp and dustopac_1.inp have different number of frequencies/opacities.\n");
			QUIT_PLUTO(EXIT_FAILURE);
		}

		//
		// Allocate memory for frequency-dependent opacities:
		//
		//		petscerrorcode = PetscNew(double[freq_nr], &opac_kabs);
		opac_kabs = malloc(freq_nr * sizeof(double));
		//		if(petscerrorcode){
		//			printf("ERROR: Unable to allocate memory in Opacity.c.\n");
		//			QUIT_PLUTO(EXIT_FAILURE);
		//		}

		//
		// Read opacities:
		//
		sprintf(filename, "%s%s", "./data/", "FrequencyDependentOpacities.out");
		if(prank == 0) OpacityTableOut = fopen(filename, "w");
		for(i = 0; i < nfdum; i++) {
			fscanf (DustOpacityFilePointer, "%lE", &opac_kabs[i]);
			if(prank == 0) fprintf(OpacityTableOut, "%e %e\n", freq_nu[i], opac_kabs[i]);
		}
		//
		// Close file:
		//
		fclose (DustOpacityFilePointer);
	}
	else{
		printf("ERROR: Can't open 'dustopac_1.inp'\n");
		QUIT_PLUTO(EXIT_FAILURE);
	}
	if(prank == 0) fclose(OpacityTableOut);

	//
	// With gas-to-dust ratio convert from the dust opacity to an opacity for the gas+dust mixture:
	//
	//	for(i = 0; i < freq_nr; i++) {
	//		opac_kabs[i] = opac_kabs[i] / GasToDustRatio;
	//	}

	//
	// Make the temperature grid for the mean opacities:
	//
	//  meanopac_ntemp = size of mean opacity table
	//  tmin, tmax = lower and upper limit of temperature column of the mean opacity table
	// NO! Handled now as pure dust opacity
	//
	double tmin = 0.1; // in [K]
	double tmax = 100000.0; // in [K]
	meanopac_ntemp = 100000;

	//
	// Allocate memory:
	//
	//	petscerrorcode = PetscNew(double[meanopac_ntemp], &meanopac_temp);
	meanopac_temp = malloc(meanopac_ntemp * sizeof(double));
	//	if(petscerrorcode){
	//		printf("ERROR: Unable to allocate memory in Opacity.c.\n");
	//		QUIT_PLUTO(EXIT_FAILURE);
	//	}

	//
	// Calculate temperatures:
	//
	for(itemp = 0; itemp < meanopac_ntemp; itemp++) {
		meanopac_temp[itemp] = tmin * pow((tmax/tmin), itemp/(meanopac_ntemp-1.0));
	}

	//
	// Make the Planck opacities:
	//
	double *dumarr1, *dumarr2;
	sprintf(filename, "%s%s", "./data/", "PlanckMeanOpacities.out");
	if(prank == 0) OpacityTableOut = fopen(filename, "w");

	//
	// Allocate memory:
	//
	//	petscerrorcode = PetscNew(double[freq_nr], &dumarr1);
	dumarr1 = malloc(freq_nr * sizeof(double));
	//	if(petscerrorcode){
	//		printf("ERROR: Unable to allocate memory in Opacity.c.\n");
	//		QUIT_PLUTO(EXIT_FAILURE);
	//	}
	//	petscerrorcode = PetscNew(double[freq_nr], &dumarr2);
	dumarr2 = malloc(freq_nr * sizeof(double));
	//	if(petscerrorcode){
	//		printf("ERROR: Unable to allocate memory in Opacity.c.\n");
	//		QUIT_PLUTO(EXIT_FAILURE);
	//	}
	//	petscerrorcode = PetscNew(double[meanopac_ntemp], &meanopac_kapplanck);
	meanopac_kapplanck = malloc(meanopac_ntemp * sizeof(double));
	//	if(petscerrorcode){
	//		printf("ERROR: Unable to allocate memory in Opacity.c.\n");
	//		QUIT_PLUTO(EXIT_FAILURE);
	//	}

	//
	// 'Integrate' over the frequency-dependent Opacities for each point in the Temperature grid.
	//
	for(itemp = 0; itemp < meanopac_ntemp; itemp++) {
		for(inu = 0; inu < freq_nr; inu++) {
			dumarr1[inu] = bplanck(meanopac_temp[itemp], freq_nu[inu]);
			dumarr2[inu] = dumarr1[inu] * opac_kabs[inu];
		}
		double dum1 = integrate(freq_nr, freq_nu, dumarr1);
		double dum2 = integrate(freq_nr, freq_nu, dumarr2);
		meanopac_kapplanck[itemp] = dum2/dum1;
		if(prank == 0) fprintf(OpacityTableOut, "%e %e\n", meanopac_temp[itemp], meanopac_kapplanck[itemp]);
	}
	if(prank == 0) fclose(OpacityTableOut);


	//
	// Make the Rosseland opacities:
	//
	sprintf(filename, "%s%s", "./data/", "RosselandMeanOpacities.out");
	if(prank == 0) OpacityTableOut = fopen(filename, "w");

	//
	// Allocate memory:
	//
	//	petscerrorcode = PetscNew(double[meanopac_ntemp], &meanopac_kapross);
	meanopac_kapross = malloc(meanopac_ntemp * sizeof(double));
	//	if(petscerrorcode){
	//		printf("ERROR: Unable to allocate memory in Opacity.c.\n");
	//		QUIT_PLUTO(EXIT_FAILURE);
	//	}

	//
	// 'Integrate' over the frequency-dependent Opacities for each point in the Temperature grid.
	//
	for(itemp = 0; itemp < meanopac_ntemp; itemp++) {
		for(inu = 0; inu < freq_nr; inu++) {
			dumarr1[inu] = bplanckdt(meanopac_temp[itemp], freq_nu[inu]);
			dumarr2[inu] = dumarr1[inu] / opac_kabs[inu];
		}
		double dum1 = integrate(freq_nr, freq_nu, dumarr1);
		double dum2 = integrate(freq_nr, freq_nu, dumarr2);
		meanopac_kapross[itemp] = dum1/dum2;
		if(prank == 0) fprintf(OpacityTableOut, "%e %e\n", meanopac_temp[itemp], meanopac_kapross[itemp]);
	}
	if(prank == 0) fclose(OpacityTableOut);

	//
	// Free memory:
	//
	free(dumarr1);
	free(dumarr2);

	return 0;
}



// *****************************************************************
// *            THE BLACKBODY PLANCK FUNCTION B_nu(T)              *
// *                                                               *
// * This function computes the Blackbody function                 *
// *                                                               *
// *                2 h nu^3 / c^2                                 *
// *    B_nu(T)  = -------------------    [ erg / cm^2 s ster Hz ] *
// *               exp(-h nu / kT) - 1                             *
// *                                                               *
// *  ARGUMENTS:                                                   *
// *     nu    [Hz]            = Frequency                         *
// *     temp  [K]             = Electron temperature              *
// *****************************************************************
//
double bplanck (double temp, double nu) {
    if(temp != 0)
    	return 1.47455E-47 * pow(nu, 3) / (exp(4.7989E-11 * nu / temp)-1.) + 1.E-290;
    else
    	return 0;
}



// *************************************************************
// *       THE TEMPERATURE DERIVATIVE OF PLANCK FUNCTION       *
// *                                                           *
// *  This function computes the temperature derivative of the *
// *  Blackbody function                                       *
// *                                                           *
// *     dB_nu(T)     2 h^2 nu^4      exp(h nu / kT)        1  *
// *     --------   = ---------- ------------------------  --- *
// *        dT          k c^2    [ exp(h nu / kT) - 1 ]^2  T^2 *
// *                                                           *
// *   ARGUMENTS:                                              *
// *       nu    [Hz]            = Frequency                   *
// *       temp  [K]             = Electron temperature        *
// *************************************************************
//
double bplanckdt (double temp, double nu) {
    double theexp = exp(4.7989E-11*nu/temp);
    if (theexp < 1E+33)
    	return 7.07661334104E-58 * pow(nu, 4) * theexp / ( pow((theexp-1.), 2) * pow(temp, 2) ) + 1E-290;
    else
    	return 7.07661334104E-58 * pow(nu, 4) /( theexp * pow(temp, 2) ) + 1.E-290;
}



// ***********************
// * INTEGRATION ROUTINE *
// ***********************
//
double integrate(int n, double * x, double * f){

	int i;
    double integral = 0;

    for(i = 1; i < n; i++) {
    	integral +=
		0.5 * (f[i] + f[i-1])
		*
		(x[i] - x[i-1]);
    }
    if(x[n] > x[1])
    	return integral;
    else
    	return -integral;
}



// *********************************
// ** FrequencyDependentOpacity() **
// *********************************
//
// Returns the (linear interpolated) frequency dependent opacity for given input frequency [in code units].
//
double FrequencyDependentOpacity(double Frequency){

	double Opacity;
//  double frac;
	int index;

	//
	// Convert input frequency into cgs units [Hz]:
	//
	Frequency *= ReferenceFrequency;

	//
	// Determine position in frequency table:
	//
	index = hunt(freq_nu, freq_nr, Frequency);

	//
	// Return the 'extreme' opacities for 'extreme' frequencies:
	//
	if(index == -1)             Opacity = opac_kabs[0];
	else if(index == freq_nr-1) Opacity = opac_kabs[freq_nr-1];
	else{
		//
		// Linear interpolation of neighboring opacities:
		//
        if(freq_nu[index] < freq_nu[index+1])
            Opacity =
                LinearInterpolation1D(
                                      freq_nu[index]  , freq_nu[index+1],
                                      opac_kabs[index], opac_kabs[index+1],
                                      Frequency
                                      );
        else
            Opacity =
                LinearInterpolation1D(
                                      freq_nu[index+1]  , freq_nu[index],
                                      opac_kabs[index+1], opac_kabs[index],
                                      Frequency
                                      );
//		frac =
//		(Frequency        - freq_nu[index])
//		/
//		(freq_nu[index+1] - freq_nu[index]);
//
//		Opacity =
//		(1-frac) * opac_kabs[index]
//		+
//		frac     * opac_kabs[index+1];
	}

#if DEBUGGING == 1
	if(index == -1 || index == freq_nr-1){
		printf("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING\n");
		printf("WARNING\n");
		printf("WARNING The frequency input into 'FrequencyDependentOpacity' is\n");
		printf("WARNING out of the range of the specified frequency/opacity table.\n");
		printf("WARNING  Frequency input = %e Hz\n", Frequency);
		printf("WARNING  Frequency[0] =    %e Hz\n", freq_nu[0]);
		printf("WARNING  Frequency[imax] = %e Hz\n", freq_nu[freq_nr-1]);
		printf("WARNING\n");
		printf("WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING\n");
		QUIT_PLUTO(EXIT_FAILURE);
	}
#endif

	//
	// Convert Opacity to code units:
	//
	Opacity /= ReferenceOpacity;

	return Opacity;
}



// *************************
// ** PlanckMeanOpacity() **
// *************************
//
// TODO: Should be using RadiationTemperature or DustTemperature as argument to opacity tables?
//       -> Probably both
//
double PlanckMeanDustOpacity(double RadiationTemperature, double DustTemperature, double GasDensity){

	double PlanckMeanOpacity, frac;
	int index;
	int Typ = 2;  // für Semenov et al. (2003), nur Staub; 1 = Rosseland, 2 = Planck

	//
	// Convert input quantities into cgs units:
	//
	RadiationTemperature *= ReferenceTemperature;
	DustTemperature      *= ReferenceTemperature;
	GasDensity           *= ReferenceDensity;

	PlanckMeanOpacity = -1.0;

	switch(DustOpacityFlag){
		case 0:
			PlanckMeanOpacity = 0.0;
			break;
		case 1:
			PlanckMeanOpacity = ConstantDustOpacity;
			break;
		case 2:
		case 3:
		case 4:
		case 5:
		case 7:
		case 8:
		case 9:
			//
			// Determine position in meanopac_temp table:
			//
			index = hunt(meanopac_temp, meanopac_ntemp, RadiationTemperature);

			//
			// Return the 'extreme' opacities for 'extreme' temperatures:
			//
			if(index == -1)	                   PlanckMeanOpacity = meanopac_kapplanck[0];
			else if(index == meanopac_ntemp-1) PlanckMeanOpacity = meanopac_kapplanck[meanopac_ntemp-1];
			else{
				//
				// Linear interpolation of neighboring opacities:
				//
// TODO: Why does the routine call make a difference (e.g. in the "Pascucci_Tau=100.0_Full" Radiation Transport test)
//                PlanckMeanOpacity =
//                    LinearInterpolation1D(
//                                          meanopac_temp[index]     , meanopac_temp[index+1],
//                                          meanopac_kapplanck[index], meanopac_kapplanck[index+1],
//                                          RadiationTemperature
//                                          );

//                PlanckMeanOpacity = (meanopac_kapplanck[index+1] - meanopac_kapplanck[index]) / (meanopac_temp[index+1] - meanopac_temp[index]) * (RadiationTemperature - meanopac_temp[index]) + meanopac_kapplanck[index];

				frac =
				(RadiationTemperature   - meanopac_temp[index])
				/
				(meanopac_temp[index+1] - meanopac_temp[index]);

#if DEBUGGING == 1
				if( (frac < 0) || (frac > 1) ){
					printf("ERROR during calculation of Planck Mean Opacity: frac < 0 or frac > 1\n");
					printf("trad =                %e\n", RadiationTemperature);
					printf("meanopac_temp[%4d] = %e\n", index, meanopac_temp[index]);
					printf("meanopac_temp[%4d] = %e\n", index+1, meanopac_temp[index+1]);
					printf("frac =                  %e\n", frac);
					QUIT_PLUTO(EXIT_FAILURE);
				}
#endif

				PlanckMeanOpacity =
				(1-frac) * meanopac_kapplanck[index]
				+
				frac     * meanopac_kapplanck[index+1];
			}
			break;
		case 6:
			// Bell & Lin (1994) Planck mean opacities not available
			//   WARNING We use Rosseland mean here too!
			//   WARNING This is wrong especially because of the density dependence (very weak for Planck)
			// TODO: Use Trad? Sollte wohl… Aber dann Evaporation nur nach TStaub. Aber hier egal da extrem falsch
//#warning "###########################################################"
//#warning "##  Using Rosseland mean as Planck for Bell & Lin 1994!  ##"
//#warning "###########################################################"
			PlanckMeanOpacity = FittedOpacities(DustTemperature, GasDensity);
			break;
		//case 10:
			// TODO: Use Trad? Sollte wohl… Aber dann Evaporation nur nach TStaub -> Semenov ändern
			//semenov03_staubig_(&Typ, &GasDensity, &DustTemperature, &PlanckMeanOpacity);
			//break;
		default:
			PetscFPrintf(PETSC_COMM_WORLD, LogFile, "### ERROR: DustOpacityFlag = %d is not in allowed range                      ###\n", DustOpacityFlag);
			QUIT_PLUTO(EXIT_FAILURE);
			break;
		case 11:
			PlanckMeanOpacity = 0.02*pow(DustTemperature/10.0, 2);
			break;
	}

	//
	// Convert PlanckMeanOpacity into code units:
	//
	PlanckMeanOpacity /= ReferenceOpacity;

	return PlanckMeanOpacity;
}



// ****************************
// ** RosselandMeanOpacity() **
// ****************************
//
// TODO: Should be using RadiationTemperature or DustTemperature as argument to opacity tables?
//
double RosselandMeanDustOpacity(double DustTemperature, double GasDensity){

	double RosselandMeanOpacity;
	int index;
	int Typ = 2;  // für Semenov et al. (2003), nur Staub; 1 = Rosseland, 2 = Planck

	//
	// Convert input quantities from code units to cgs units:
	//
	GasDensity     *= ReferenceDensity;
	DustTemperature *= ReferenceTemperature;

	RosselandMeanOpacity = -1.0;

	switch(DustOpacityFlag){
		case 0:
			RosselandMeanOpacity = 0.0;
			break;
		case 1:
			RosselandMeanOpacity = ConstantDustOpacity;
			break;
		case 2:
		case 3:
		case 4:
		case 5:
		case 7:
		case 8:
		case 9:
			//
			// Determine position in meanopac_temp table:
			//
			index = hunt(meanopac_temp, meanopac_ntemp, DustTemperature);
			//
			// Return the 'extreme' opacities for 'extreme' temperatures:
			//
			if(index == -1)                    RosselandMeanOpacity = meanopac_kapross[0];
			else if(index == meanopac_ntemp-1) RosselandMeanOpacity = meanopac_kapross[meanopac_ntemp-1];
			else{
				//
				// Linear interpolation of neighbouring opacities:
				//
				RosselandMeanOpacity =
							LinearInterpolation1D(
									  meanopac_temp[index]   , meanopac_temp[index+1],
									  meanopac_kapross[index], meanopac_kapross[index+1],
									  DustTemperature
									  );
                //				frac =
                //	    		(DustTemperature        - meanopac_temp[index])
                //	    		/
                //	    		(meanopac_temp[index+1] - meanopac_temp[index]);
                //#if DEBUG > 0
                //				if( (frac < 0) || (frac > 1) ){
                //					printf("ERROR during calculation of Rosseland Mean Opacity: frac < 0 or frac > 1\n");
                //					printf("DustTemperature =      %e\n", DustTemperature);
                //					printf("meanopac_temp[%4d] = %e\n", index,   meanopac_temp[index]);
                //					printf("meanopac_temp[%4d] = %e\n", index+1, meanopac_temp[index+1]);
                //					printf("frac =                %e\n", frac);
                //					QUIT_PLUTO(EXIT_FAILURE);
                //				}
                //#endif
                //				RosselandMeanOpacity =
                //	        	(1-frac) * meanopac_kapross[index]
                //	        	+
                //	        	frac     * meanopac_kapross[index+1];
			}
			break;
		case 6:
			RosselandMeanOpacity = FittedOpacities(DustTemperature, GasDensity);
			break;
		//case 10:
		//	semenov03_staubig_(&Typ, &GasDensity, &DustTemperature, &RosselandMeanOpacity);
		//	break;
		default:
			PetscFPrintf(PETSC_COMM_WORLD, LogFile, "### ERROR: DustOpacityFlag = %d is not in allowed range                      ###\n", DustOpacityFlag);
			QUIT_PLUTO(EXIT_FAILURE);
			break;
		case 11:
			RosselandMeanOpacity = 0.02*pow(DustTemperature/10.0, 2);
			break;
	}

	//
	// Convert RosselandMeanOpacity output from cgs units to code units:
	//
	RosselandMeanOpacity /= ReferenceOpacity;

	return RosselandMeanOpacity;
}



// *********************
// ** FittedOpacities **
// *********************
//
// Get opacity from fitted values. See Bell & Lin (1994).
//
double FittedOpacities(double DustTemperature, double GasDensity){
//
//	The opacity table is divided into eight regions
//	region 1 is dominated by ice grains
//	region 2 is where ice grains melts
//	region 3 is dominated by metal grains
//	region 4 is where metal grains melts
//	region 5 is dominated by molecules
//	region 6 is dominated by H-
//	region 7 is dominated by Kramer's law
//	region 8 is dominated by electron scattering
//
	double k;
	double RosselandMeanOpacity;
	double temp, rho;
	double num1 = 1.e22, num2 = 6.5e-5, num3;
	double k1,k2,k3,k4,k5,k6,k7,k8;
	double k11,k22,k33,k44,k66,k77,k88;

	// Opacity coefficients ki
	// Cannot handle 2.e81 in c[3] --> 2.e21
#if BL94_ONLY_DUST == YES
	double c[8] = {2.e-4, 2.e16, 0.1, 2.e21,   0.0,    0.0,    0.0,   0.0};
#else
	double c[8] = {2.e-4, 2.e16, 0.1, 2.e21, 1.e-8, 1.e-36, 1.5e20, 0.348};
#endif

	// Density exponent a
	double a[8] = {0.0, 0.0, 0.0, 1.0, 2/3., 1/3., 1.0, 0.0};

	// Temperature exponent b
	double b[8] = {2.0, -7.0, 1./2., -24.0, 3.0, 10.0, -5/2., 0.0};

	// Different regions..
	double power[3] = {2.8369e-2, 1.1464e-2, 2.2667e-1};
	double coef[3] = {1.46e3, 4.51e3, 2.37e6};

	//
	// Convert input quantities into cgs units:
	//
	GasDensity      *= ReferenceDensity;
	DustTemperature *= ReferenceTemperature;

	temp = DustTemperature;
	rho = GasDensity;

	if(temp < coef[0]*(pow(rho,power[0])) ){
		k1 = c[0]*pow(rho,a[0])*pow(temp,b[0]);
		k2 = c[1]*pow(rho,a[1])*pow(temp,b[1]);
		k3 = c[2]*pow(rho,a[2])*pow(temp,b[2]);
		k11 = pow(k1,2.);
		k22 = pow(k2,2.);
		k = pow( pow((k11*k22/(k11+k22)),2.) + pow(k3/(1.+(num1/pow(temp,10.))),4) ,1/4.);
		RosselandMeanOpacity = k;
	}
	else {
		if (temp < (coef[1]*pow(rho,power[1])) ){
			k3 = c[2]*pow(rho,a[2])*pow(temp,b[2]);
			k4 = c[3]*pow(rho,a[3])*pow(temp,b[3])*(1.e20)*(1.e20)*(1.e20);
			k5 = c[4]*pow(rho,a[4])*pow(temp,b[4]);
			k33 = pow(k3,4);
			k44 = pow(k4,4);
			num3 = pow(1 - 4*temp,8);
			k = pow( (k33*k44/(k33+k44)) + pow(k5/(1.+(num2/(num3*pow(rho,2/3.)*100.))),4) ,1/4.);
			RosselandMeanOpacity = k;
		}
		else {
			if ( temp < (coef[2]*pow(rho,power[2])) ){
				k5 = c[4]*pow(rho,a[4])*pow(temp,b[4]);
				k6 = c[5]*pow(rho,a[5])*pow(temp,b[5]);
				k7 = c[6]*pow(rho,a[6])*pow(temp,b[6]);
				k66 = pow(k6,2);
				k77 = pow(k7,2);
				num3 = (1-4*temp);
				k = pow( pow(k66*k77/(k66+k77),2) + pow(k5/(1.+pow(num3/(1.1*pow(rho,0.04762)),10)),4) ,1/4.);
				RosselandMeanOpacity = k;
			}
			else {
				k7 = c[6]*pow(rho,a[6])*pow(temp,b[6]);
				k8 = c[7]*pow(rho,a[7])*pow(temp,b[7]);
				k77 = pow(k7,2);
				k88 = pow(k8,2);
				k = pow(k77 + k88,1/2.);
				RosselandMeanOpacity = k;
			}
		}
	}

	return RosselandMeanOpacity;
}



// *****************************************************************************
// ** Gas opacities by Helling et al. (2000) (as used by Semenov et al. 2003) **
// *****************************************************************************
//
void InitializeGasOpacities(){

	FILE *GasOpacityFile;
	int i, ix, iy;

	//
	// Allocate memory for the gas opacity density and temperature grid values:
	//
	GasOpacityArrayDensityValues =     malloc(nx * sizeof(double));
	GasOpacityArrayTemperatureValues = malloc(ny * sizeof(double));

	//
	// Set the values of the density - temperature grid:
	//
	GasOpacityArrayDensityValues[0]  = 0.2364E-06;
	GasOpacityArrayDensityValues[1]  = 0.1646E-06;
	GasOpacityArrayDensityValues[2]  = 0.1147E-06;
	GasOpacityArrayDensityValues[3]  = 0.7985E-07;
	GasOpacityArrayDensityValues[4]  = 0.5560E-07;
	GasOpacityArrayDensityValues[5]  = 0.3872E-07;
	GasOpacityArrayDensityValues[6]  = 0.2697E-07;
	GasOpacityArrayDensityValues[7]  = 0.1878E-07;
	GasOpacityArrayDensityValues[8]  = 0.1308E-07;
	GasOpacityArrayDensityValues[9]  = 0.9107E-08;
	GasOpacityArrayDensityValues[10] = 0.6342E-08;
	GasOpacityArrayDensityValues[11] = 0.4417E-08;
	GasOpacityArrayDensityValues[12] = 0.3076E-08;
	GasOpacityArrayDensityValues[13] = 0.2142E-08;
	GasOpacityArrayDensityValues[14] = 0.1492E-08;
	GasOpacityArrayDensityValues[15] = 0.1039E-08;
	GasOpacityArrayDensityValues[16] = 0.7234E-09;
	GasOpacityArrayDensityValues[17] = 0.5038E-09;
	GasOpacityArrayDensityValues[18] = 0.3508E-09;
	GasOpacityArrayDensityValues[19] = 0.2443E-09;
	GasOpacityArrayDensityValues[20] = 0.1701E-09;
	GasOpacityArrayDensityValues[21] = 0.1185E-09;
	GasOpacityArrayDensityValues[22] = 0.8252E-10;
	GasOpacityArrayDensityValues[23] = 0.5746E-10;
	GasOpacityArrayDensityValues[24] = 0.4002E-10;
	GasOpacityArrayDensityValues[25] = 0.2787E-10;
	GasOpacityArrayDensityValues[26] = 0.1941E-10;
	GasOpacityArrayDensityValues[27] = 0.1352E-10;
	GasOpacityArrayDensityValues[28] = 0.9412E-11;
	GasOpacityArrayDensityValues[29] = 0.6554E-11;
	GasOpacityArrayDensityValues[30] = 0.4565E-11;
	GasOpacityArrayDensityValues[31] = 0.3179E-11;
	GasOpacityArrayDensityValues[32] = 0.2214E-11;
	GasOpacityArrayDensityValues[33] = 0.1542E-11;
	GasOpacityArrayDensityValues[34] = 0.1074E-11;
	GasOpacityArrayDensityValues[35] = 0.7476E-12;
	GasOpacityArrayDensityValues[36] = 0.5206E-12;
	GasOpacityArrayDensityValues[37] = 0.3626E-12;
	GasOpacityArrayDensityValues[38] = 0.2525E-12;
	GasOpacityArrayDensityValues[39] = 0.1758E-12;
	GasOpacityArrayDensityValues[40] = 0.1225E-12;
	GasOpacityArrayDensityValues[41] = 0.8528E-13;
	GasOpacityArrayDensityValues[42] = 0.5939E-13;
	GasOpacityArrayDensityValues[43] = 0.4136E-13;
	GasOpacityArrayDensityValues[44] = 0.2880E-13;
	GasOpacityArrayDensityValues[45] = 0.2006E-13;
	GasOpacityArrayDensityValues[46] = 0.1397E-13;
	GasOpacityArrayDensityValues[47] = 0.9727E-14;
	GasOpacityArrayDensityValues[48] = 0.6774E-14;
	GasOpacityArrayDensityValues[49] = 0.4717E-14;
	GasOpacityArrayDensityValues[50] = 0.3285E-14;
	GasOpacityArrayDensityValues[51] = 0.2288E-14;
	GasOpacityArrayDensityValues[52] = 0.1593E-14;
	GasOpacityArrayDensityValues[53] = 0.1109E-14;
	GasOpacityArrayDensityValues[54] = 0.7726E-15;
	GasOpacityArrayDensityValues[55] = 0.5381E-15;
	GasOpacityArrayDensityValues[56] = 0.3747E-15;
	GasOpacityArrayDensityValues[57] = 0.2609E-15;
	GasOpacityArrayDensityValues[58] = 0.1817E-15;
	GasOpacityArrayDensityValues[59] = 0.1265E-15;
	GasOpacityArrayDensityValues[60] = 0.8813E-16;
	GasOpacityArrayDensityValues[61] = 0.6137E-16;
	GasOpacityArrayDensityValues[62] = 0.4274E-16;
	GasOpacityArrayDensityValues[63] = 0.2976E-16;
	GasOpacityArrayDensityValues[64] = 0.2073E-16;
	GasOpacityArrayDensityValues[65] = 0.1443E-16;
	GasOpacityArrayDensityValues[66] = 0.1005E-16;
	GasOpacityArrayDensityValues[67] = 0.7000E-17;
	GasOpacityArrayDensityValues[68] = 0.4875E-17;
	GasOpacityArrayDensityValues[69] = 0.3395E-17;
	GasOpacityArrayDensityValues[70] = 0.2364E-17;

	GasOpacityArrayTemperatureValues[0]  =   500.00;
	GasOpacityArrayTemperatureValues[1]  =   521.86;
	GasOpacityArrayTemperatureValues[2]  =   544.68;
	GasOpacityArrayTemperatureValues[3]  =   568.50;
	GasOpacityArrayTemperatureValues[4]  =   593.35;
	GasOpacityArrayTemperatureValues[5]  =   619.30;
	GasOpacityArrayTemperatureValues[6]  =   646.38;
	GasOpacityArrayTemperatureValues[7]  =   674.64;
	GasOpacityArrayTemperatureValues[8]  =   704.14;
	GasOpacityArrayTemperatureValues[9]  =   734.93;
	GasOpacityArrayTemperatureValues[10] =   767.06;
	GasOpacityArrayTemperatureValues[11] =   800.60;
	GasOpacityArrayTemperatureValues[12] =   835.61;
	GasOpacityArrayTemperatureValues[13] =   872.15;
	GasOpacityArrayTemperatureValues[14] =   910.28;
	GasOpacityArrayTemperatureValues[15] =   950.08;
	GasOpacityArrayTemperatureValues[16] =   991.63;
	GasOpacityArrayTemperatureValues[17] =  1034.99;
	GasOpacityArrayTemperatureValues[18] =  1080.24;
	GasOpacityArrayTemperatureValues[19] =  1127.47;
	GasOpacityArrayTemperatureValues[20] =  1176.77;
	GasOpacityArrayTemperatureValues[21] =  1228.23;
	GasOpacityArrayTemperatureValues[22] =  1281.93;
	GasOpacityArrayTemperatureValues[23] =  1337.99;
	GasOpacityArrayTemperatureValues[24] =  1396.49;
	GasOpacityArrayTemperatureValues[25] =  1457.55;
	GasOpacityArrayTemperatureValues[26] =  1521.28;
	GasOpacityArrayTemperatureValues[27] =  1587.80;
	GasOpacityArrayTemperatureValues[28] =  1657.23;
	GasOpacityArrayTemperatureValues[29] =  1729.69;
	GasOpacityArrayTemperatureValues[30] =  1805.32;
	GasOpacityArrayTemperatureValues[31] =  1884.26;
	GasOpacityArrayTemperatureValues[32] =  1966.65;
	GasOpacityArrayTemperatureValues[33] =  2052.64;
	GasOpacityArrayTemperatureValues[34] =  2142.39;
	GasOpacityArrayTemperatureValues[35] =  2236.07;
	GasOpacityArrayTemperatureValues[36] =  2333.84;
	GasOpacityArrayTemperatureValues[37] =  2435.89;
	GasOpacityArrayTemperatureValues[38] =  2542.40;
	GasOpacityArrayTemperatureValues[39] =  2653.56;
	GasOpacityArrayTemperatureValues[40] =  2769.59;
	GasOpacityArrayTemperatureValues[41] =  2890.69;
	GasOpacityArrayTemperatureValues[42] =  3017.09;
	GasOpacityArrayTemperatureValues[43] =  3149.01;
	GasOpacityArrayTemperatureValues[44] =  3286.70;
	GasOpacityArrayTemperatureValues[45] =  3430.41;
	GasOpacityArrayTemperatureValues[46] =  3580.41;
	GasOpacityArrayTemperatureValues[47] =  3736.96;
	GasOpacityArrayTemperatureValues[48] =  3900.36;
	GasOpacityArrayTemperatureValues[49] =  4070.91;
	GasOpacityArrayTemperatureValues[50] =  4248.91;
	GasOpacityArrayTemperatureValues[51] =  4434.69;
	GasOpacityArrayTemperatureValues[52] =  4628.60;
	GasOpacityArrayTemperatureValues[53] =  4830.98;
	GasOpacityArrayTemperatureValues[54] =  5042.22;
	GasOpacityArrayTemperatureValues[55] =  5262.69;
	GasOpacityArrayTemperatureValues[56] =  5492.80;
	GasOpacityArrayTemperatureValues[57] =  5732.98;
	GasOpacityArrayTemperatureValues[58] =  5983.65;
	GasOpacityArrayTemperatureValues[59] =  6245.29;
	GasOpacityArrayTemperatureValues[60] =  6518.36;
	GasOpacityArrayTemperatureValues[61] =  6803.38;
	GasOpacityArrayTemperatureValues[62] =  7100.86;
	GasOpacityArrayTemperatureValues[63] =  7411.34;
	GasOpacityArrayTemperatureValues[64] =  7735.41;
	GasOpacityArrayTemperatureValues[65] =  8073.64;
	GasOpacityArrayTemperatureValues[66] =  8426.66;
	GasOpacityArrayTemperatureValues[67] =  8795.12;
	GasOpacityArrayTemperatureValues[68] =  9179.68;
	GasOpacityArrayTemperatureValues[69] =  9581.07;
	GasOpacityArrayTemperatureValues[70] = 10000.00;


	//
	// Allocate memory for the gas opacity array:
	//
	RosselandMeanGasOpacityArray = malloc(ny*sizeof(double));
	for (i = 0; i < ny; i++) {
		RosselandMeanGasOpacityArray[i] = malloc(nx*sizeof(double));
	}

	//
	// Read the gas opacity file into the gas opacity array:
	//
	char path2file[255];
	strcpy(path2file, STRINGIZE(BELT_DIR));
	strcat(path2file,"/src/Makemake2.1/Opacities/");
	strcat(path2file,"Hellingetal2000/kR_h2001.dat");
	GasOpacityFile = fopen(path2file, "r");

	// Read the bla bla value:
	fscanf(GasOpacityFile, "%le", &RosselandMeanGasOpacityArray[0][0]);

	// Read the 2D array:
	for(ix = 0; ix < nx; ix++){
		for(iy = 0; iy < ny; iy++){
			fscanf(GasOpacityFile, "%le", &RosselandMeanGasOpacityArray[iy][ix]);
			//printf("%d %d: %e\n", iy, ix, RosselandMeanGasOpacityArray[iy][ix]);
		}
	}
	fclose(GasOpacityFile);

	//
	// Allocate memory for the gas opacity array:
	//
	PlanckMeanGasOpacityArray = malloc(ny*sizeof(double));
	for (i = 0; i < ny; i++) {
		PlanckMeanGasOpacityArray[i] = malloc(nx*sizeof(double));
	}

	//
	// Read the gas opacity file into the gas opacity array:
	//
	strcpy(path2file, STRINGIZE(BELT_DIR));
	strcat(path2file,"/src/Makemake2.1/Opacities/");
	strcat(path2file,"Hellingetal2000/kP_h2001.dat");
	GasOpacityFile = fopen(path2file, "r");

	// Read the bla bla value:
	fscanf(GasOpacityFile, "%le", &PlanckMeanGasOpacityArray[0][0]);

	// Read the 2D array:
	for(ix = 0; ix < nx; ix++){
		for(iy = 0; iy < ny; iy++){
			fscanf(GasOpacityFile, "%le", &PlanckMeanGasOpacityArray[iy][ix]);
			//printf("%d %d: %e\n", iy, ix, PlanckMeanGasOpacityArray[iy][ix]);
		}
	}
	fclose(GasOpacityFile);
}



void FinalizeGasOpacities(){
	free(RosselandMeanGasOpacityArray);
	free(PlanckMeanGasOpacityArray);
	free(GasOpacityArrayDensityValues);
	free(GasOpacityArrayTemperatureValues);
}



double GasOpacity_Hellingetal2000(double *GasDensities,    int NGasDensities,
								double *GasTemperatures, int NGasTemperatures,
								double **D,
								double GasDensity, double GasTemperature){

	int nlx, nly;

	//
	// Search the nearest lower grid point 'nly' to the 'GasTemperature':
	//
	nly = hunt(GasTemperatures, NGasTemperatures, GasTemperature);

	//
	// Search the nearest lower grid point 'nlx' to the 'GasDensity':
	//
	nlx = hunt(GasDensities, NGasDensities, GasDensity);

	//
	// Check, if GasTemperature and GasDensity are inside the given grid range:
	//
	if(GasDensity <= GasDensities[NGasDensities-1]){
		// use the opacity regarding the lowest gas density:
		nlx = NGasDensities-2;
		GasDensity = GasDensities[NGasDensities-1];
	}
	if(GasDensity > GasDensities[0]){
		//printf("ERROR: Gas opacity for densities > %e g cm^-3 not available yet.\n",GasDensities[0]);
		//QUIT_PLUTO(EXIT_FAILURE);
//#warning "*************************************************"
//#warning "***  konst. kappa fuer H00 bei rho > 2.36e-7  ***"
//#warning "*************************************************"
		// use the opacity regarding the highest gas density:
		nlx = 0;
		GasDensity = GasDensities[0];
	}
	if(GasTemperature < GasTemperatures[0]){
		// use the opacity regarding the lowest gas temperature:
		nly = 0;
		GasTemperature = GasTemperatures[0];
	}
	if(GasTemperature >= GasTemperatures[NGasTemperatures-1]){
		//printf("ERROR: Gas opacity for temperatures > 10000 K not available yet.\n");
		//QUIT_PLUTO(EXIT_FAILURE);
		// use the opacity regarding the highest gas temperature:
		nly = NGasTemperatures-2;
		GasTemperature = GasTemperatures[NGasTemperatures-1];
	}

	//
	// Interpolate:
	//
	return LinearInterpolation2D(GasDensities[nlx+1],  GasDensities[nlx],
								 GasTemperatures[nly], GasTemperatures[nly+1],
								 D[nly][nlx+1], D[nly][nlx], D[nly+1][nlx+1], D[nly+1][nlx],
								 GasDensity, GasTemperature);
}



double PlanckMeanGasOpacity(double RadiationTemperature, double GasTemperature, double GasDensity){

	double PlanckMeanOpacity;

	//
	// Convert to cgs:
	//
	RadiationTemperature *= ReferenceTemperature;
	GasTemperature       *= ReferenceTemperature;
	GasDensity           *= ReferenceDensity;

	PlanckMeanOpacity = -1.0;

	switch(GasOpacityFlag){
		case 0:
			PlanckMeanOpacity = 0.0;
			break;
		case 1:
			PlanckMeanOpacity = ConstantGasOpacity;
			break;
		case 2:
			PlanckMeanOpacity = GasOpacity_Hellingetal2000(GasOpacityArrayDensityValues, nx, GasOpacityArrayTemperatureValues, ny,
								       PlanckMeanGasOpacityArray, GasDensity, GasTemperature);

			//
			// Convert the opacity of the power exponent to the linear scale:
			//
			PlanckMeanOpacity = pow(10.0, PlanckMeanOpacity);
			break;
		case 3:
			// TODO wenn man 2T-Strahlungstransport macht aber 1T-Opazitäten rechnet,
			//      wäre es besser, RadiationTemperature zu benutzen!
			// TODO: Auch für Staub?
			//
			PlanckMeanOpacity = NewPlanckMeanGasOpacity(GasDensity, GasTemperature);
			break;
		case 4:
			PlanckMeanOpacity = NewTwoTemPlanckMeanGasOpacity(GasTemperature, GasDensity, RadiationTemperature);
			break;
		case 5:
			PlanckMeanOpacity = 0.02*pow(GasTemperature/10.0, 2);
			break;
		default:
			PetscFPrintf(PETSC_COMM_WORLD, LogFile, "### ERROR: GasOpacityFlag = %d is not in allowed range              ###\n", GasOpacityFlag);
			QUIT_PLUTO(EXIT_FAILURE);
	}

	//
	// Convert to code units:
	//
	PlanckMeanOpacity /= ReferenceOpacity;

	return PlanckMeanOpacity;
}



double RosselandMeanGasOpacity(double GasTemperature, double GasDensity){

	double RosselandMeanOpacity;

	//
	// Convert to cgs:
	//
	GasDensity     *= ReferenceDensity;
	GasTemperature *= ReferenceTemperature;

	RosselandMeanOpacity = -1.0;

	switch(GasOpacityFlag){
		case 0:
		RosselandMeanOpacity = 0.0;
			break;
		case 1:
			RosselandMeanOpacity = ConstantGasOpacity;
			break;
		case 2:
			RosselandMeanOpacity = GasOpacity_Hellingetal2000(GasOpacityArrayDensityValues, nx, GasOpacityArrayTemperatureValues, ny, RosselandMeanGasOpacityArray, GasDensity, GasTemperature);

			//
			// Convert the opacity of the power exponent to the linear scale:
			//
			RosselandMeanOpacity = pow(10.0, RosselandMeanOpacity);

			break;
		case 3:
		case 4:
			RosselandMeanOpacity = NewRosselandMeanGasOpacity(GasDensity, GasTemperature);
			break;
		case 5:
			RosselandMeanOpacity = 0.02*pow(GasTemperature/10.0, 2);
			break;
		default:
			PetscFPrintf(PETSC_COMM_WORLD, LogFile, "### ERROR: GasOpacityFlag = %d is not in allowed range              ###\n", GasOpacityFlag);
			QUIT_PLUTO(EXIT_FAILURE);
	}

	//
	// Convert to code units:
	//
	RosselandMeanOpacity /= ReferenceOpacity;

	return RosselandMeanOpacity;
}



//
// **********************************************************************
// ** Gas opacities by Malygin et al. (2013) Two-Temperature Planck    **
// **********************************************************************
//
void NewInitializeGasOpacitiesWithTTPL(){

	FILE *GasOpacityFile;
	char path2file[255];
	int iRAD, iRHO, iTEM;

	//
	// Allocate memory for the gas opacity density and temperature grid values:
	//
	if(nRHO<0 || nTEM <0 || nTRAD <0){
		printf("[NewInitializeGasOpacitiesWithTTPL()] Problem! nRHO = %d oder nTEM = %d oder nTRAD = %d <0\n",nRHO,nTEM,nTRAD);
		QUIT_PLUTO(EXIT_FAILURE);
	}
	NewGasOpacityArrayDensityValues        = malloc(nRHO  * sizeof(double));
	NewGasOpacityArrayTemperatureValues    = malloc(nTEM  * sizeof(double));
	NewGasOpacityArrayRadiationTemperature = malloc(nTRAD * sizeof(double));

	//
	// Read the gas opacity density values:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.rhovalues.dat");
	else
		strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.rhovalues.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 1D array:
	for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayDensityValues[iRHO]);
	}
	fclose(GasOpacityFile);

	//
	// Read the gas opacity temperature values:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.Tvalues.dat");
	else
		strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.Tvalues.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 1D array:
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
		fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayTemperatureValues[iTEM]);
	}
	fclose(GasOpacityFile);


	//
	// Allocate memory for the Rosseland gas opacity array:
	//
	NewRosselandMeanGasOpacityArray = malloc (nTEM*sizeof(double *)) ;
	for (iTEM = 0; iTEM < nTEM; iTEM++) {
		NewRosselandMeanGasOpacityArray[iTEM] = malloc (nRHO*sizeof(double)) ;
	}

	//
	// Read the gas opacity file into the Rosseland gas opacity array:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat (path2file,"Malyginetal2014/kR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.dat");
	else
		strcat (path2file,"Malyginetal2014/kR_92tx99rho_2013.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	//
	// Read the 2D array:
	//
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
	    for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewRosselandMeanGasOpacityArray[iTEM][iRHO]);
	    }
	}
	fclose(GasOpacityFile);

	//
	// Allocate memory for the 2-D Single-Temperature Planck Mean Opacity array:
	//
	NewPlanckMeanGasOpacityArray = malloc(nTEM*sizeof(double *));
	for (iTEM = 0; iTEM < nTEM; iTEM++) {
		NewPlanckMeanGasOpacityArray[iTEM] = malloc(nRHO*sizeof(double));
	}

	//
	// Read the 2-D Single-Temperature Planck gas opacity file into the Planck gas opacity array:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kP_p00_T700-21544_rho1e-04_1e-18_70tx99rho.dat");
	else
		strcat(path2file, "Malyginetal2014/kP_92tx99rho_2013.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// read the 2D array:
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
	    for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewPlanckMeanGasOpacityArray[iTEM][iRHO]);
	    }
	}
	fclose(GasOpacityFile);

	//
	// Read the 2-T Planck gas opacity radiation temperature values:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1) {
		// TODO missing: ttpl extended table
		printf(" *** Extended table (to rho=1e-4) for two-T is not yet available! Change gas flag ***\n");
		QUIT_PLUTO(EXIT_FAILURE);
	}
	else
		strcat(path2file, "Malyginetal2014/ttpl_10tradx99rhox92tgas_2013.Tradvalues.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 1D array:
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
		fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayRadiationTemperature[iTEM]);
	}
	fclose(GasOpacityFile);

	//
	// Allocate memory for the 3D gas two-temperature Planck opacity array:
	//
	NewTwoTemPlanckMeanGasOpacityArray = malloc(nTRAD*sizeof(double **));
	for (iRAD = 0; iRAD < nTRAD; iRAD++) {
	    NewTwoTemPlanckMeanGasOpacityArray[iRAD] = malloc(nRHO*sizeof(double *));
	    for (iRHO = 0; iRHO < nRHO; iRHO++) {
		NewTwoTemPlanckMeanGasOpacityArray[iRAD][iRHO] = malloc(nTEM*sizeof(double));
	    }
	}

	//
	// Read the gas opacity file into the gas opacity array:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	// TODO missing: ttpl extended table
	if (UseTableDensityTo1em4==1){
		printf(" *** Extended table (to rho=1e-4) for two-T is not yet available! Change gas flag ***\n");
		QUIT_PLUTO(EXIT_FAILURE);
	}
	else
		strcat(path2file, "Malyginetal2014/ttpl_10tradx99rhox92tgas_2013.dat");

	GasOpacityFile = fopen(path2file, "r");
	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacitiesWithTTPL()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 3D array:
	for(iRAD = 0; iRAD < nTRAD; ++iRAD){
	    for(iRHO = 0; iRHO < nRHO; ++iRHO){
		for(iTEM = 0; iTEM < nTEM; ++iTEM){
		    fscanf(GasOpacityFile, "%lf", &NewTwoTemPlanckMeanGasOpacityArray[iRAD][iRHO][iTEM]);
		}
	    }
	}
	fclose(GasOpacityFile);
}

void NewFinalizeGasOpacitiesWithTTPL(){
	free(NewRosselandMeanGasOpacityArray);
	free(NewTwoTemPlanckMeanGasOpacityArray);
	free(NewPlanckMeanGasOpacityArray);
	free(NewGasOpacityArrayDensityValues);
	free(NewGasOpacityArrayTemperatureValues);
}


//
// **********************************************************************
// ** Gas opacities by Malygin et al. (2012) Single-TemperaturePlanck **
// **********************************************************************
//
void NewInitializeGasOpacities(){

	FILE *GasOpacityFile;
	char path2file[255];
	int i, iRHO, iTEM;

	//
	// Allocate memory for the gas opacity density and temperature grid values:
	//
	if(nRHO<0 || nTEM <0){
		printf("[NewInitializeGasOpacities()] Problem! nRHO = %d oder nTEM = %d <0\n",nRHO,nTEM);
		QUIT_PLUTO(EXIT_FAILURE);
	}
	NewGasOpacityArrayDensityValues     = malloc(nRHO * sizeof(double));
	NewGasOpacityArrayTemperatureValues = malloc(nTEM * sizeof(double));

	//
	// Read the gas opacity density values:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.rhovalues.dat");
	else
		strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.rhovalues.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacities()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 1D array:
	for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayDensityValues[iRHO]);
	}
	fclose(GasOpacityFile);

	//
	// Read the gas opacity temperature values:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kPR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.Tvalues.dat");
	else
		strcat(path2file, "Malyginetal2014/kPR_92tx99rho_2013.Tvalues.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacities()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 1D array:
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
		fscanf(GasOpacityFile, "%lf", &NewGasOpacityArrayTemperatureValues[iTEM]);
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
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kR_p00_T700-21544_rho1e-04_1e-18_70tx99rho.dat");
	else
		strcat(path2file, "Malyginetal2014/kR_92tx99rho_2013.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacities()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	// Read the 2D array:
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
	    for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewRosselandMeanGasOpacityArray[iTEM][iRHO]);
		//printf("%d %d: %e\n", iTEM, iRHO, NewRosselandMeanGasOpacityArray[iTEM][iRHO]);
	    }
	}
	fclose(GasOpacityFile);

	//
	// Allocate memory for the gas Planck opacity array:
	//
	NewPlanckMeanGasOpacityArray = malloc(nTEM*sizeof(double *));  // (double *) vorne entfernt (siehe oben)
	for (i = 0; i < nTEM; i++) {
		NewPlanckMeanGasOpacityArray[i] = malloc(nRHO*sizeof(double));
	}

	//
	// Read the gas Planck opacity file into the gas opacity array:
	//
	strcpy (path2file, STRINGIZE(BELT_DIR));
	strcat (path2file,"/src/Makemake2.1/Opacities/");
	if (UseTableDensityTo1em4==1)
		strcat(path2file, "Malyginetal2014/kP_p00_T700-21544_rho1e-04_1e-18_70tx99rho.dat");
	else
		strcat(path2file, "Malyginetal2014/kP_92tx99rho_2013.dat");

	GasOpacityFile = fopen(path2file, "r");

	if (GasOpacityFile == NULL) {
	    printf("[NewInitializeGasOpacities()] Problem mit `%s'\n", path2file);
	    QUIT_PLUTO(EXIT_FAILURE);
	}

	//
	// Read the 2D array:
	//
	for(iTEM = 0; iTEM < nTEM; ++iTEM){
	    for(iRHO = 0; iRHO < nRHO; ++iRHO){
		fscanf(GasOpacityFile, "%lf", &NewPlanckMeanGasOpacityArray[iTEM][iRHO]);
	    }
	}
	fclose(GasOpacityFile);
}

void NewFinalizeGasOpacities(){
	free(NewRosselandMeanGasOpacityArray);
	free(NewPlanckMeanGasOpacityArray);
	free(NewGasOpacityArrayDensityValues);
	free(NewGasOpacityArrayTemperatureValues);
}


double TwoTemPlanckGasOpacity_Malygin2013(double *GasDensities,    int NGasDensities,
			      double *GasTemperatures, int NGasTemperatures,
			      double *RadiationTemperatures, int NRadiationTemperatures,
			      double ***D,
			      double GasDensity, double GasTemperature, double RadiationTemperature){
	//
	// D - is a 3D array with two-temperature Planck opacities
	//
	int nlrho, nltem, nlrad;
	//
	// Search the nearest lower grid point 'nltem' to the 'GasTemperature':
	//
	nltem = hunt(GasTemperatures, NGasTemperatures, GasTemperature);
	//
	// Search the nearest lower grid point 'nlrho' to the 'GasDensity':
	//
	nlrho = hunt(GasDensities, NGasDensities, GasDensity);
	//
	// Search the nearest lower grid point 'nlrad' to the 'RadiationTemperature':
	//
	nlrad = hunt(RadiationTemperatures, NRadiationTemperatures, RadiationTemperature);
	//
	// Check, if GasTemperature GasDensity and RadiationTemperature are inside the given grid range:
	//
	if(GasDensity <= GasDensities[nRHO-1]){
		nlrho = nRHO-2;
		// use the opacity of the lowest gas density or not:
		if(MalyginLinExtrapol==0) GasDensity = GasDensities[nRHO-1];
	}
	if(GasDensity > GasDensities[0]){
		//printf("ERROR: KAPPA GAS OPACITY: Gas opacity for densities > %0.4E g cm^-3 are not available yet.\n", GasDensities[0]);
		//QUIT_PLUTO(EXIT_FAILURE);
		nlrho = 0;
		// use the opacity of the highest gas density or not:
		if(MalyginLinExtrapol==0) GasDensity = GasDensities[0];
	}
	if(GasTemperature < GasTemperatures[0]){
		nltem = 0;
		// use the opacity of the lowest gas temperature or not:
		if(MalyginLinExtrapol==0) GasTemperature = GasTemperatures[0];
	}
	if(GasTemperature >= GasTemperatures[nTEM-1]){
		//printf("ERROR: KAPPA GAS OPACITY: Gas opacity for temperature %f > %f K is not available yet.\n", GasTemperature, GasTemperatures[nTEM-1]);
		//QUIT_PLUTO(EXIT_FAILURE);
		nltem = nTEM-2;
		// use the opacity of the highest gas temperature or not:
		if(MalyginLinExtrapol==0) GasTemperature = GasTemperatures[nTEM-1];
	}
	if(RadiationTemperature < RadiationTemperatures[0]){
		nlrad = 0;
		// use the opacity of the lowest radiation temperature or not:
		if(MalyginLinExtrapol==0) RadiationTemperature = RadiationTemperatures[0];
	}
	if(RadiationTemperature > RadiationTemperatures[nTRAD-1]){
		//printf("ERROR: KAPPA GAS OPACITY: Gas opacity for radiation temperature %f > %f K is not available yet.\n", RadiationTemperature, RadiationTemperatures[nTRAD-1]);
		//QUIT_PLUTO(EXIT_FAILURE);
		nlrad = nTRAD-2;
		// use the opacity of the highest gas temperature or not:
		if(MalyginLinExtrapol==0) RadiationTemperature = RadiationTemperatures[nTRAD-1];
	}
	//
	// Interpolate:
	//
	return LinearInterpolation3D(RadiationTemperatures[nlrad], RadiationTemperatures[nlrad+1],
				GasDensities[nlrho+1] ,  GasDensities[nlrho]    ,
				GasTemperatures[nltem], GasTemperatures[nltem+1],
				D[nlrad][nlrho+1][nltem]  , D[nlrad+1][nlrho+1][nltem]  ,
				D[nlrad][nlrho]  [nltem]  , D[nlrad+1][nlrho]  [nltem]  ,
				D[nlrad][nlrho+1][nltem+1], D[nlrad+1][nlrho+1][nltem+1],
				D[nlrad][nlrho]  [nltem+1], D[nlrad+1][nlrho]  [nltem+1],
				RadiationTemperature, GasDensity, GasTemperature);
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


double NewRosselandMeanGasOpacity(double GasDensity, double GasTemperature){

	double RosselandMeanOpacity;

	RosselandMeanOpacity =
	GasOpacity_Malygin2012(NewGasOpacityArrayDensityValues, nRHO,
			       NewGasOpacityArrayTemperatureValues, nTEM,
			       NewRosselandMeanGasOpacityArray,
			       GasDensity, GasTemperature);

	//
	// Convert the opacity of the power exponent to the linear scale:
	//
	RosselandMeanOpacity = pow(10.0, RosselandMeanOpacity);

	return RosselandMeanOpacity;
}

double NewPlanckMeanGasOpacity(double GasDensity, double RadiationTemperature){

	double PlanckMeanOpacity;

	PlanckMeanOpacity =
	GasOpacity_Malygin2012(NewGasOpacityArrayDensityValues, nRHO,
			       NewGasOpacityArrayTemperatureValues, nTEM,
			       NewPlanckMeanGasOpacityArray,
			       GasDensity, RadiationTemperature);

	//
	// Convert the opacity of the power exponent to the linear scale:
	//
	PlanckMeanOpacity = pow(10.0, PlanckMeanOpacity);

	return PlanckMeanOpacity;
}

double NewTwoTemPlanckMeanGasOpacity(double GasTemperature, double GasDensity, double RadiationTemperature){

	double PlanckMeanOpacity;

	// TODO ZUTUN: Sollten wir nicht immer (aus numerischen und vielleicht praktischen Gründen)
	//             TwoTemPlanckGasOpacity_Malygin2013() aufrufen, auch wenn T = Trad?
	//             Prüfen, aber, daß die Interpolation gut genug ist.
	//             Oder nur dann die 1-T-Opazitäten benutzen, wenn Trad nicht im Gitter ist?
	// TODO Aufpassen, wenn das (rho,T)-Gitter für beide Opazitäten nicht gleich ist!
	//             Wäre ein Problem für breite Tabelle
	if (GasTemperature != RadiationTemperature) {
	    PlanckMeanOpacity =
		TwoTemPlanckGasOpacity_Malygin2013(NewGasOpacityArrayDensityValues, nRHO,
			NewGasOpacityArrayTemperatureValues, nTEM,
			NewGasOpacityArrayRadiationTemperature, nTRAD,
			NewTwoTemPlanckMeanGasOpacityArray,
			GasDensity, GasTemperature, RadiationTemperature);
	} else {
	    PlanckMeanOpacity = GasOpacity_Malygin2012(NewGasOpacityArrayDensityValues, nRHO,
			NewGasOpacityArrayTemperatureValues, nTEM,
			NewPlanckMeanGasOpacityArray,
			GasDensity, RadiationTemperature);
	}

	//
	// Convert the opacity of the power exponent to the linear scale:
	//
	PlanckMeanOpacity = pow(10.0, PlanckMeanOpacity);

	return PlanckMeanOpacity;
}
