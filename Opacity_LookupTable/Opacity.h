#ifndef OPACITY_H_
#define OPACITY_H_

extern int DustOpacityFlag;
//extern char SemenovDustModelTopShape[7];
extern int OpacityExponent;
extern double ConstantDustOpacity;

extern int GasOpacityFlag;
extern double ConstantGasOpacity;

//
// Frequency dependent Opacities:
//
// TODO: should be not required; currently used in Irradiation.c
extern int freq_nr;
extern double *freq_nu;


// ***********************
// ** External Routines **
// ***********************
//
int InitializeOpacity(void);
int FinalizeOpacity(void);

double FrequencyDependentOpacity(double);
// Achtung: The dust opacity data of Semenov et al. (2003) and Bell & Lin (1994),
//          the only ones currently available, require the *gas* density
//
double PlanckMeanDustOpacity(double RadiationTemperature, double DustTemperature, double GasDensity);
double RosselandMeanDustOpacity(double DustTemperature, double GasDensity);
double PlanckMeanGasOpacity(double RadiationTemperature, double GasTemperature, double GasDensity);
double RosselandMeanGasOpacity(double GasTemperature, double GasDensity);

#endif /*OPACITY_H_*/
