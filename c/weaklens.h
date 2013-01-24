/*
  code to do weak lensing power spectrum computations
  
  Matthew Becker, UofC 2011
*/

#include <gsl/gsl_spline.h>

#ifndef _WEAKLENS_
#define _WEAKLENS_

typedef struct {
  int wlNum;
  double lmin;
  double lmax;
  double tmin;
  double tmax;
} weaklensData;

extern weaklensData wlData;

typedef struct {
  int initFlag;
  int currCosmoNum;
  int currWLNum;
  double zs1;
  double zs2;
  double chis1;
  double chis2;
  double chiLim;
  double sn;
  double ell; //temp var for integrals
  gsl_spline *spline;
  gsl_interp_accel *accel;
} *lensPowerSpectra,_lensPowerSpectra;

typedef struct {
  int initFlag;
  int currCosmoNum;
  int currWLNum;
  lensPowerSpectra lps;
  double theta; //temp var for integrals
  gsl_spline *splineP;
  gsl_interp_accel *accelP;
  gsl_spline *splineM;
  gsl_interp_accel *accelM;
} *lensCorrFunc,_lensCorrFunc;

//functions to compute lensing power and cross spectra
double nonlinear_powspec_for_lens(double k, double a);
double lens_power_spectrum(double ell, lensPowerSpectra lps);
lensPowerSpectra init_lens_power_spectrum(double zs1, double zs2);
void free_lens_power_spectrum(lensPowerSpectra lps);

//functions to compute lens corr func.
double lens_corr_func_minus(double theta, lensCorrFunc lcf);
double lens_corr_func_plus(double theta, lensCorrFunc lcf);
lensCorrFunc init_lens_corr_func(lensPowerSpectra lps);
void free_lens_corr_func(lensCorrFunc lcf);

#endif /* _WEAKLENS_ */
