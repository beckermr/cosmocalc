#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "cosmocalc.h"

#define COSMOCALC_PEAKHEIGHT_TABLE_LENGTH 400
#define R_MIN 1e-3
#define R_MAX 1e3

double sigmaRtophat_exact(double topHatRad, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm;
  double gf;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
    }
  
  gf = growth_function(a);
  
  return sqrt(linear_powspec_norm*tophatradnorm_linear_powspec_exact_nonorm(topHatRad))*gf;
}

double sigmaRtophat(double topHatRad, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm;
  static gsl_spline *cosmocalc_sigma_spline = NULL;
  static gsl_interp_accel *cosmocalc_sigma_acc = NULL;
  
  double sigmar_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  double topHatRad_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  long i;
  double rad,gf;
    
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
      
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        {
          rad = log(R_MAX/R_MIN)/(COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1.0)*((double) i) + log(R_MIN);
          topHatRad_table[i] = rad;
	  sigmar_table[i] = log(sqrt(linear_powspec_norm*tophatradnorm_linear_powspec_exact_nonorm(exp(rad))));
	}
      
      //init the spline and accelerators
      if(cosmocalc_sigma_spline != NULL)
        gsl_spline_free(cosmocalc_sigma_spline);
      cosmocalc_sigma_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_sigma_spline,topHatRad_table,sigmar_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      if(cosmocalc_sigma_acc != NULL)
        gsl_interp_accel_reset(cosmocalc_sigma_acc);
      else
        cosmocalc_sigma_acc = gsl_interp_accel_alloc();
    }

  gf = growth_function(a);
  
  return exp(gsl_spline_eval(cosmocalc_sigma_spline,log(topHatRad),cosmocalc_sigma_acc))*gf;
}

double sigmaMtophat(double m, double a)
{
  return sigmaRtophat(pow(m/(4.0/3.0*M_PI*RHO_CRIT*cosmoData.OmegaM),1.0/3.0),a);
}

double inverse_sigmaRtophat(double sigmaR, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm;
  static gsl_spline *cosmocalc_sigma_spline = NULL;
  static gsl_interp_accel *cosmocalc_sigma_acc = NULL;
  
  double sigmar_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  double topHatRad_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  long i;
  double rad,gf;
    
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
      
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        {
          rad = log(R_MAX/R_MIN)/(COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1.0)*((double) i) + log(R_MIN);
          topHatRad_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1-i] = rad;
	  sigmar_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1-i] = log(sqrt(linear_powspec_norm*tophatradnorm_linear_powspec_exact_nonorm(exp(rad))));
	}
            
      //init the spline and accelerators
      if(cosmocalc_sigma_spline != NULL)
        gsl_spline_free(cosmocalc_sigma_spline);
      cosmocalc_sigma_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_sigma_spline,sigmar_table,topHatRad_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      if(cosmocalc_sigma_acc != NULL)
        gsl_interp_accel_reset(cosmocalc_sigma_acc);
      else
        cosmocalc_sigma_acc = gsl_interp_accel_alloc();
    }
  
  gf = growth_function(a);
    
  return exp(gsl_spline_eval(cosmocalc_sigma_spline,log(sigmaR/gf),cosmocalc_sigma_acc));
}

double inverse_sigmaMtophat(double sigmaR, double a)
{
  return (4.0/3.0*M_PI*RHO_CRIT*cosmoData.OmegaM)*pow(inverse_sigmaRtophat(sigmaR,a),3.0);
}

double inverse_nuMtophat(double nu, double a)
{
  return (4.0/3.0*M_PI*RHO_CRIT*cosmoData.OmegaM)*pow(inverse_sigmaRtophat(DELTAC/nu,a),3.0);
}

double inverse_nuRtophat(double nu, double a)
{
  return inverse_sigmaRtophat(DELTAC/nu,a);
}

#undef COSMOCALC_PEAKHEIGHT_TABLE_LENGTH
#undef R_MIN
#undef R_MAX
