#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "cosmocalc.h"

static double linear_powspec_norm;
static gsl_spline *cosmocalc_sigmaR_spline = NULL;
static gsl_interp_accel *cosmocalc_sigmaR_acc = NULL;
static gsl_spline *cosmocalc_Rsigma_spline = NULL;
static gsl_interp_accel *cosmocalc_Rsigma_acc = NULL;

static void init_cosmocalc_peakheight_table(void);

static void init_cosmocalc_peakheight_table(void) 
{
  static int initFlag = 1;
  static int currCosmoNum;
  
  double sigmar_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  double topHatRad_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  size_t sindex[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  double tmpDouble[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH];
  long i;
  double rad;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);

      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        {
          rad = log(PH_R_MAX/PH_R_MIN)/(COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1.0)*((double) i) + log(PH_R_MIN);
          topHatRad_table[i] = rad;
          sigmar_table[i] = log(sqrt(linear_powspec_norm*tophatradnorm_linear_powspec_exact_nonorm(exp(rad))));
        }
      
      //init the spline and accelerators
      if(cosmocalc_sigmaR_spline != NULL)
        gsl_spline_free(cosmocalc_sigmaR_spline);
      cosmocalc_sigmaR_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_sigmaR_spline,topHatRad_table,sigmar_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      if(cosmocalc_sigmaR_acc != NULL)
        gsl_interp_accel_reset(cosmocalc_sigmaR_acc);
      else
        cosmocalc_sigmaR_acc = gsl_interp_accel_alloc();
      
      //reverse the table
      gsl_sort_index(sindex,sigmar_table,(size_t) 1,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        tmpDouble[i] = sigmar_table[sindex[i]];
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        sigmar_table[i] = tmpDouble[i];
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        tmpDouble[i] = topHatRad_table[sindex[i]];
      for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
        topHatRad_table[i] = tmpDouble[i];
      
      //init the spline and accelerators
      if(cosmocalc_Rsigma_spline != NULL)
        gsl_spline_free(cosmocalc_Rsigma_spline);
      cosmocalc_Rsigma_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_Rsigma_spline,sigmar_table,topHatRad_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
      if(cosmocalc_Rsigma_acc != NULL)
        gsl_interp_accel_reset(cosmocalc_Rsigma_acc);
      else
        cosmocalc_Rsigma_acc = gsl_interp_accel_alloc();
    }
}

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
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      init_cosmocalc_peakheight_table();
    }
  
  return exp(gsl_spline_eval(cosmocalc_sigmaR_spline,log(topHatRad),cosmocalc_sigmaR_acc))*growth_function(a);
}

double sigmaMtophat(double m, double a)
{
  return sigmaRtophat(pow(m/(4.0/3.0*M_PI*RHO_CRIT*cosmoData.OmegaM),1.0/3.0),a);
}

double inverse_sigmaRtophat(double sigmaR, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      init_cosmocalc_peakheight_table();
    }
  
  return exp(gsl_spline_eval(cosmocalc_Rsigma_spline,log(sigmaR/growth_function(a)),cosmocalc_Rsigma_acc));
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

