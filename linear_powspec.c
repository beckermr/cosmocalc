#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_fit.h>

#include "cosmocalc.h"

#define COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH 1000
#define K_MIN 1e-9
#define K_MAX 1e20
#define COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH 20

static double fourierTransformTopHat(double y);
static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p);
static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p);

static double fourierTransformTopHat(double y)
{
  if(y < 1e-3)
    return 1.0;
  else
    return 3.0/y/y/y*(sin(y) - y*cos(y));
}

static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p)
{
  // topHatRad = (*((double*)p));
  double k = exp(lnk);
  return fourierTransformTopHat(k*(*((double*)p)))*fourierTransformTopHat(k*(*((double*)p)))
    *transfer_function(k)*transfer_function(k)*pow(k,cosmoData.SpectralIndex)
    *k*k/2.0/M_PI/M_PI*k;
}

static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p)
{
  // topHatRad = (*((double*)p));
  return fourierTransformTopHat(k*(*((double*)p)))*fourierTransformTopHat(k*(*((double*)p)))
    *transfer_function(k)*transfer_function(k)*pow(k,cosmoData.SpectralIndex)
    *k*k/2.0/M_PI/M_PI;
}

double tophatradnorm_linear_powspec_exact_nonorm(double topHatRad)
{
  double I0,I1;
  double abserr;
  double epsrel,epsabs;
  gsl_integration_workspace *workspace;
  gsl_function F;
    
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-8
#define RELERR 0.0 
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  F.params = &(topHatRad);
  if(topHatRad > 1e-4)
    {
      epsabs = 1e-20;
      epsrel = 1e-6;
      F.function = &tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0;
      gsl_integration_qags(&F,0.0,2.0*M_PI/topHatRad,epsabs,epsrel,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
      gsl_integration_qagiu(&F,2.0*M_PI/topHatRad,epsabs,epsrel,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
    }
  else
    {
      F.function = &tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0;
      gsl_integration_qagil(&F,log(2.0*M_PI/topHatRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
      gsl_integration_qagiu(&F,log(2.0*M_PI/topHatRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
    }
  
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return I0 + I1;
}

double linear_powspec_exact(double k, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm;
  double gf = growth_function(a);
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;

      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
    }
  
  return transfer_function(k)*transfer_function(k)*pow(k,cosmoData.SpectralIndex)*gf*gf*linear_powspec_norm;
}

double linear_powspec(double k, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm;
  static gsl_spline *cosmocalc_linear_powspec_spline = NULL;
  static gsl_interp_accel *cosmocalc_linear_powspec_acc = NULL; 
  static double c0,c1;
  
  double linear_powspec_table[COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH];
  double k_table[COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH];
  long i;
  double gf,cov00,cov01,cov11,sumsq;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      for(i=0;i<COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  k_table[i] = log(K_MAX/K_MIN)/(COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH-1.0)*((double) i) + log(K_MIN);
	  linear_powspec_table[i] = log(transfer_function(exp(k_table[i]))*transfer_function(exp(k_table[i]))*pow(exp(k_table[i]),cosmoData.SpectralIndex));
	}
            
      //init the spline and accelerators
      if(cosmocalc_linear_powspec_spline != NULL)
	gsl_spline_free(cosmocalc_linear_powspec_spline);
      cosmocalc_linear_powspec_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_linear_powspec_spline,k_table,linear_powspec_table,(size_t) (COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH));
      if(cosmocalc_linear_powspec_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_linear_powspec_acc);
      else
	cosmocalc_linear_powspec_acc = gsl_interp_accel_alloc();
      
      gsl_fit_linear(k_table+COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH-COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		     linear_powspec_table+COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH-COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		     (size_t) COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH,&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
      
      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
    }
  
  gf = growth_function(a);
  
  if(k < K_MIN)
    return transfer_function(k)*transfer_function(k)*pow(k,cosmoData.SpectralIndex)*linear_powspec_norm*gf*gf;
  else if(k < K_MAX)
    return exp(gsl_spline_eval(cosmocalc_linear_powspec_spline,log(k),cosmocalc_linear_powspec_acc))*linear_powspec_norm*gf*gf;
  else
    return exp(c0+c1*log(k))*linear_powspec_norm*gf*gf;
}

#undef COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH
#undef K_MIN
#undef K_MAX
#undef COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH
