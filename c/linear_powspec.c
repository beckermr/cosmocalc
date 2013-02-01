#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_roots.h>

#include "cosmocalc.h"

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
  double ft = fourierTransformTopHat(k*(*((double*)p)));
  double tf = transfer_function(k);
  return ft*ft*tf*tf*pow(k,cosmoData.SpectralIndex)*k*k*k/2.0/M_PI/M_PI*k;
}

static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p)
{
  // topHatRad = (*((double*)p));
  double ft = fourierTransformTopHat(k*(*((double*)p)));
  double tf = transfer_function(k);
  return ft*ft*tf*tf*pow(k,cosmoData.SpectralIndex)*k*k/2.0/M_PI/M_PI;
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
  static double linear_powspec_norm = 1.0;
  double gf = growth_function(a);
  double tf = transfer_function(k);
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;

      linear_powspec_norm = cosmoData.Sigma8*cosmoData.Sigma8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
    }
  
  return tf*tf*pow(k,cosmoData.SpectralIndex)*gf*gf*linear_powspec_norm;
}

double linear_powspec(double k, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double linear_powspec_norm = 1.0;
  static gsl_spline *cosmocalc_linear_powspec_spline = NULL;
  static gsl_interp_accel *cosmocalc_linear_powspec_acc = NULL; 
  static double c0,c1;
  
  double linear_powspec_table[COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH];
  double k_table[COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH];
  long i;
  double gf,cov00,cov01,cov11,sumsq,tf;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      for(i=0;i<COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  k_table[i] = log(K_MAX/K_MIN)/(COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH-1.0)*((double) i) + log(K_MIN);
	  tf = transfer_function(exp(k_table[i]));
	  linear_powspec_table[i] = log(tf*tf*pow(exp(k_table[i]),cosmoData.SpectralIndex));
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
    {
      tf = transfer_function(k);
      return tf*tf*pow(k,cosmoData.SpectralIndex)*linear_powspec_norm*gf*gf;
    }
  else if(k < K_MAX)
    return exp(gsl_spline_eval(cosmocalc_linear_powspec_spline,log(k),cosmocalc_linear_powspec_acc))*linear_powspec_norm*gf*gf;
  else
    return exp(c0+c1*log(k))*linear_powspec_norm*gf*gf;
}

double tophatnorm_linear_powspec(double topHatRad)
{
#define NL_RTOPHAT_MIN 0.0001
#define NL_RTOPHAT_MAX 100.0
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-8
#define RELERR 0.0 
  
  static int initFlag = 1;
  static int currCosmoNum;
  static gsl_spline *spline = NULL;
  static gsl_interp_accel *accel = NULL;
  
  double I0,I1;
  double abserr,epsabs,epsrel;
  gsl_integration_workspace *workspace;
  gsl_function F;
  int i;
  double xtab[COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH];
  double ytab[COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH];
  double dlnr;
  double lnrmin;
  double lnr;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;

      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
      
      dlnr = log(NL_RTOPHAT_MAX/NL_RTOPHAT_MIN)/(COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH-1.0);
      lnrmin = log(NL_RTOPHAT_MIN);
      
      for(i=0;i<COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH;++i)
        {
          lnr = dlnr*i + lnrmin;
          xtab[i] = exp(lnr);
	  
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
	  
          xtab[i] = lnr;
          ytab[i] = log(I0+I1);
        }
      
      gsl_integration_workspace_free(workspace);
      
      if(spline != NULL)
        gsl_spline_free(spline);
      spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH));
      gsl_spline_init(spline,xtab,ytab,(size_t) (COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH));
      if(accel != NULL)
        gsl_interp_accel_reset(accel);
      else
        accel = gsl_interp_accel_alloc();
      
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
#undef NL_RTOPHAT_MIN
#undef NL_RTOPHAT_MAX
    }
  
  return exp(gsl_spline_eval(spline,log(topHatRad),accel));
}

static double linear_tophatnorm_scale_funct(double rad, void *p)
{
  double gf = ((double*)p)[0];
  
  return tophatnorm_linear_powspec(rad)*gf*gf-1.0;
}

double get_linear_tophatnorm_scale(double a)
{
  double gf = growth_function(a);
  double Rsigma,Rlow=0.001,Rhigh=10.0;
  int itr,maxItr=1000,status;
  
#define ABSERR 1e-6
#define RELERR 1e-6
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
       
  F.function = &linear_tophatnorm_scale_funct;
  F.params = &gf;
  
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s,&F,Rlow,Rhigh);
  itr = 0;
  
  do
    {
      itr++;
      status = gsl_root_fsolver_iterate(s);
      Rsigma = gsl_root_fsolver_root(s);
      Rlow = gsl_root_fsolver_x_lower(s);
      Rhigh = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(Rlow,Rhigh,ABSERR,RELERR);
    }
  while(status == GSL_CONTINUE && itr < maxItr);
  
#undef ABSERR
#undef RELERR

  gsl_root_fsolver_free(s);
  
  return Rsigma;
}
