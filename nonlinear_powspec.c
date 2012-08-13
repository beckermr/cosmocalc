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

static double gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
static double onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
static double twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
static double gaussiannorm_linear_powspec_exact(double gaussRad);
static double onederiv_gaussiannorm_linear_powspec_exact(double gaussRad);
static double twoderiv_gaussiannorm_linear_powspec_exact(double gaussRad);
static double nonlinear_gaussnorm_scale_funct(double gaussR, void *p);
static double get_nonlinear_gaussnorm_scale(double a);

static double gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = ((double*)p)[0];
  double k = exp(lnk);
  return linear_powspec_exact(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad);
}

static double gaussiannorm_linear_powspec_exact(double gaussRad)
{
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
    
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0 
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  F.params = &(gaussRad);
  F.function = &gaussiannorm_linear_powspec_exact_lnk_integ_funct;
  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  //gsl_integration_qagil(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  //gsl_integration_qagiu(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return I0 + I1;
}

static double onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = ((double*)p)[0];
  double k = exp(lnk);
  return linear_powspec_exact(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad)*(-1.0*k*k*2.0*gaussRad);
}

static double onederiv_gaussiannorm_linear_powspec_exact(double gaussRad)
{
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
    
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0 
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  F.params = &(gaussRad);
  F.function = &onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct;
  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  //gsl_integration_qagil(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  //gsl_integration_qagiu(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return I0 + I1;
}

static double twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = ((double*)p)[0];
  double k = exp(lnk);
  return linear_powspec_exact(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad)*(-2.0*k*k + 4.0*k*k*k*k*gaussRad*gaussRad);
}

static double twoderiv_gaussiannorm_linear_powspec_exact(double gaussRad)
{
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
    
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0 
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  F.params = &(gaussRad);
  F.function = &twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct;
  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  //gsl_integration_qagil(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  //gsl_integration_qagiu(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
  
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return I0 + I1;
}

static double nonlinear_gaussnorm_scale_funct(double gaussR, void *p)
{
  double gf = ((double*)p)[0];
  
  return gaussiannorm_linear_powspec_exact(gaussR)*gf*gf-1.0;
}

static double get_nonlinear_gaussnorm_scale(double a)
{
  double gf = growth_function_exact(a);
  double Rsigma,Rlow=0.001,Rhigh=10.0;
  int itr,maxItr=1000,status;
  
#define ABSERR 1e-6
#define RELERR 1e-6
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
       
  F.function = &nonlinear_gaussnorm_scale_funct;
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

double nonlinear_powspec_exact(double k, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double aint;
  static double Rsigma,C,neff,ksigma,sigma2;
  static double an,bn,cn,alphan,gamman,betan,mun,nun;
  static double f1,f2,f3;
  double DeltakNL,dsigma2dR,d2sigma2d2R,PkNL,PkL;
  double y,DeltakL,fy,DeltakQ,DeltakHprime,DeltakH;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum || aint != a)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      aint = a;
      
      Rsigma = get_nonlinear_gaussnorm_scale(aint);
      ksigma = 1.0/Rsigma;
      sigma2 = gaussiannorm_linear_powspec_exact(Rsigma);
      dsigma2dR = onederiv_gaussiannorm_linear_powspec_exact(Rsigma);
      d2sigma2d2R = twoderiv_gaussiannorm_linear_powspec_exact(Rsigma);
      neff = -1.0*Rsigma/sigma2*dsigma2dR - 3.0;
      C = -1.0*(d2sigma2d2R*Rsigma*Rsigma/sigma2 + dsigma2dR*Rsigma/sigma2 - dsigma2dR*dsigma2dR*Rsigma*Rsigma/sigma2/sigma2);
      
      //FIXME fprintf(stderr,"a = %f, sigma = %f, ksigma = %f, neff = %f, C = %f\n",a,sqrt(sigma2),ksigma,neff,C);
      
      an = pow(10.0,1.4861 + 1.8369*neff + 1.6762*neff*neff + 0.7940*neff*neff*neff + 0.1670*neff*neff*neff*neff - 0.6206*C);
      bn = pow(10.0,0.9463 + 0.9466*neff + 0.3084*neff*neff - 0.9400*C);
      cn = pow(10.0,-0.2807 + 0.6669*neff + 0.3214*neff*neff - 0.0793*C);
      gamman = 0.8649 + 0.2989*neff + 0.1631*C;
      alphan = 1.3884 + 0.3700*neff - 0.1452*neff*neff;
      betan = 0.8291 + 0.9854*neff + 0.3401*neff*neff;
      mun = pow(10.0,-3.5442 + 0.1908*neff);
      nun = pow(10.0,0.9589 + 1.2857*neff);
      
      f1 = pow(cosmoData.OmegaM,-0.0307);
      f2 = pow(cosmoData.OmegaM,-0.0585);
      f3 = pow(cosmoData.OmegaM,0.0743);
    }
  
  PkL = linear_powspec_exact(k,a);
  y = k/ksigma;
  fy = y/4.0 + y*y/8.0;
  DeltakL = PkL*k*k*k/2.0/M_PI/M_PI;
  
  DeltakQ = DeltakL*pow(1.0 + DeltakL,betan)/(1.0 + alphan*DeltakL)*exp(-1.0*fy);
  
  DeltakHprime = an*pow(y,3.0*f1)/(1.0 + bn*pow(y,f2) + pow(cn*f3*y,3.0 - gamman));
  DeltakH = DeltakHprime/(1.0 + mun/y + nun/y/y);
  
  DeltakNL = DeltakQ + DeltakH;
  PkNL = DeltakNL/(k*k*k/2.0/M_PI/M_PI);
  
  return PkNL;
}

#define COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH 1000
#define K_MIN 1e-7
#define K_MAX 1e5
#define COSMOCALC_NONLINEAR_POWSPEC_FIT_LENGTH 20

double nonlinear_powspec(double k, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double aint;
  static gsl_spline *cosmocalc_nonlinear_powspec_spline = NULL;
  static gsl_interp_accel *cosmocalc_nonlinear_powspec_acc = NULL; 
    
  double nonlinear_powspec_table[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH];
  double k_table[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH];
  long i;
    
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum || a != aint)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      aint = a;
      
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  k_table[i] = log(K_MAX/K_MIN)/(COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH-1.0)*((double) i) + log(K_MIN);
	  nonlinear_powspec_table[i] = log(nonlinear_powspec_exact(exp(k_table[i]),aint));
	}
      
      //init the spline and accelerators
      if(cosmocalc_nonlinear_powspec_spline != NULL)
	gsl_spline_free(cosmocalc_nonlinear_powspec_spline);
      cosmocalc_nonlinear_powspec_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_nonlinear_powspec_spline,k_table,nonlinear_powspec_table,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      if(cosmocalc_nonlinear_powspec_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_nonlinear_powspec_acc);
      else
	cosmocalc_nonlinear_powspec_acc = gsl_interp_accel_alloc();
  
    }
  
  if(k < K_MIN)
    return nonlinear_powspec_exact(k,aint);
  else if(k < K_MAX)
    return exp(gsl_spline_eval(cosmocalc_nonlinear_powspec_spline,log(k),cosmocalc_nonlinear_powspec_acc));
  else
    return nonlinear_powspec_exact(k,aint);
  }

#undef COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH
#undef K_MIN
#undef K_MAX
#undef COSMOCALC_NONLINEAR_POWSPEC_FIT_LENGTH
