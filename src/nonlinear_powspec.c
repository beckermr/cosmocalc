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
static double nonlinear_gaussnorm_scale_funct(double gaussR, void *p);

static double gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = (*(double*)p);
  double k = exp(lnk);
  return linear_powspec(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad);
}

static double onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = (*(double*)p);
  double k = exp(lnk);
  return linear_powspec(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad)*(-1.0*k*k*2.0*gaussRad);
}

static double twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  double gaussRad = (*(double*)p);
  double k = exp(lnk);
  return  linear_powspec(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad)*(-2.0*k*k + 4.0*k*k*k*k*gaussRad*gaussRad);
}

//uses Takahashi et al. (2012) arXiv:1208.2701 unless macro below is set to use true Smith+03
//#define SMITH03

double nonlinear_powspec(double k, double a) 
{
  static int initFlag = 1;
  static int currCosmoNum;
  static gsl_spline *spline[4];
  static gsl_interp_accel *accel[4];
  int i;
  double xtab[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH],ytab[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH];
  double Rsigma,C,neff,ksigma,sigma2;
  double an,bn,cn,alphan,gamman,betan,mun,nun;
  double f1,f2,f3;
  double DeltakNL,dsigma2dR,d2sigma2d2R,PkNL,PkL;
  double y,DeltakL,fy,DeltakQ,DeltakHprime,DeltakH;
  //double t;
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
  double gaussRad;
  
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      currCosmoNum = cosmoData.cosmoNum;
      
      if(initFlag)
	{
	  for(i=0;i<4;++i)
	    spline[i] = gsl_spline_alloc(gsl_interp_akima,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
	  for(i=0;i<4;++i)
	    accel[i] = gsl_interp_accel_alloc();
          
	  initFlag = 0;
	}
      else
	{
	  for(i=0;i<4;++i)
	    gsl_spline_free(spline[i]);
	  for(i=0;i<4;++i)
	    spline[i] = gsl_spline_alloc(gsl_interp_akima,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
	  for(i=0;i<4;++i)
	    gsl_interp_accel_reset(accel[i]);
	}
      
      //t = -wtime();
      
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  xtab[i] = i*(1.0-AEXPN_MIN_NONLINPK)/(COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH-1.0) + AEXPN_MIN_NONLINPK;
	  ytab[i] = get_nonlinear_gaussnorm_scale(xtab[i]);
	}
      gsl_spline_init(spline[0],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      
      gsl_sort(ytab,(size_t) 1,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  xtab[i] = ytab[i];
	  ytab[i] = gaussiannorm_linear_powspec(xtab[i]);
	}
      gsl_spline_init(spline[1],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      
      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
      F.params = &gaussRad;
    
      F.function = &onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct;
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  gaussRad = xtab[i];
	  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
	  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
	  ytab[i] = I0+I1;
	}
      gsl_spline_init(spline[2],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      
      F.function = &twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct;
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  gaussRad = xtab[i];
	  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
	  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
	  ytab[i] = I0+I1;
	}
      gsl_spline_init(spline[3],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      
      gsl_integration_workspace_free(workspace);

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
      
      //t += wtime();
      //fprintf(stderr,"comp of non-linear Pk took %f seconds.\n",t);
    }

  Rsigma = gsl_spline_eval(spline[0],a,accel[0]);
  sigma2 = gsl_spline_eval(spline[1],Rsigma,accel[1]);
  dsigma2dR = gsl_spline_eval(spline[2],Rsigma,accel[2]);
  d2sigma2d2R = gsl_spline_eval(spline[3],Rsigma,accel[3]);
  
  ksigma = 1.0/Rsigma;
  neff = -1.0*Rsigma/sigma2*dsigma2dR - 3.0;
  C = -1.0*(d2sigma2d2R*Rsigma*Rsigma/sigma2 + dsigma2dR*Rsigma/sigma2 - dsigma2dR*dsigma2dR*Rsigma*Rsigma/sigma2/sigma2);
  
  double ha,weffa,omegaMz,omegaDEwz;
  
  ha = hubble_noscale(a);
  weffa = weff(a);
  omegaMz = cosmoData.OmegaM/a/a/a/ha/ha;
  omegaDEwz = (1.0-cosmoData.OmegaM)/ha/ha/pow(a,3.0*(1.0 + weffa));
  
#ifdef SMITH03
  an = pow(10.0,1.4861 + 1.8369*neff + 1.6762*neff*neff + 0.7940*neff*neff*neff + 0.1670*neff*neff*neff*neff - 0.6206*C);
  bn = pow(10.0,0.9463 + 0.9466*neff + 0.3084*neff*neff - 0.9400*C);
  cn = pow(10.0,-0.2807 + 0.6669*neff + 0.3214*neff*neff - 0.0793*C);
  gamman = 0.8649 + 0.2989*neff + 0.1631*C;
  alphan = 1.3884 + 0.3700*neff - 0.1452*neff*neff;
  betan = 0.8291 + 0.9854*neff + 0.3401*neff*neff;
  mun = pow(10.0,-3.5442 + 0.1908*neff);
  nun = pow(10.0,0.9589 + 1.2857*neff);
  
  f1 = pow(OmegaMz,-0.0307);
  f2 = pow(OmegaMz,-0.0585);
  f3 = pow(OmegaMz,0.0743);
#else
  an = pow(10.0,1.5222 + 2.8553*neff + 2.3706*neff*neff + 0.9903*neff*neff*neff + 0.2250*neff*neff*neff*neff - 0.6038*C + 0.1749*omegaDEwz*(1.0 + weffa));
  bn = pow(10.0,-0.5642 + 0.5864*neff + 0.5716*neff*neff - 1.5474*C + 0.2279*omegaDEwz*(1.0 + weffa));
  cn = pow(10.0,0.3698 + 2.0404*neff + 0.8161*neff*neff + 0.5869*C);
  gamman = 0.1971 - 0.0843*neff + 0.8460*C;
  alphan = fabs(6.0835 + 1.3373*neff - 0.1959*neff*neff - 5.5274*C);
  betan = 2.0379 - 0.7354*neff + 0.3157*neff*neff + 1.2490*neff*neff*neff + 0.3980*neff*neff*neff*neff - 0.1682*C;
  mun = 0.0;
  nun = pow(10.0,5.2105 + 3.6902*neff);
  
  f1 = pow(omegaMz,-0.0307);
  f2 = pow(omegaMz,-0.0585);
  f3 = pow(omegaMz,0.0743);
#endif
  
  PkL = linear_powspec(k,a);
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

double gaussiannorm_linear_powspec(double gaussRad)
{
#define PNL_RGAUSS_MIN 0.0001
#define PNL_RGAUSS_MAX 100.0
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0 
  
  static int initFlag = 1;
  static int currCosmoNum;
  static gsl_spline *spline = NULL;
  static gsl_interp_accel *accel = NULL;
  
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
  int i;
  double xtab[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH];
  double ytab[COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH];
  double dlnr;
  double lnrmin;
  double lnr;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;

      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
      F.function = &gaussiannorm_linear_powspec_exact_lnk_integ_funct;
      
      dlnr = log(PNL_RGAUSS_MAX/PNL_RGAUSS_MIN)/(COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH-1.0);
      lnrmin = log(PNL_RGAUSS_MIN);
      
      for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
	{
	  lnr = dlnr*i + lnrmin;
	  xtab[i] = exp(lnr);
	  F.params = &(xtab[i]);

	  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
	  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
	  
	  //gsl_integration_qagil(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
	  //gsl_integration_qagiu(&F,log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);

	  xtab[i] = lnr;
	  ytab[i] = log(I0+I1);
	}
      
      gsl_integration_workspace_free(workspace);
      
      if(spline != NULL)
	gsl_spline_free(spline);
      spline = gsl_spline_alloc(GSL_SPLINE_TYPE,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      gsl_spline_init(spline,xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      if(accel != NULL)
	gsl_interp_accel_reset(accel);
      else
	accel = gsl_interp_accel_alloc();
      
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
#undef PNL_RGAUSS_MIN
#undef PNL_RGAUSS_MAX
    }
  
  return exp(gsl_spline_eval(spline,log(gaussRad),accel));
}

static double nonlinear_gaussnorm_scale_funct(double gaussR, void *p)
{
  double gf = ((double*)p)[0];
  
  return gaussiannorm_linear_powspec(gaussR)*gf*gf-1.0;
}

double get_nonlinear_gaussnorm_scale(double a)
{
  double gf = growth_function(a);
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


