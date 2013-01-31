#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>

#include "cosmocalc.h"
#include "weaklens.h"

double chiLim;

static void comp_lens_power_spectrum(lensPowerSpectra lps);

double lens_power_spectrum(double ell, lensPowerSpectra lps)
{
  if(lps->initFlag == 1 || lps->currCosmoNum != cosmoData.cosmoNum || lps->currWLNum != wlData.wlNum)
    {
      lps->initFlag = 0;
      lps->currCosmoNum = cosmoData.cosmoNum;
      lps->currWLNum = wlData.wlNum;
      
      comp_lens_power_spectrum(lps);
    }
  
  return exp(gsl_spline_eval(lps->spline,log(ell),lps->accel));
}

double nonlinear_powspec_for_lens(double k, double a)
{
  return nonlinear_powspec(k,a);
}

static double lenskern(double chi, double chis)
{
  if(chi > chis)
    return 0.0;
  else
    return 1.5*cosmoData.OmegaM*(100.0/CSOL)*(100.0/CSOL)/acomvdist(chi)*(chis-chi)/chis;
}

static double lenspk_integrand(double chi, void *p)
{
  lensPowerSpectra lps = (lensPowerSpectra) p;
  double sn = lps->sn;
  
  if(chi == 0.0 || chi < chiLim)
    return 0.0;
  else
    return lenskern(chi,lps->chis1)*lenskern(chi,lps->chis2)*(nonlinear_powspec_for_lens(lps->ell/chi,acomvdist(chi)) + sn);
}

static void comp_lens_power_spectrum(lensPowerSpectra lps)
{
#define WORKSPACE_NUM 100000
#define ABSERR 1e-12
#define RELERR 1e-12
#define TABLE_LENGTH 1000
  
  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr;
  double logltab[TABLE_LENGTH];
  double logpkltab[TABLE_LENGTH];
  double chimax;
  int i;
  
  //fill in bin information
  chiLim = lps->chiLim;
  if(lps->chis1 > lps->chis2)
    chimax = lps->chis1;
  else
    chimax = lps->chis2;
  
  fprintf(stderr,"doing lens pk - chiLim = %lf, chiMax = %lf\n",chiLim,chimax);
  
  //init
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  F.function = &lenspk_integrand;
  F.params = lps;
  
  //make table
  double lnlmin = log(wlData.lmin);
  double lnlmax = log(wlData.lmax);
  for(i=0;i<TABLE_LENGTH;++i)
    {
      logltab[i] = i*(lnlmax-lnlmin)/(TABLE_LENGTH-1) + lnlmin;
      
      lps->ell = exp(logltab[i]);
      gsl_integration_qag(&F,0.0,chimax,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
      
      logpkltab[i] = log(result);
    }
  
  //free
  gsl_integration_workspace_free(workspace);
  
  //init splines and accels
  if(lps->spline != NULL)
    gsl_spline_free(lps->spline);
  lps->spline = gsl_spline_alloc(gsl_interp_akima,(size_t) (TABLE_LENGTH));
  gsl_spline_init(lps->spline,logltab,logpkltab,(size_t) (TABLE_LENGTH));
  if(lps->accel != NULL)
    gsl_interp_accel_reset(lps->accel);
  else
    lps->accel = gsl_interp_accel_alloc();
    
#undef TABLE_LENGTH
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
}

lensPowerSpectra init_lens_power_spectrum(double zs1, double zs2)
{
  lensPowerSpectra lps;
  
  lps = (lensPowerSpectra)malloc(sizeof(_lensPowerSpectra));
  assert(lps != NULL);
  
  lps->initFlag = 1;
  lps->zs1 = zs1;
  lps->zs2 = zs2;
  lps->chis1 = comvdist(1.0/(1.0 + zs1));
  lps->chis2 = comvdist(1.0/(1.0 + zs2));
  lps->spline = NULL;
  lps->accel = NULL;
  
  return lps;
}

void free_lens_power_spectrum(lensPowerSpectra lps)
{
  if(lps->spline != NULL)
    gsl_spline_free(lps->spline);
  if(lps->accel != NULL)
    gsl_interp_accel_free(lps->accel);
  free(lps);
}

////////////////////////////////////////
// corr. funcs!
//////////////////////////////////////

static double lenscfp_integrand(double ell, void *p)
{
  lensCorrFunc lcf = (lensCorrFunc) p;
  
  return ell/2.0/M_PI*lens_power_spectrum(ell,lcf->lps)*gsl_sf_bessel_J0(ell*lcf->theta/60.0/180.0*M_PI);
  
}

static double lenscfm_integrand(double ell, void *p)
{
  lensCorrFunc lcf = (lensCorrFunc) p;
  
  return ell/2.0/M_PI*lens_power_spectrum(ell,lcf->lps)*gsl_sf_bessel_Jn(4,ell*lcf->theta/60.0/180.0*M_PI);
}

static void comp_lens_corr_funcs(lensCorrFunc lcf)
{
#define WORKSPACE_NUM 100000
#define ABSERR 1e-12
#define RELERR 1e-12
#define TABLE_LENGTH 1000
  
  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr;
  double logttab[TABLE_LENGTH];
  double logcfptab[TABLE_LENGTH];
  double logcfmtab[TABLE_LENGTH];
  int i;
  double lntmin;
  double lntmax;
  
  //init
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  F.params = lcf;
  lntmin = log(wlData.tmin);
  lntmax = log(wlData.tmax);
  
  //make tables
  F.function = &lenscfp_integrand;
  for(i=0;i<TABLE_LENGTH;++i)
    {
      logttab[i] = i*(lntmax-lntmin)/(TABLE_LENGTH-1) + lntmin;
      
      lcf->theta = exp(logttab[i]);
      gsl_integration_qag(&F,wlData.lmin,wlData.lmax,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
      
      logcfptab[i] = log(result);
    }
  
  F.function = &lenscfm_integrand;
  for(i=0;i<TABLE_LENGTH;++i)
    {
      logttab[i] = i*(lntmax-lntmin)/(TABLE_LENGTH-1) + lntmin;
      
      lcf->theta = exp(logttab[i]);
      gsl_integration_qag(&F,wlData.lmin,wlData.lmax,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
      
      logcfmtab[i] = log(result);
    }
    
  //free
  gsl_integration_workspace_free(workspace);
  
  //init splines and accels
  if(lcf->splineP != NULL)
    gsl_spline_free(lcf->splineP);
  lcf->splineP = gsl_spline_alloc(gsl_interp_akima,(size_t) (TABLE_LENGTH));
  gsl_spline_init(lcf->splineP,logttab,logcfptab,(size_t) (TABLE_LENGTH));
  if(lcf->accelP != NULL)
    gsl_interp_accel_reset(lcf->accelP);
  else
    lcf->accelP = gsl_interp_accel_alloc();
  
  if(lcf->splineM != NULL)
    gsl_spline_free(lcf->splineM);
  lcf->splineM = gsl_spline_alloc(gsl_interp_akima,(size_t) (TABLE_LENGTH));
  gsl_spline_init(lcf->splineM,logttab,logcfmtab,(size_t) (TABLE_LENGTH));
  if(lcf->accelM != NULL)
    gsl_interp_accel_reset(lcf->accelM);
  else
    lcf->accelM = gsl_interp_accel_alloc();
      
#undef TABLE_LENGTH
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM  
}

double lens_corr_func_minus(double theta, lensCorrFunc lcf)
{
  if(lcf->initFlag == 1 || lcf->currCosmoNum != cosmoData.cosmoNum || lcf->currWLNum != wlData.wlNum)
    {
      lcf->initFlag = 0;
      lcf->currCosmoNum = cosmoData.cosmoNum;
      lcf->currWLNum = wlData.wlNum;

      comp_lens_corr_funcs(lcf);
    }

  return exp(gsl_spline_eval(lcf->splineM,log(theta),lcf->accelM));
}

double lens_corr_func_plus(double theta, lensCorrFunc lcf)
{
  if(lcf->initFlag == 1 || lcf->currCosmoNum != cosmoData.cosmoNum || lcf->currWLNum != wlData.wlNum)
    {
      lcf->initFlag = 0;
      lcf->currCosmoNum = cosmoData.cosmoNum;
      lcf->currWLNum = wlData.wlNum;

      comp_lens_corr_funcs(lcf);
    }

  return exp(gsl_spline_eval(lcf->splineP,log(theta),lcf->accelP));
}

lensCorrFunc init_lens_corr_func(lensPowerSpectra lps)
{
  lensCorrFunc lcf;
  
  lcf = (lensCorrFunc)malloc(sizeof(_lensCorrFunc));
  assert(lcf != NULL);
  
  lcf->initFlag = 1;
  lcf->lps = lps;
  lcf->splineM = NULL;
  lcf->accelM = NULL;
  lcf->splineP = NULL;
  lcf->accelP = NULL;
  
  return lcf;
}

void free_lens_corr_func(lensCorrFunc lcf)
{
  if(lcf->splineM != NULL)
    gsl_spline_free(lcf->splineM);
  if(lcf->accelM != NULL)
    gsl_interp_accel_free(lcf->accelM);
    
  if(lcf->splineP != NULL)
    gsl_spline_free(lcf->splineP);
  if(lcf->accelP != NULL)
    gsl_interp_accel_free(lcf->accelP);
  
  free(lcf);
}
