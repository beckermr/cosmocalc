#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "cosmocalc.h"

gsl_spline *cosmocalc_aexpn2comvdist_spline = NULL;
gsl_interp_accel *cosmocalc_aexpn2comvdist_acc = NULL; 
gsl_spline *cosmocalc_comvdist2aexpn_spline = NULL;
gsl_interp_accel *cosmocalc_comvdist2aexpn_acc = NULL; 
double DH_sqrtok;

/* function for integration using gsl integration */
static double comvdist_integ_funct(double a, void *p)
{
  return 1.0/a/a/hubble_noscale(a);
}

double comvdist_exact(double a)
{
#define WORKSPACE_NUM 100000
#define ABSERR 0.0
#define RELERR 1e-8

  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr;
  
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  F.function = &comvdist_integ_funct;
  gsl_integration_qag(&F,a,1.0,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
  
  gsl_integration_workspace_free(workspace);

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return result*DH;
}

/* init function  - some help from Gadget-2 applied here */
void init_cosmocalc_comvdist_table(void)
{
#define WORKSPACE_NUM 100000
#define ABSERR 0.0
#define RELERR 1e-8

  static int initFlag = 1;
  static int currCosmoNum;
  
  gsl_integration_workspace *workspace;
  gsl_function F;
  long i;
  double result,abserr,afact;
  double comvdist_table[COSMOCALC_COMVDIST_TABLE_LENGTH];
  double aexpn_table[COSMOCALC_COMVDIST_TABLE_LENGTH];
  size_t sindex[COSMOCALC_COMVDIST_TABLE_LENGTH];
  double tmpDouble[COSMOCALC_COMVDIST_TABLE_LENGTH];
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      if(cosmoData.OmegaK != 0.0)
	DH_sqrtok = DH/sqrt(cosmoData.OmegaK);
      else
	DH_sqrtok = 1.0;
      
      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
      for(i=0;i<COSMOCALC_COMVDIST_TABLE_LENGTH-1;++i)
	{
	  afact = (AEXPN_MAX - AEXPN_MIN)/(COSMOCALC_COMVDIST_TABLE_LENGTH-1.0)*((double) i) + AEXPN_MIN;
	  F.function = &comvdist_integ_funct;
	  F.params = &(cosmoData.OmegaM);
	  gsl_integration_qag(&F,afact,1.0,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
	  aexpn_table[i] = afact;
	  comvdist_table[i] = result*DH;
	}
      aexpn_table[i] = 1.0;
      comvdist_table[i] = 0.0;
      
      gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
      
      //init the spline and accelerators
      if(cosmocalc_aexpn2comvdist_spline != NULL)
	gsl_spline_free(cosmocalc_aexpn2comvdist_spline);
      cosmocalc_aexpn2comvdist_spline = gsl_spline_alloc(GSL_SPLINE_TYPE,(size_t) (COSMOCALC_COMVDIST_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_aexpn2comvdist_spline,aexpn_table,comvdist_table,(size_t) (COSMOCALC_COMVDIST_TABLE_LENGTH));
      if(cosmocalc_aexpn2comvdist_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_aexpn2comvdist_acc);
      else
	cosmocalc_aexpn2comvdist_acc = gsl_interp_accel_alloc();
      
      //sort the distances properly
      gsl_sort_index(sindex,comvdist_table,(size_t) 1,(size_t) (COSMOCALC_COMVDIST_TABLE_LENGTH));
      for(i=0;i<COSMOCALC_COMVDIST_TABLE_LENGTH;++i)
	tmpDouble[i] = comvdist_table[sindex[i]];
      for(i=0;i<COSMOCALC_COMVDIST_TABLE_LENGTH;++i)
	comvdist_table[i] = tmpDouble[i];
      for(i=0;i<COSMOCALC_COMVDIST_TABLE_LENGTH;++i)
	tmpDouble[i] = aexpn_table[sindex[i]];
      for(i=0;i<COSMOCALC_COMVDIST_TABLE_LENGTH;++i)
	aexpn_table[i] = tmpDouble[i];
      
      if(cosmocalc_comvdist2aexpn_spline != NULL)
	gsl_spline_free(cosmocalc_comvdist2aexpn_spline);
      cosmocalc_comvdist2aexpn_spline = gsl_spline_alloc(GSL_SPLINE_TYPE,(size_t) (COSMOCALC_COMVDIST_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_comvdist2aexpn_spline,comvdist_table,aexpn_table,(size_t) (COSMOCALC_COMVDIST_TABLE_LENGTH));
      if(cosmocalc_comvdist2aexpn_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_comvdist2aexpn_acc);
      else
	cosmocalc_comvdist2aexpn_acc = gsl_interp_accel_alloc();
    }
}

double acomvdist(double dist)
{
  static int initFlag = 1;
  static int currCosmoNum;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      init_cosmocalc_comvdist_table();
    }
  
  return gsl_spline_eval(cosmocalc_comvdist2aexpn_spline,dist,cosmocalc_comvdist2aexpn_acc);
}

double comvdist(double a)
{
  static int initFlag = 1;
  static int currCosmoNum;

  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      init_cosmocalc_comvdist_table();
    }
  
  return gsl_spline_eval(cosmocalc_aexpn2comvdist_spline,a,cosmocalc_aexpn2comvdist_acc);
}

//NOTE: you must call comvdist *FIRST* to init DH_sqrtok!
double angdist(double a)
{
  double cd;
  cd = comvdist(a);

  if(cosmoData.OmegaK > 0.0)
    return DH_sqrtok*sinh(cd/DH_sqrtok)*a;
  else if(cosmoData.OmegaK < 0)
    return DH_sqrtok*sin(cd/DH_sqrtok)*a;
  else
    return cd*a;
}

double lumdist(double a)
{
  double cd;
  cd = comvdist(a);
  
  if(cosmoData.OmegaK > 0.0)
    return DH_sqrtok*sinh(cd/DH_sqrtok)/a;
  else if(cosmoData.OmegaK < 0)
    return DH_sqrtok*sin(cd/DH_sqrtok)/a;
  else
    return cd/a;
}

double angdistdiff(double amin, double amax)
{
  double cdmin,cdmax;
  assert(amin <= amax);
  cdmin = comvdist(amin);
  cdmax = comvdist(amax);
  
  if(cosmoData.OmegaK > 0.0)
    return DH_sqrtok*sinh((cdmin-cdmax)/DH_sqrtok)*amin;
  else if(cosmoData.OmegaK < 0)
    return DH_sqrtok*sin((cdmin-cdmax)/DH_sqrtok)*amin;
  else
    return (cdmin-cdmax)*amin;
}

