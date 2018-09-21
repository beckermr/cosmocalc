#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "cosmocalc.h"

gsl_spline *cosmocalc_aexpn2age_spline = NULL;
gsl_interp_accel *cosmocalc_aexpn2age_acc = NULL; 
double AGE = 0;

/* function for integration using gsl integration */
static double age_integ_funct(double a, void *p)
{
  return 1.0/a/hubble_noscale(a);
}

/* init function  - some help from Gadget-2 applied here */
void init_cosmocalc_age_table(void)
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
  double age_table[COSMOCALC_AGE_TABLE_LENGTH];
  double aexpn_table[COSMOCALC_AGE_TABLE_LENGTH];
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);

      aexpn_table[0] = 0.0;
      age_table[0] = 0.0;  
      for(i=1;i<COSMOCALC_AGE_TABLE_LENGTH;++i)
	{
	  afact = (1.0 - 0.0)/(COSMOCALC_AGE_TABLE_LENGTH-1.0)*((double) i);
	  F.function = &age_integ_funct;
	  gsl_integration_qag(&F,1e-12,afact,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
	  aexpn_table[i] = afact;
	  age_table[i] = result*TH/cosmoData.h;
	}
      
      gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
      
      //init the spline and accelerators
      if(cosmocalc_aexpn2age_spline != NULL)
	gsl_spline_free(cosmocalc_aexpn2age_spline);
      cosmocalc_aexpn2age_spline = gsl_spline_alloc(GSL_SPLINE_TYPE,(size_t) (COSMOCALC_AGE_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_aexpn2age_spline,aexpn_table,age_table,(size_t) (COSMOCALC_AGE_TABLE_LENGTH));
      if(cosmocalc_aexpn2age_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_aexpn2age_acc);
      else
	cosmocalc_aexpn2age_acc = gsl_interp_accel_alloc();
      
      AGE = gsl_spline_eval(cosmocalc_aexpn2age_spline,1.0,cosmocalc_aexpn2age_acc);
    }
}

double age(double a)
{
  static int initFlag = 1;
  static int currCosmoNum;

  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      init_cosmocalc_age_table();
    }
  
  return gsl_spline_eval(cosmocalc_aexpn2age_spline,a,cosmocalc_aexpn2age_acc);
}

//call age first to set AGE
double lookback(double a)
{
  double age_a = age(a);
  return AGE - age_a;
}
