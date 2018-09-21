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

static double linear_corrfunc_integ_funct(double k, void *p);

static double linear_corrfunc_integ_funct(double k, void *p)
{
  double r = ((double*)p)[0];
  double a = ((double*)p)[1];
  return linear_powspec(k,a)*k/r;
}

double linear_corrfunc_exact(double r, double a)
{
  double I0;
  double abserr,p[2];
  gsl_integration_workspace *workspace;
  gsl_function F;
  gsl_integration_qawo_table *wf;
  gsl_integration_workspace *cycle_workspace;
  
#define WORKSPACE_NUM 70
#define ABSERR 1e-3
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  cycle_workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  wf = gsl_integration_qawo_table_alloc(r,1.0,GSL_INTEG_SINE,(size_t) WORKSPACE_NUM);
  
  F.function = &linear_corrfunc_integ_funct;
  p[0] = r;
  p[1] = a;
  F.params = p;
  gsl_integration_qawf(&F,1e-8,ABSERR,(size_t) WORKSPACE_NUM,workspace,cycle_workspace,wf,&I0,&abserr);
  
  gsl_integration_qawo_table_free(wf);
  gsl_integration_workspace_free(cycle_workspace);
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef WORKSPACE_NUM
  
  return I0/2.0/M_PI/M_PI;
}

double linear_corrfunc(double r, double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static gsl_spline *cosmocalc_linear_corrfunc_spline = NULL;
  static gsl_interp_accel *cosmocalc_linear_corrfunc_acc = NULL; 
  
  double linear_corrfunc_table[COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH];
  double r_table[COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH];
  long i;
  double gf;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      for(i=0;i<COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH;++i)
	{
	  r_table[i] = log(CF_R_MAX/CF_R_MIN)/(COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH-1.0)*((double) i) + log(CF_R_MIN);
	  linear_corrfunc_table[i] = linear_corrfunc_exact(exp(r_table[i]),1.0);
	}
            
      //init the spline and accelerators
      if(cosmocalc_linear_corrfunc_spline != NULL)
	gsl_spline_free(cosmocalc_linear_corrfunc_spline);
      cosmocalc_linear_corrfunc_spline = gsl_spline_alloc(GSL_SPLINE_TYPE,(size_t) (COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_linear_corrfunc_spline,r_table,linear_corrfunc_table,(size_t) (COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH));
      if(cosmocalc_linear_corrfunc_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_linear_corrfunc_acc);
      else
	cosmocalc_linear_corrfunc_acc = gsl_interp_accel_alloc();
    }
  
  gf = growth_function(a);
  
  if(r < CF_R_MIN)
    return 0.0;
  else if(r < CF_R_MAX)
    return gsl_spline_eval(cosmocalc_linear_corrfunc_spline,log(r),cosmocalc_linear_corrfunc_acc)*gf*gf;
  else
    return 0.0;
}

