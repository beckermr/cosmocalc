#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "cosmocalc.h"

#define ODEGROWTH //the integral formulation is onyl correct for LCDM, any wCDM cosmology will not work

#ifdef ODEGROWTH

#define COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH 50
#define AEXPN_MIN (1.0/31.0)
#define AEXPN_MAX 1.0

#define LNAINITGROWTH (log(AEXPN_MIN))
#define ABSERR 0.0
#define RELERR 1e-8
#define HSTART 0.0001

static int growth_ode_sys_w0wa(double t, const double y[], double dydt[], void *params);

//ODE for growth from Komatsu et al. 2009, ApJS, 180, 330 (WMAP5 paper)
static int growth_ode_sys_w0wa(double t, const double y[], double dydt[], void *params)
{
  double w0 = cosmoData.w0;
  double wa = cosmoData.wa;
  double a = exp(t); //t = log(a)
  
  double weffa;
  if(t != 0.0)
    weffa = w0 + wa - wa*(a - 1.0)/t;
  else
    weffa = w0;
  
  double ha = hubble_noscale(a);
  double omegaDEa = (1.0-cosmoData.OmegaM)/ha/ha/pow(a,3.0*(1.0 + weffa));
  
  dydt[0] = y[1];
  dydt[1] = (1.5*weffa*omegaDEa - 2.5)*y[1] - 1.5*(1.0 - weffa)*omegaDEa*y[0];
  
  return GSL_SUCCESS;
}

double growth_function_exact(double a)
{
  gsl_odeiv2_system odesys = {growth_ode_sys_w0wa, NULL, 2, NULL};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&odesys,gsl_odeiv2_step_rk8pd,HSTART,ABSERR,RELERR);
  
  double lna_init;
  double lna_final,y[2];
  double ga,gnorm;
  int status;
  
  //do g(a) unormalized
  y[0] = 1.0;
  y[1] = 0.0;
  lna_init = LNAINITGROWTH;
  lna_final = log(a);
  
  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
  ga = y[0]*a;
  if(status != GSL_SUCCESS)
    {
      fprintf(stderr,"error in integrating growth function! a = %lf\n",a);
      assert(status == GSL_SUCCESS);
    }

  //do g(a = 1.0) to get normalization
  y[0] = 1.0;
  y[1] = 0.0;
  lna_init = LNAINITGROWTH;
  lna_final = 0.0;
  
  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
  gnorm = y[0];
  if(status != GSL_SUCCESS)
    {
      fprintf(stderr,"error in integrating growth function! a = %lf\n",a);
      assert(status == GSL_SUCCESS);
    }
    
  gsl_odeiv2_driver_free(d);
  
  return ga/gnorm;
}

double growth_function(double a)
{
  static int initFlag = 1;
  static int currCosmoNum;
  static double growth_function_norm;
  static gsl_spline *cosmocalc_growth_function_spline = NULL;
  static gsl_interp_accel *cosmocalc_growth_function_acc = NULL;
  
  double growth_function_table[COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH];
  double a_table[COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH];
  long i;
    
  gsl_odeiv2_system odesys = {growth_ode_sys_w0wa, NULL, 2, NULL};
  gsl_odeiv2_driver *d;
  double afact,lna_init;
  double lna_final,y[2];
  int status;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      //init ODE integrator
      d = gsl_odeiv2_driver_alloc_y_new(&odesys,gsl_odeiv2_step_rk8pd,HSTART,ABSERR,RELERR);
      
      //init ODE
      y[0] = 1.0;
      y[1] = 0.0;
      lna_init = LNAINITGROWTH;
      
      for(i=0;i<COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH;++i)
        {
	  //get scale factor
	  afact = (AEXPN_MAX - AEXPN_MIN)/(COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH-1.0)*((double) i) + AEXPN_MIN;
	  a_table[i] = afact;
	  
	  //do growth function integration
	  lna_final = log(afact);
  	  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
	  growth_function_table[i] = y[0]*afact;
	  if(status != GSL_SUCCESS)
	    {
	      fprintf(stderr,"error in integrating growth function! a = %lf\n",a);
	      assert(status == GSL_SUCCESS);
	    }
	  
	  //step to next interval
          lna_init = lna_final;
	}
      
      //get normalization
      growth_function_norm = growth_function_table[COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH-1];
            
      //free ODE stuff
      gsl_odeiv2_driver_free(d);

      //init the spline and accelerators
      if(cosmocalc_growth_function_spline != NULL)
        gsl_spline_free(cosmocalc_growth_function_spline);
      cosmocalc_growth_function_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_growth_function_spline,a_table,growth_function_table,(size_t) (COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH));
      if(cosmocalc_growth_function_acc != NULL)
        gsl_interp_accel_reset(cosmocalc_growth_function_acc);
      else
        cosmocalc_growth_function_acc = gsl_interp_accel_alloc();
    }

  return gsl_spline_eval(cosmocalc_growth_function_spline,a,cosmocalc_growth_function_acc)/growth_function_norm;
}

#else

#define COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH 25
#define AEXPN_MIN 0.001
#define AEXPN_MAX 1.0

static double growth_function_integ_funct(double a, void *p);

/* function for integration using gsl integration */
static double growth_function_integ_funct(double a, void *p)
{
  double hubblefac = hubble_noscale(a);
  double adot = a*hubblefac;
  double alim = (*((double*)p));
  
  return 1.0/adot/adot/adot*hubble_noscale(alim);
}

double growth_function_exact(double a)
{
#define WORKSPACE_NUM 100000
#define ABSERR 0.0
#define RELERR 1e-8

  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr,afact,norm;
  
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  
  afact = a;
  F.function = &growth_function_integ_funct;
  F.params = &(afact);
  gsl_integration_qag(&F,0.0,afact,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
  
  afact = 1.0;
  F.params = &(afact);
  gsl_integration_qag(&F,0.0,afact,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&norm,&abserr);
  
  gsl_integration_workspace_free(workspace);

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return result/norm;
}

double growth_function(double a)
{
#define WORKSPACE_NUM 100000
#define ABSERR 1e-8
#define RELERR 1e-8
  
  static int initFlag = 1;
  static int currCosmoNum;
  static double growth_function_norm;
  static gsl_spline *cosmocalc_growth_function_spline = NULL;
  static gsl_interp_accel *cosmocalc_growth_function_acc = NULL;
  
  double growth_function_table[COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH];
  double a_table[COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH];
  long i;
  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr,afact;
  
  if(initFlag == 1 || currCosmoNum != cosmoData.cosmoNum)
    {
      initFlag = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
      F.function = &growth_function_integ_funct;
            
      for(i=0;i<COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH;++i)
	{
	  afact = (AEXPN_MAX - AEXPN_MIN)/(COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH-1.0)*((double) i) + AEXPN_MIN;
	  a_table[i] = afact;
	  F.params = &(afact);
	  gsl_integration_qag(&F,0.0,afact,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
	  growth_function_table[i] = result;
	}
      
      afact = 1.0;
      F.params = &(afact);
      gsl_integration_qag(&F,0.0,1.0,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&growth_function_norm,&abserr);
      
      gsl_integration_workspace_free(workspace);
      
      //init the spline and accelerators
      if(cosmocalc_growth_function_spline != NULL)
	gsl_spline_free(cosmocalc_growth_function_spline);
      cosmocalc_growth_function_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH));
      gsl_spline_init(cosmocalc_growth_function_spline,a_table,growth_function_table,(size_t) (COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH));
      if(cosmocalc_growth_function_acc != NULL)
	gsl_interp_accel_reset(cosmocalc_growth_function_acc);
      else
	cosmocalc_growth_function_acc = gsl_interp_accel_alloc();
    }

  return gsl_spline_eval(cosmocalc_growth_function_spline,a,cosmocalc_growth_function_acc)/growth_function_norm;

#undef WORKSPACE_NUM
#undef ABSERR
#undef RELERR
}
#endif

#undef ODEGROWTH
#undef COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH
#undef AEXPN_MIN
#undef AEXPN_MAX
