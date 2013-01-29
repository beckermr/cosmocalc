#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "w0wacosmo.h"

#define ABSERR 0.0
#define RELERR 1e-8
#define HSTART 0.0001

struct growth_params {
  class Hubble *h;
  class w0wa_CosmoData *cd;
};

//ODE for growth from Komatsu et al. 2009, ApJS, 180, 330 (WMAP5 paper)
static int growth_ode_sys_w0wa(double t, const double y[], double dydt[], void *params)
{
  class growth_params *gp = (struct growth_params*)params;
    
  double w0 = gp->cd->w0;
  double wa = gp->cd->wa;
  double a = exp(t); //t = log(a)
  
  double weffa;
  if(t != 0.0)
    weffa = w0 + wa - wa*(a - 1.0)/t;
  else
    weffa = w0;
  
  double ha = (*(gp->h))(a);
  double omegaDEa = gp->cd->ol/ha/ha/pow(a,3.0*(1.0 + weffa));
  double omegaKa = gp->cd->ok/a/a/ha/ha;
  
  dydt[0] = y[1];
  dydt[1] = (1.5*weffa*omegaDEa - 0.5*omegaKa - 2.5)*y[1] - (2.0*omegaKa + 1.5*(1.0 - weffa)*omegaDEa)*y[0];
  
  return GSL_SUCCESS;
}

double w0wa_GrowthFunction::growth_function_exact(double k, double a, class Hubble& h)
{
  class w0wa_CosmoData cd = this->cosmology();
  struct growth_params gp = {&h, &cd};
  
  gsl_odeiv2_system odesys = {growth_ode_sys_w0wa, NULL, 2, &gp};
  gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(&odesys,gsl_odeiv2_step_rk8pd,HSTART,ABSERR,RELERR);
  cosmocalc_assert(d != NULL,"could not alloc GSL ODE driver for exact growth function computation!");
  
  double lna_init;
  double lna_final,y[2];
  double ga,gnorm;
  int status;
  
  //do g(a) unormalized
  y[0] = 1.0;
  y[1] = 0.0;
  lna_init = LOG_AEXPN_MIN;
  lna_final = log(a);
  
  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
  ga = y[0]*a;
  cosmocalc_assert(status == GSL_SUCCESS,"error in integrating growth function to a = %lf for exact growth function computation!",a);
  
  //do g(a = 1.0) to get normalization
  y[0] = 1.0;
  y[1] = 0.0;
  lna_init = LOG_AEXPN_MIN;
  lna_final = 0.0;
  
  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
  gnorm = y[0];
  cosmocalc_assert(status == GSL_SUCCESS,"error getting growth function norm for exact growth function computation!");
  
  gsl_odeiv2_driver_free(d);
  
  return ga/gnorm;
}

void w0wa_GrowthFunction::init_growth_function_table(class Hubble &h)
{
  double *growth_function_table = (double*)malloc(sizeof(double)*GROWTH_FUNCTION_TABLE_LENGTH);
  double *a_table = (double*)malloc(sizeof(double)*GROWTH_FUNCTION_TABLE_LENGTH);
  long i;
  
  cosmocalc_assert(growth_function_table != NULL,"out of memory for growth function table!");
  cosmocalc_assert(a_table != NULL,"out of memory for growth function table!");
  
  class w0wa_CosmoData cd = this->cosmology();
  struct growth_params gp = {&h, &cd};
  
  gsl_odeiv2_system odesys = {growth_ode_sys_w0wa, NULL, 2, &gp};
  gsl_odeiv2_driver *d;
  double afact,lna_init;
  double lna_final,y[2];
  int status;
  double da,amin;
  double amax = AEXPN_MAX+1e-3;
  
  da = (amax - AEXPN_MIN)/(GROWTH_FUNCTION_TABLE_LENGTH-1.0);
  amin = AEXPN_MIN;
  
  //get normalization
  d = gsl_odeiv2_driver_alloc_y_new(&odesys,gsl_odeiv2_step_rk8pd,HSTART,ABSERR,RELERR);
  cosmocalc_assert(d != NULL,"could not alloc GSL ODE driver for growth function table!");

  y[0] = 1.0;
  y[1] = 0.0;
  lna_init = LOG_AEXPN_MIN;
  lna_final = 0.0;
  
  status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
  _growth_function_norm = y[0];
  cosmocalc_assert(status == GSL_SUCCESS,"error getting growth function norm for growth function table!");
  
  gsl_odeiv2_driver_free(d);
  
  //#ifdef _OPENMP
  //if(_num_threads > 0) omp_set_num_threads(_num_threads);
  //#endif

#pragma omp parallel default(none) \
  private(d,y,lna_init,status,lna_final,i,afact)	\
  shared(odesys,a_table,growth_function_table,gsl_odeiv2_step_rk8pd,da,amin,stderr)
  {
    //init ODE integrator
    d = gsl_odeiv2_driver_alloc_y_new(&odesys,gsl_odeiv2_step_rk8pd,HSTART,ABSERR,RELERR);
    cosmocalc_assert(d != NULL,"could not alloc GSL ODE driver for growth function table!");

#pragma omp for
    for(i=0;i<GROWTH_FUNCTION_TABLE_LENGTH;++i)
      {
	//get scale factor
	afact = da*i+amin;
	a_table[i] = afact;
	
	//do growth function integration
	y[0] = 1.0;
	y[1] = 0.0;
	lna_init = LOG_AEXPN_MIN;
	lna_final = log(afact);
	status = gsl_odeiv2_driver_apply(d,&lna_init,lna_final,y);
	growth_function_table[i] = y[0]*afact/_growth_function_norm;
	cosmocalc_assert(status == GSL_SUCCESS,"error in integrating growth function to a = %lf for growth function table!",afact);
      }
    
    //free ODE stuff
    gsl_odeiv2_driver_free(d);
  }
  
  //init the spline and accelerators
  if(growth_function_spline != NULL)
    gsl_spline_free(growth_function_spline);
  growth_function_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (GROWTH_FUNCTION_TABLE_LENGTH));
  cosmocalc_assert(growth_function_spline != NULL,"could not alloc spline for growth function table!");
  gsl_spline_init(growth_function_spline,a_table,growth_function_table,(size_t) (GROWTH_FUNCTION_TABLE_LENGTH));
  if(growth_function_acc != NULL)
    gsl_interp_accel_reset(growth_function_acc);
  else
    {
      growth_function_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(growth_function_acc != NULL,"could not alloc accel for growth function table!");
    }
  
  //free mem
  free(a_table);
  free(growth_function_table);
}

#undef ABSERR
#undef RELERR
#undef HSTART
