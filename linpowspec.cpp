#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_fit.h>

#include "linpowspec.h"

//#define TIMELPS
#ifdef TIMELPS
#include <sys/time.h>
static double wtime(void)
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ((double) (tp.tv_sec)) + ((double) (tp.tv_usec))/1e6;
}
#endif

static double fourierTransformTopHat(double y);
static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p);
static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p);

static double fourierTransformTopHat(double y)
{
  if(y < 1e-3)
    return 1.0;
  else
    return 3.0/y/y/y*(sin(y) - y*cos(y));
}

struct linear_powspec_params {
  double topHatRad;
  double ns;
  class TransferFunction *tf;
  class LinearPowerSpectrum *lp;
};

static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double k = exp(lnk);
  double ft = fourierTransformTopHat(k*lpp->topHatRad);
  
  //double tf = (*(lpp->tf))(k);
  //return ft*ft*tf*tf*pow(k,lpp->ns)*k*k*k/2.0/M_PI/M_PI;
  
  return ft*ft*(*(lpp->lp))(k,1.0)*k*k*k/2.0/M_PI/M_PI;
}

static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double ft = fourierTransformTopHat(k*lpp->topHatRad);

  //double tf = (*(lpp->tf))(k);
  //return ft*ft*tf*tf*pow(k,lpp->ns)*k*k/2.0/M_PI/M_PI;
  
  return ft*ft*(*(lpp->lp))(k,1.0)*k*k/2.0/M_PI/M_PI;
}

double LinearPowerSpectrum::tophatradnorm_linear_powspec_exact_nonorm(double topHatRad)
{
  double I0,I1;
  double abserr;
  double epsrel,epsabs;
  gsl_integration_workspace *workspace;
  gsl_function F;
  struct linear_powspec_params lpp;

  lpp.topHatRad = topHatRad;
  lpp.ns = _ns;
  lpp.lp = this;
  lpp.tf = _tf;
  
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-8
#define RELERR 0.0 
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  cosmocalc_assert(workspace != NULL,"could not alloc workspace for GSL integration for top hat norm of linear powspec!");
  
  F.params = &lpp;
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

void LinearPowerSpectrum::init_linear_powspec_table(void)
{
#ifdef TIMELPS
  double tt,tn;
  tt = - wtime();
#endif
  
  double *ln_linear_powspec_table = (double*)malloc(sizeof(double)*LINEAR_POWSPEC_TABLE_LENGTH);
  double *lnk_table = (double*)malloc(sizeof(double)*LINEAR_POWSPEC_TABLE_LENGTH);
  long i;
  double cov00,cov01,cov11,sumsq;
  double lntf,dlnk,lnkmin;
    
  cosmocalc_assert(ln_linear_powspec_table != NULL,"out of memory for linear powspec table!");
  cosmocalc_assert(lnk_table != NULL,"out of memory for linear powspec");
  
  dlnk = log(PL_K_MAX/PL_K_MIN)/(LINEAR_POWSPEC_TABLE_LENGTH-1.0);
  lnkmin = log(PL_K_MIN);

  //#ifdef _OPENMP
  //if(_num_threads > 0) omp_set_num_threads(_num_threads);
  //#endif

#pragma omp parallel for default(none) private(lntf,i) shared(lnk_table,ln_linear_powspec_table,dlnk,lnkmin)
  for(i=0;i<LINEAR_POWSPEC_TABLE_LENGTH;++i)
    {
      lnk_table[i] = dlnk*i + lnkmin;
      lntf = log((*_tf)(exp(lnk_table[i])));
      ln_linear_powspec_table[i] = 2.0*lntf + _ns*lnk_table[i];
    }
            
  //init the spline and accelerators
  if(linear_powspec_spline != NULL)
    gsl_spline_free(linear_powspec_spline);
  linear_powspec_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (LINEAR_POWSPEC_TABLE_LENGTH));
  cosmocalc_assert(linear_powspec_spline != NULL,"could not alloc spline for linear powspec table!");
  gsl_spline_init(linear_powspec_spline,lnk_table,ln_linear_powspec_table,(size_t) (LINEAR_POWSPEC_TABLE_LENGTH));
  if(linear_powspec_acc != NULL)
    gsl_interp_accel_reset(linear_powspec_acc);
  else
    {
      linear_powspec_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(linear_powspec_acc != NULL,"could not alloc accel for linear powspec!");
    }
  
  gsl_fit_linear(lnk_table+LINEAR_POWSPEC_TABLE_LENGTH-LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		 ln_linear_powspec_table+LINEAR_POWSPEC_TABLE_LENGTH-LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		 (size_t) LINEAR_POWSPEC_FIT_LENGTH,&_linear_powspec_c0,&_linear_powspec_c1,&cov00,&cov01,&cov11,&sumsq);
  
  //get norm
#ifdef TIMELPS
  tt += wtime();
  tn = -wtime();
#endif
  _linear_powspec_norm = 1.0;
  _linear_powspec_norm = _s8*_s8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
  
  /*double lnnorm = log(_linear_powspec_norm);
  for(i=0;i<LINEAR_POWSPEC_TABLE_LENGTH;++i)
    ln_linear_powspec_table[i] += lnnorm;
  
  //init the spline and accelerators
  if(linear_powspec_spline != NULL)
    gsl_spline_free(linear_powspec_spline);
  linear_powspec_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (LINEAR_POWSPEC_TABLE_LENGTH));
  cosmocalc_assert(linear_powspec_spline != NULL,"could not alloc spline for linear powspec table!");
  gsl_spline_init(linear_powspec_spline,lnk_table,ln_linear_powspec_table,(size_t) (LINEAR_POWSPEC_TABLE_LENGTH));
  if(linear_powspec_acc != NULL)
    gsl_interp_accel_reset(linear_powspec_acc);
  else
    {
      linear_powspec_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(linear_powspec_acc != NULL,"could not alloc accel for linear powspec!");
    }
  
  gsl_fit_linear(lnk_table+LINEAR_POWSPEC_TABLE_LENGTH-LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		 ln_linear_powspec_table+LINEAR_POWSPEC_TABLE_LENGTH-LINEAR_POWSPEC_FIT_LENGTH,(size_t) 1,
		 (size_t) LINEAR_POWSPEC_FIT_LENGTH,&_linear_powspec_c0,&_linear_powspec_c1,&cov00,&cov01,&cov11,&sumsq);
  
  */
  free(ln_linear_powspec_table);
  free(lnk_table);
#ifdef TIMELPS
  tn += wtime();
  fprintf(stderr,"time table,norm = %f|%f seconds.\n",tt,tn);
#endif
}

