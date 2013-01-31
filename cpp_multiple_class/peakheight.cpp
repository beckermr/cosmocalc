#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "peakheight.h"

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
  class PowerSpectrum *pkl;
};

static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double k = exp(lnk);
  double ft = fourierTransformTopHat(k*lpp->topHatRad);
  
  return ft*ft*(*(lpp->pkl))(k,1.0)*k*k*k/2.0/M_PI/M_PI;
}

static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double ft = fourierTransformTopHat(k*lpp->topHatRad);

  return ft*ft*(*(lpp->pkl))(k,1.0)*k*k/2.0/M_PI/M_PI;
}

void PeakHeight::init_peakheight_table(void)
{
  double *sigmar_table = (double*)malloc(sizeof(double)*PEAKHEIGHT_TABLE_LENGTH);
  double *topHatRad_table = (double*)malloc(sizeof(double)*PEAKHEIGHT_TABLE_LENGTH);
  double *tmpDouble = (double*)malloc(sizeof(double)*PEAKHEIGHT_TABLE_LENGTH);
  long i;
  double lnR,dlnR,lnRmin,topHatRad;
  
  cosmocalc_assert(sigmar_table != NULL,"out of memory for peak height tables!");
  cosmocalc_assert(topHatRad_table != NULL,"out of memory for peak height tables!");
  cosmocalc_assert(tmpDouble != NULL,"out of memory for peak height tables!");
    
  double I0,I1;
  double abserr;
  double epsrel,epsabs;
  gsl_integration_workspace *workspace;
  gsl_function F;
  struct linear_powspec_params lpp;

#define WORKSPACE_NUM 10000000
#define ABSERR 1e-8
#define RELERR 0.0

  //#ifdef _OPENMP
  //if(_num_threads > 0) omp_set_num_threads(_num_threads);
  //#endif

  dlnR = log(PH_R_MAX/PH_R_MIN)/(PEAKHEIGHT_TABLE_LENGTH-1.0);
  lnRmin = log(PH_R_MIN);
  
#pragma omp parallel default(none) \
  private(workspace,lpp,F,epsabs,epsrel,abserr,I0,I1,topHatRad,i,lnR)	\
  shared(sigmar_table,topHatRad_table,dlnR,lnRmin)
  {
    workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
    cosmocalc_assert(workspace != NULL,"could not alloc workspace for GSL integration for peak height table!");
    
    lpp.pkl = _pkl;
    F.params = &lpp;

#pragma omp for
    for(i=0;i<PEAKHEIGHT_TABLE_LENGTH;++i)
      {
	lnR = dlnR*i + lnRmin;
	topHatRad = exp(lnR);
	lpp.topHatRad = topHatRad;

	if(topHatRad > 1e-4)
	  {
	    epsabs = 1e-15;
	    epsrel = 1e-5;
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
	
	topHatRad_table[i] = lnR;
	sigmar_table[i] = log(sqrt((I0+I1)));
      }

    gsl_integration_workspace_free(workspace);
  }

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  //init the spline and accelerators
  if(R2sigma_spline != NULL)
    gsl_spline_free(R2sigma_spline);
  R2sigma_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (PEAKHEIGHT_TABLE_LENGTH));
  cosmocalc_assert(R2sigma_spline != NULL,"could not alloc spline for sigma(R) for peak height!");
  gsl_spline_init(R2sigma_spline,topHatRad_table,sigmar_table,(size_t) (PEAKHEIGHT_TABLE_LENGTH));
  if(R2sigma_acc != NULL)
    gsl_interp_accel_reset(R2sigma_acc);
  else
    {
      R2sigma_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(R2sigma_acc!= NULL,"could not alloc accel for sigma(R) for peak height!");
    }
  
  //sort the peak heights properly
  for(i=0;i<PEAKHEIGHT_TABLE_LENGTH;++i)
    tmpDouble[i] = sigmar_table[PEAKHEIGHT_TABLE_LENGTH-1-i];
  for(i=0;i<PEAKHEIGHT_TABLE_LENGTH;++i)
    sigmar_table[i] = tmpDouble[i];
  for(i=0;i<PEAKHEIGHT_TABLE_LENGTH;++i)
    tmpDouble[i] = topHatRad_table[PEAKHEIGHT_TABLE_LENGTH-1-i];
  for(i=0;i<PEAKHEIGHT_TABLE_LENGTH;++i)
    topHatRad_table[i] = tmpDouble[i];
  
  //init the spline and accelerators
  if(sigma2R_spline != NULL)
    gsl_spline_free(sigma2R_spline);
  sigma2R_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (PEAKHEIGHT_TABLE_LENGTH));
  cosmocalc_assert(sigma2R_spline != NULL,"could not alloc spline for R(sigma) for peak height!");
  gsl_spline_init(sigma2R_spline,sigmar_table,topHatRad_table,(size_t) (PEAKHEIGHT_TABLE_LENGTH));
  if(sigma2R_acc != NULL)
    gsl_interp_accel_reset(sigma2R_acc);
  else
    {
      sigma2R_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(sigma2R_acc!= NULL,"could not alloc accel for R(sigma) for peak height!");
    }
  
  free(tmpDouble);
  free(topHatRad_table);
  free(sigmar_table);
}

