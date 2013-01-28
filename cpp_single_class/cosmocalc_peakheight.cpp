#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include "cosmocalc.h"
#include "cosmocalc_assert.h"

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
  cosmoCalc *cd;
};

static double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double k = exp(lnk);
  double ft = fourierTransformTopHat(k*lpp->topHatRad);
  double tf = lpp->cd->transfer_function(k);
  
  return ft*ft*tf*tf*pow(k,lpp->cd->ns())*k*k*k/2.0/M_PI/M_PI;
}

static double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p)
{
  struct linear_powspec_params *lpp = (struct linear_powspec_params*)p;
  double ft = fourierTransformTopHat(k*lpp->topHatRad);
  double tf = lpp->cd->transfer_function(k);

  return ft*ft*tf*tf*pow(k,lpp->cd->ns())*k*k/2.0/M_PI/M_PI;
}

void cosmoCalc::init_cosmocalc_peakheight_table(void)
{
  
  double *sigmar_table = (double*)malloc(sizeof(double)*COSMOCALC_PEAKHEIGHT_TABLE_LENGTH);
  double *topHatRad_table = (double*)malloc(sizeof(double)*COSMOCALC_PEAKHEIGHT_TABLE_LENGTH);
  double *tmpDouble = (double*)malloc(sizeof(double)*COSMOCALC_PEAKHEIGHT_TABLE_LENGTH);
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

#ifdef _OPENMP
  if(_num_threads > 0) omp_set_num_threads(_num_threads);
#endif

  dlnR = log(PH_R_MAX/PH_R_MIN)/(COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1.0);
  lnRmin = log(PH_R_MIN);
  
#pragma omp parallel default(none) \
  private(workspace,lpp,F,epsabs,epsrel,abserr,I0,I1,topHatRad,i,lnR)	\
  shared(sigmar_table,topHatRad_table,dlnR,lnRmin)
  {
    workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
    cosmocalc_assert(workspace != NULL,"could not alloc workspace for GSL integration for peak height table!");
    
    lpp.cd = this;
    F.params = &lpp;

#pragma omp for
    for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
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
	sigmar_table[i] = log(sqrt(_linear_powspec_norm*(I0+I1)));
      }

    gsl_integration_workspace_free(workspace);
  }

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM

  //init the spline and accelerators
  if(cosmocalc_R2sigma_spline != NULL)
    gsl_spline_free(cosmocalc_R2sigma_spline);
  cosmocalc_R2sigma_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
  cosmocalc_assert(cosmocalc_R2sigma_spline != NULL,"could not alloc spline for sigma(R) for peak height!");
  gsl_spline_init(cosmocalc_R2sigma_spline,topHatRad_table,sigmar_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
  if(cosmocalc_R2sigma_acc != NULL)
    gsl_interp_accel_reset(cosmocalc_R2sigma_acc);
  else
    {
      cosmocalc_R2sigma_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(cosmocalc_R2sigma_acc!= NULL,"could not alloc accel for sigma(R) for peak height!");
    }
  
  //sort the peak heights properly
  for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
    tmpDouble[i] = sigmar_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1-i];
  for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
    sigmar_table[i] = tmpDouble[i];
  for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
    tmpDouble[i] = topHatRad_table[COSMOCALC_PEAKHEIGHT_TABLE_LENGTH-1-i];
  for(i=0;i<COSMOCALC_PEAKHEIGHT_TABLE_LENGTH;++i)
    topHatRad_table[i] = tmpDouble[i];
  
  //init the spline and accelerators
  if(cosmocalc_sigma2R_spline != NULL)
    gsl_spline_free(cosmocalc_sigma2R_spline);
  cosmocalc_sigma2R_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
  cosmocalc_assert(cosmocalc_sigma2R_spline != NULL,"could not alloc spline for R(sigma) for peak height!");
  gsl_spline_init(cosmocalc_sigma2R_spline,sigmar_table,topHatRad_table,(size_t) (COSMOCALC_PEAKHEIGHT_TABLE_LENGTH));
  if(cosmocalc_sigma2R_acc != NULL)
    gsl_interp_accel_reset(cosmocalc_sigma2R_acc);
  else
    {
      cosmocalc_sigma2R_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(cosmocalc_sigma2R_acc!= NULL,"could not alloc accel for R(sigma) for peak height!");
    }
  
  free(tmpDouble);
  free(topHatRad_table);
  free(sigmar_table);
}

