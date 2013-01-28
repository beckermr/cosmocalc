#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "w0wacosmo.h"

#define WORKSPACE_NUM 100000
#define ABSERR 0.0
#define RELERR 1e-8

/* function for integration using gsl integration */
static double comvdist_integ_funct(double a, void *p)
{
  class Hubble *h = (class Hubble*)p;
  return 1.0/a/a/(*h)(a);
}

/* init function  - some help from Gadget-2 applied here */
void w0wa_Distances::init_comvdist_table(class Hubble& h)
{
  gsl_integration_workspace *workspace;
  gsl_function F;
  long i;
  double result,abserr,afact;
  double *comvdist_table = (double*)malloc(sizeof(double)*COMVDIST_TABLE_LENGTH);
  double *aexpn_table = (double*)malloc(sizeof(double)*COMVDIST_TABLE_LENGTH);
  double *tmpDouble = (double*)malloc(sizeof(double)*COMVDIST_TABLE_LENGTH);
  double da,amin;
  
  cosmocalc_assert(comvdist_table != NULL,"out of memory for distances table!");
  cosmocalc_assert(aexpn_table != NULL,"out of memory for distances table!");
  cosmocalc_assert(tmpDouble != NULL,"out of memory for distances table!");
  
  F.function = &comvdist_integ_funct;
  F.params = (void*)(&h);
  amin = AEXPN_MIN;
  da = (AEXPN_MAX - AEXPN_MIN)/(COMVDIST_TABLE_LENGTH-1.0);

  //#ifdef _OPENMP
  //if(_num_threads > 0) omp_set_num_threads(_num_threads);
  //#endif
  
#pragma omp parallel default(none) \
  private(i,workspace,result,abserr,afact)		\
  shared(amin,da,F,aexpn_table,comvdist_table,stderr)
  {  
    workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
    cosmocalc_assert(workspace != NULL,"could not alloc GSL integration workspace for distances table!");
    
#pragma omp for 
    for(i=0;i<COMVDIST_TABLE_LENGTH-1;++i)
      {
	afact = da*i + amin;
	gsl_integration_qag(&F,afact,1.0,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
	aexpn_table[i] = afact;
	comvdist_table[i] = result;
      }
    
    gsl_integration_workspace_free(workspace);
  }
  aexpn_table[COMVDIST_TABLE_LENGTH-1] = 1.0;
  comvdist_table[COMVDIST_TABLE_LENGTH-1] = 0.0;
  for(i=0;i<COMVDIST_TABLE_LENGTH-1;++i)
    comvdist_table[i] *= DH;
  
  //init the spline and accelerators
  if(aexpn2comvdist_spline != NULL)
    gsl_spline_free(aexpn2comvdist_spline);
  aexpn2comvdist_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COMVDIST_TABLE_LENGTH));
  cosmocalc_assert(aexpn2comvdist_spline != NULL,"could not alloc a->comvdist spline for distances table!");
  gsl_spline_init(aexpn2comvdist_spline,aexpn_table,comvdist_table,(size_t) (COMVDIST_TABLE_LENGTH));
  if(aexpn2comvdist_acc != NULL)
    gsl_interp_accel_reset(aexpn2comvdist_acc);
  else
    {
      aexpn2comvdist_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(aexpn2comvdist_acc != NULL,"could not alloc a->comvdist accel for distances table!");
    }
      
  //sort the distances properly
  for(i=0;i<COMVDIST_TABLE_LENGTH;++i)
    tmpDouble[i] = comvdist_table[COMVDIST_TABLE_LENGTH-1-i];
  for(i=0;i<COMVDIST_TABLE_LENGTH;++i)
    comvdist_table[i] = tmpDouble[i];
  for(i=0;i<COMVDIST_TABLE_LENGTH;++i)
    tmpDouble[i] = aexpn_table[COMVDIST_TABLE_LENGTH-1-i];
  for(i=0;i<COMVDIST_TABLE_LENGTH;++i)
    aexpn_table[i] = tmpDouble[i];
  
  //init other splines and accelerators
  if(comvdist2aexpn_spline != NULL)
    gsl_spline_free(comvdist2aexpn_spline);
  comvdist2aexpn_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COMVDIST_TABLE_LENGTH));
  cosmocalc_assert(comvdist2aexpn_spline != NULL,"could not alloc comvdist->a spline for distances table!");
  gsl_spline_init(comvdist2aexpn_spline,comvdist_table,aexpn_table,(size_t) (COMVDIST_TABLE_LENGTH));
  if(comvdist2aexpn_acc != NULL)
    gsl_interp_accel_reset(comvdist2aexpn_acc);
  else
    {
      comvdist2aexpn_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(comvdist2aexpn_acc != NULL,"could not alloc comvdist->a accel for distances table!");
    }
  
  //free mem
  free(comvdist_table);
  free(aexpn_table);
  free(tmpDouble);
}

double w0wa_Distances::comvdist_exact(double a, class Hubble& h)
{
  gsl_integration_workspace *workspace;
  gsl_function F;
  double result,abserr;
  
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  cosmocalc_assert(workspace != NULL,"could not alloc GSL integration workspace for exact comvdist computation!");
  
  F.function = &comvdist_integ_funct;
  F.params = (void*)(&h);
  gsl_integration_qag(&F,a,1.0,ABSERR,RELERR,(size_t) WORKSPACE_NUM,GSL_INTEG_GAUSS51,workspace,&result,&abserr);
  
  gsl_integration_workspace_free(workspace);
  
  return result*DH;
}

#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
      
