#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "cosmocalc.h"
#include "haloprofs.h"
#include "weaklens.h"

double eps_fun(double x, void *p)
{
  double *v = (double*)p;
  return NFWprof_menc(x,v[0],v[1],v[2]) - v[3];
}

int main(int argc, char **argv)
{
  //init cosmology
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.27;
  cosmoData.OmegaL = 0.73;
  cosmoData.OmegaB = 0.045;
  cosmoData.OmegaK = 0.0;
  cosmoData.OmegaNu = 0.0;
  cosmoData.h = 0.7;
  cosmoData.SpectralIndex = 0.96;
  cosmoData.Sigma8 = 0.8;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
    
  cosmoData.useSmoothTransFunc = 0;
    
  //prints mass function in bins to stdout
  double a = atof(argv[1]);
  double np = atof(argv[2]);
  double mmin = 1e7;
  double mmax = 1e16;
  long Nm = 2000;
  double m,dlnm = log(mmax/mmin)/Nm;
  long i;
  double rvir,rs,c;
  
  cosmoData.delta = 200.0*RHO_CRIT*hubble_noscale(a)*hubble_noscale(a)/(cosmoData.OmegaM*RHO_CRIT/a/a/a);
  //fprintf(stderr,"delta = %f\n",cosmoData.delta);
  
  int status;
  int iter,max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double x_lo,x_hi;
  gsl_function F;
  double p[4];
  double r;
     
  F.function = &eps_fun;
  F.params = p;
  
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
    
  fprintf(stdout,"# m rs rvir rp\n");
  for(i=0;i<Nm;++i)
    {
      m = exp(dlnm*i)*mmin;
      rvir = pow(m/(cosmoData.delta*RHO_CRIT*cosmoData.OmegaM*4.0/3.0*M_PI),1.0/3.0);
      c = duffy2008_concNFW(m,a);
      if(c < 4.0)
	c = 4.0;
      rs = rvir/c;
      
      iter = 0;
      r = 1.0;
      x_lo = 0.0;
      x_hi = 5.0;
      p[0] = m;
      p[1] = rvir;
      p[2] = c;
      p[3] = m/np;
      gsl_root_fsolver_set(s,&F,x_lo,x_hi);
      do
	{
	  iter++;
	  status = gsl_root_fsolver_iterate(s);
	  r = gsl_root_fsolver_root(s);
	  x_lo = gsl_root_fsolver_x_lower(s);
	  x_hi = gsl_root_fsolver_x_upper(s);
	  status = gsl_root_test_interval(x_lo, x_hi,0,0.001);
	}
      while (status == GSL_CONTINUE && iter < max_iter);
            
      fprintf(stdout,"%e %e %e %e\n",m,rs,rvir,r);
    }
  
  gsl_root_fsolver_free(s);
  
  return 0;
}
