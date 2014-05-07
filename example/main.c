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

double solve_scale_nfw(gsl_root_fsolver *s, gsl_function F, double *p)
{
  int status;
  int iter,max_iter = 100;
  double x_lo = 0.0,x_hi = 5.0;
  double r;
  
  iter = 0;
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
  
  return r;
}

int main(int argc, char **argv)
{
  //init
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.27;
  cosmoData.OmegaL = 0.73;
  cosmoData.OmegaB = 0.045;
  cosmoData.OmegaNu = 0.0;
  cosmoData.OmegaK = 0.0;
  cosmoData.h = 0.7;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.SpectralIndex = 0.96;
  cosmoData.Sigma8 = 0.8;
  cosmoData.useSmoothTransFunc = 0;
    
  //prints mass function in bins to stdout
  double a = atof(argv[1]);
  double mp = atof(argv[2]);
  double np;
  double mmin = 1e1;
  double mmax = 1e18;
  long Nm = 2000;
  double m,dlnm = log(mmax/mmin)/Nm;
  long i;
  double rvir,rs,c;
  
  cosmoData.delta = 200.0*RHO_CRIT*hubble_noscale(a)*hubble_noscale(a)/(cosmoData.OmegaM*RHO_CRIT/a/a/a);
  //fprintf(stderr,"delta = %f\n",cosmoData.delta);
  
  gsl_root_fsolver *s;
  gsl_function F;
  double p[4],r1,rN,rN2;
  
  F.function = &eps_fun;
  F.params = p;
  s = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    
  fprintf(stdout,"# m rs rvir r1 rN rN2\n");
  for(i=0;i<Nm;++i)
    {
      m = exp(dlnm*i)*mmin;
      rvir = pow(m/(cosmoData.delta*RHO_CRIT*cosmoData.OmegaM*4.0/3.0*M_PI),1.0/3.0);
      c = duffy2008_concNFW(m,a);
      //if(c < 4.0)
      //c = 4.0;
      rs = rvir/c;
      
      np = m/mp;
      
      if(np >= 1.0)
	{
	  p[0] = m;
	  p[1] = rvir;
	  p[2] = c;
	  p[3] = mp;
	  r1 = solve_scale_nfw(s,F,p);
	}
      else
	r1 = -1.0;
      
      if(np > 100)
	{
	  p[0] = m;
	  p[1] = rvir;
	  p[2] = c;
	  p[3] = mp*100.0;
	  rN = solve_scale_nfw(s,F,p);
	}
      else
	rN = -1.0;
      
      if(np > 10)
	{
	  p[0] = m;
	  p[1] = rvir;
	  p[2] = c;
	  p[3] = mp*10.0;
	  rN2 = solve_scale_nfw(s,F,p);
	}
      else
	rN2 = -1.0;
      
      fprintf(stdout,"%e %e %e %e %e %e\n",m,rs,rvir,r1,rN,rN2);
    }
  
  gsl_root_fsolver_free(s);
  
  return 0;
}

