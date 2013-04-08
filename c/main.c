#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"
#include "haloprofs.h"
#include "weaklens.h"

int main(int argc, char **argv)
{
  //init cosmology
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.272187;
  cosmoData.OmegaL = 0.727813;
  cosmoData.OmegaB = 0.0455998;
  cosmoData.OmegaK = 0.0;
  cosmoData.OmegaNu = 0.0;
  cosmoData.h = 0.704;
  cosmoData.As = 2.449e-9;
  cosmoData.As_pivot = 0.002;
  cosmoData.SpectralIndex = 0.963;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.delta = 200.0;
  
  cosmoData.useSmoothTransFunc = 0;
  
  cosmoData.Sigma8 = convert_cmbnorm2sigma8();
  fprintf(stderr,"for As(k=%g Mpc^-1) = %g, sigma8 = %g\n",cosmoData.As_pivot,cosmoData.As,cosmoData.Sigma8);
  
  //prints mass function in bins to stdout
  double a = 1.0;
  
  double mmin = 1e9;
  double mmax = 1e16;
  long Nm = 2000;
  double m,dlnm = log(mmax/mmin)/Nm;
  long i;
  double sigma;
  double dm;
  double dlnsiginvdm;

  fprintf(stdout,"# m t08 t10 s ds\n");
  for(i=0;i<Nm;++i)
    {
      m = exp(dlnm*i)*mmin;
      sigma = sigmaMtophat(m,a);
      dm = 1e-6*m;
      dlnsiginvdm = log(sigmaMtophat(m-dm/2.0,a)/sigmaMtophat(m+dm/2.0,a))/dm;
      fprintf(stdout,"%e %e %e %e %e\n",m,tinker2008_mass_function(m,a,cosmoData.delta),
              tinker2010_mass_function(m,a,cosmoData.delta),sigma,dlnsiginvdm);
    }
  
  return 0;
}
