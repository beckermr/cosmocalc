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
  cosmoData.OmegaM = 0.29;
  cosmoData.OmegaL = 0.71;
  cosmoData.OmegaB = 0.047;
  cosmoData.OmegaK = 0.0;
  cosmoData.OmegaNu = 0.0;
  cosmoData.h = 0.7;
  cosmoData.Sigma8 = 0.8;
  cosmoData.SpectralIndex = 0.96;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.delta = 200.0;
  
  cosmoData.useSmoothTransFunc = 0;
  
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
