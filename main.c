#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

cosmocalcData cosmoData;

int main(int argc, char **argv)
{
  
#ifndef TEST_CODE
  //init cosmology
  cosmoData.cosmoNum = 1;
  cosmoData.delta = 200.0; //delta > 0.0 Tinker+10, delta < 0.0 Tinker+08, always wrt to mean density
  cosmoData.OmegaM = 0.25;
  cosmoData.h = 0.7;
  cosmoData.Sigma8 = 0.8;
  cosmoData.SpectralIndex = 1.0;
  cosmoData.OmegaB = 0.045;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  
  if(0)
    {
      double lnkmin = log(1e-4);
      double lnkmax = log(1e2);
      double dlnk,k,a;
      int N = 10000,i;
      
      dlnk = (lnkmax - lnkmin)/N;
      
      a = atof(argv[2]);
      for(i=0;i<N;++i)
	{
	  k = exp(dlnk*i + lnkmin);
	  fprintf(stdout,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
	}
    }
  else
    {
      double amin = 1.0/31.0;
      double amax = 1.0;
      int Na = 100;
      double da = (amax-amin)/(Na-1.0);
      double a;
      int i;
            
      for(i=0;i<Na;++i)
	{
	  a = amin + da*i;
	  fprintf(stdout,"%.20e\t%.20e\n",a,growth_function(a));
	}
    }
    
#else
  test_nonlinear_corrfunc();
  test_nonlinear_powspec();
  test_linxi_fftlog();
  test_linxi();
  test_biasfunc();
  test_massfunc();
  test_peakheight();
  test_gf();
  test_distances();
  test_transfunct();
  test_linpk();
#endif
    
  return 0;
}
