#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"
#include "weaklens.h"

cosmocalcData cosmoData;
weaklensData wlData;

int main(int argc, char **argv)
{
  
#ifndef TEST_CODE
  //init cosmology
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.25;
  cosmoData.h = 0.7;
  cosmoData.Sigma8 = 0.8;
  cosmoData.SpectralIndex = 1.0;
  cosmoData.OmegaB = 0.04;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.useSmoothTransFunc = 0;
  
  double a = atof(argv[1]);
  
  double kmin = 1e-6;
  double kmax = 1e4;
  long Nk = 2000;
  double k,dlnk = log(kmax/kmin)/Nk;
  long i;

  FILE *fp;
  fp = fopen("test.dat","w");
  for(i=0;i<Nk;++i)
    {
      k = exp(dlnk*i)*kmin;
      //fprintf(fp,"%e %e %e %e %e\n",k,tf(k),tfs(k),lp(k,a),pknl(k,a));                                                                                                                                                                                                                    
      fprintf(fp,"%e %e %e %e %e\n",k,transfer_function(k),transfer_function(k),linear_powspec(k,a),nonlinear_powspec(k,a));
      //fprintf(fp,"%e %e %e %e %e\n",k,transfer_function(k),transfer_function(k),linear_powspec(k,a),nonlinear_powspec_for_lens(k,a));
    }
  fclose(fp);
  
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
