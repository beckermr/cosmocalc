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
  cosmoData.OmegaM = 0.25;
  cosmoData.OmegaL = 0.75;
  cosmoData.OmegaB = 0.045;
  cosmoData.OmegaK = 0.0;
  cosmoData.OmegaNu = 0.0;
  cosmoData.h = 0.7;
  cosmoData.Sigma8 = 0.8;
  cosmoData.SpectralIndex = 1.0;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.delta = 200.0;
  
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
  
  double m = atof(argv[2]);
  fprintf(stderr,"sigma8 = %f\n",sqrt(tophatnorm_linear_powspec(8.0)));
  fprintf(stderr,"c(%e,%f) = %f\n",m,a,concNFW(m,a));
  
  return 0;
}
