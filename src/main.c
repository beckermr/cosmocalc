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

void test_cosmo(char *path);

int main(int argc, char **argv)
{
  cosmoData.useSmoothTransFunc = 0;
  cosmoData.delta = 200.0;
  char path[4096];
  
  //init cosmology 2
  cosmoData.cosmoNum = 2;
  cosmoData.OmegaM = 0.28;
  cosmoData.OmegaL = 0.72;
  cosmoData.OmegaB = 0.05;
  cosmoData.OmegaNu = 0.0;
  cosmoData.OmegaK = 0.0;
  cosmoData.h = 0.7;
  cosmoData.w0 = -1.0;
  cosmoData.wa = 0.0;
  cosmoData.SpectralIndex = 0.96;
  cosmoData.Sigma8 = 0.8;

  sprintf(path,"./cosmo2");
  test_cosmo(path);
  
  //init cosmology 1
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.3;
  cosmoData.OmegaL = 0.68;
  cosmoData.OmegaB = 0.05;
  cosmoData.OmegaNu = 0.0;
  cosmoData.OmegaK = 1.0 - cosmoData.OmegaM - cosmoData.OmegaL - cosmoData.OmegaNu;
  cosmoData.h = 0.7;
  cosmoData.w0 = -0.8;
  cosmoData.wa = 0.2;
  cosmoData.SpectralIndex = 0.96;
  cosmoData.Sigma8 = 0.8;
  
  sprintf(path,"./cosmo1");
  test_cosmo(path);
  
  return 0;
}

void test_cosmo(char *path)
{
  char fname[4096];
  char fnameo[4096];
  char str[4096];
  FILE *fp;
  FILE *fpout;
  double z;
  
  sprintf(fname,"%s/codist.dat",path);
  sprintf(fnameo,"%s/codist.dat.matt",path);
  fp = fopen(fname,"r");
  if( fp != NULL)
    {
      fpout = fopen(fnameo,"w");
      while(fgets(str,4096,fp) != NULL)
	{
	  if(str[0] == '#')
	    continue;
	  sscanf(str,"%le %*e\n",&z);
	  
	  fprintf(fpout,"%e %e\n",z,comvdist(1.0/(1.0+z)));
	}
      fclose(fp);
      fclose(fpout);
    }
  
  sprintf(fname,"%s/Evol_z.dat",path);
  sprintf(fnameo,"%s/Evol_z.dat.matt",path);
  fp = fopen(fname,"r");
  if(fp != NULL)
    {
      fpout = fopen(fnameo,"w");
      while(fgets(str,4096,fp) != NULL)
	{
	  if(str[0] == '#')
	    continue;
	  sscanf(str,"%le %*e\n",&z);
	  
	  //fprintf(stderr,"weff(%f) = %f, h(a) = %f\n",1.0/(1.0+z),weff(1.0/(1.0+z)),hubble_noscale(1.0/(1.0+z)));
	  fprintf(fpout,"%e %e\n",z,hubble_noscale(1.0/(1.0+z)));
	}
      fclose(fp);
      fclose(fpout);
    }
}
