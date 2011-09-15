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
  
  //init cosmology
  cosmoData.cosmoNum = 1;
  cosmoData.OmegaM = 0.3;
  cosmoData.h = 0.7;
  cosmoData.Sigma8 = 0.8;
  cosmoData.SpectralIndex = 1.0;
  cosmoData.OmegaB = 0.04;
  cosmoData.delta = -1;
  
#ifdef TEST_CODE
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
#else
  //compute what you want here
  //see test_code.c for examples
  
  int pval = atoi(argv[1]);
  int i,Nk;
  double a,k,pkL,sn;
  FILE *fp;
  
  double m,mmin,dlnm,mmax,rmin,rmax,dr,intg;
  int Nm,Nr,j;
  double z,s,dndm,tf;
  
  if(pval == -1)
    {
      z = 0.325;
      m = 1e15;
      a = 1.0/(1.0+z);
      fprintf(stderr,"test for michael:\n");
      fprintf(stderr,"z = %g\n",z);
      fprintf(stderr,"Delta = %g\n",cosmoData.delta = 200.0*(cosmoData.OmegaM/a/a/a+1.0-cosmoData.OmegaM)/(pow(1.0/a,3.0)*cosmoData.OmegaM));
      cosmoData.delta = -cosmoData.delta;
      fprintf(stderr,"mf = %g\n",mass_function(m,a));
      
      //WMAP 1
      cosmoData.cosmoNum = 2;
      cosmoData.OmegaM = 0.3;
      cosmoData.h = 0.7;
      cosmoData.Sigma8 = 0.9;
      cosmoData.SpectralIndex = 1.0;
      cosmoData.OmegaB = 0.04;
      cosmoData.delta = -1;
      
      fprintf(stderr,"log_{10}(1/\\sigma(m=1e10,a=1)) = %lg\n",log10(1.0/sigmaMtophat(1e10,1.0)));
      fprintf(stderr,"log_{10}(1/\\sigma(m=1e16,a=1)) = %lg\n",log10(1.0/sigmaMtophat(1e16,1.0)));
    }
  else if(pval == -2)
    {
      fp = fopen("../data/testing_0.3.dndM","r");
      z = 0.3;
      a = 1.0/(1.0+z);
      cosmoData.delta = -200.0*(cosmoData.OmegaM/a/a/a+1.0-cosmoData.OmegaM)/(pow(1.0/a,3.0)*cosmoData.OmegaM);
      while(1)
	{
	  fscanf(fp,"%lg %lg %lg %lg\n",&m,&dndm,&s,&tf);
	  if(feof(fp))
	    break;
	  k = 1.0/pow(m/(4.0/3.0*M_PI*RHO_CRIT*cosmoData.OmegaM),1.0/3.0)/2.0/M_PI;
	  fprintf(stderr,"m = %lg, dndm = %lg (?= %lg), sigma(m) = %lg (?= %lg), Tf(%lg) = %lg (?= %lg)\n",
		  m,mass_function(m,a),dndm,sigmaMtophat(m,a),s,
		  k,transfer_function(k),tf);
	  fprintf(stdout,"%lg %lg %lg\n",m,mass_function(m,a),sigmaMtophat(m,a));
	}
      fclose(fp);
    }
  else if(pval == 1) 
    fprintf(stdout,"%.20e\n",transfunct_eh98(atof(argv[2])));
  else if(pval == 2)
    {
      fp = fopen(argv[2],"r"); 
      fscanf(fp,"%lg %*lg\n",&a);
      fprintf(stderr,"gf = %g\n",1.0/growth_function(1.0/(1.0+a)));
      
      while(1)
	{
	  fscanf(fp,"%lg %lg    %*lg %*lg\n",&k,&pkL);
	  
	  if(feof(fp))
	    break;
	  
	  fprintf(stdout,"%.20e %.20e %.20e\n",k,pkL,linear_powspec(k,1.0)*k*k*k*4.0*M_PI);
	}
      fclose(fp);
    }
  else if(pval == 3)
    fprintf(stdout,"%.20e\n",growth_function_exact(atof(argv[2])));
  else if(pval == 4)
    {
      sn = atof(argv[3]);
      fp = fopen(argv[2],"r"); 
      fscanf(fp,"%lg\n",&a);
      fscanf(fp,"%d\n",&Nk);
      
      for(i=0;i<Nk;++i)
	{
	  fscanf(fp,"%lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg %*lg\n",&k);
	  fprintf(stdout,"%.20e %.20e %.20e\n",k,(linear_powspec(k,a)+sn)*k*k*k/2.0/M_PI/M_PI,(nonlinear_powspec(k,a)+sn)*k*k*k/2.0/M_PI/M_PI);
	}
      
      fclose(fp);
    }
  else if(pval == 5)
    {
      mmin = atof(argv[2])/100.0;
      mmax = atof(argv[3])*100.0;
      rmin = atof(argv[4]);
      rmax = atof(argv[5]);
      Nm = 10000;
      dlnm = log(mmax/mmin)/Nm;
      fprintf(stderr,"a = %g|%g\n",acomvdist(rmin),acomvdist(rmax));
      fprintf(stderr,"z = %g|%g\n",1.0/acomvdist(rmin)-1.0,1.0/acomvdist(rmax)-1.0);
      fprintf(stderr,"m = %g|%g\n",mmin,mmax);
      
      a = (acomvdist(rmax)+acomvdist(rmin))/2.0;
      
      Nr = 100.0;
      dr = (rmax-rmin)/Nr;
      
      for(i=0;i<Nm;++i)
	{
	  m = exp(i*dlnm+0.5*dlnm)*mmin;
	  intg = 0;
	  for(j=0;j<Nr;++j)
	    {
	      a = acomvdist(dr*j+rmin);
	      if(strcmp(argv[6],"m") == 0)
		cosmoData.delta = 200.0;
	      else
		{
		  cosmoData.delta = 200.0*(cosmoData.OmegaM/a/a/a+1.0-cosmoData.OmegaM)/(pow(1.0/a,3.0)*cosmoData.OmegaM);
		}
	      intg += mass_function(m,a)*m*4.0*M_PI*pow(j*dr+rmin,2.0)*dr;
	      //fprintf(stderr,"m = %g, a = %g, delta = %g, mf = %g\n",m,a,cosmoData.delta,mass_function(m,a)*m);
	    }
	  intg = intg/(4.0/3.0*M_PI*(pow(rmax,3.0)-pow(rmin,3.0)));
	  
	  fprintf(stdout,"%.20e %.20e\n",m,intg);
	}
    }
  
#endif /* TEST_CODE */  
  
  return 0;
}
