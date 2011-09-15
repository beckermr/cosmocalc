#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

void test_nonlinear_corrfunc(void)
{
  double lnrmin = log(1e-5);
  double lnrmax = log(3e2);
  int N = 10000;
  double r,dlnr = (lnrmax - lnrmin)/N;
  int i;
  FILE *fp;
  double a;
  
  a = 1.0;
  fp = fopen("xinltest.txt","w");
  for(i=0;i<N;++i)
    {
      r = exp(dlnr*i + lnrmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",r,nonlinear_corrfunc(r,a),linear_corrfunc(r,a));
    }
  fclose(fp);
}

void test_nonlinear_powspec(void)
{
  double lnkmin = log(1e-4);
  double lnkmax = log(1e2);
  double dlnk,k;
  int N = 10000,i;
  FILE *fp;
  double a = 1.0;
  
  dlnk = (lnkmax - lnkmin)/N;
  
  a = 1.0;
  fp = fopen("pknltest.txt","w");
  for(i=0;i<N;++i)
    {
      k = exp(dlnk*i + lnkmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
    }
  fclose(fp);
  
  a = 1.0/(1.0 + 0.5);
  fp = fopen("pknltest1.txt","w");
  for(i=0;i<N;++i)
    {
      k = exp(dlnk*i + lnkmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
    }
  fclose(fp);
  
  a = 1.0/(1.0 + 1.0);
  fp = fopen("pknltest2.txt","w");
  for(i=0;i<N;++i)
    {
      k = exp(dlnk*i + lnkmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
    }
  fclose(fp);
  
  a = 1.0/(1.0 + 2.00);
  fp = fopen("pknltest3.txt","w");
  for(i=0;i<N;++i)
    {
      k = exp(dlnk*i + lnkmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
    }
  fclose(fp);

  a = 1.0/(1.0 + 3.00);
  fp = fopen("pknltest4.txt","w");
  for(i=0;i<N;++i)
    {
      k = exp(dlnk*i + lnkmin);
      fprintf(fp,"%.20e\t%.20e\t%.20e\n",k,nonlinear_powspec(k,a),linear_powspec(k,a));
    }
  fclose(fp);
}

void test_linxi_fftlog(void)
{
  double rmin;
  long i,N = 2048;
  double dlnr;
  FILE *fp;
  double r;
  double a = 1.0;
  
  double *pk,*xi;
  double k0,L,q,r0,k;
  double kmin = 1e-7;
  double kmax = 1e7;
  double dlnk = log(kmax/kmin)/N;
  
  L = log(kmax/kmin);
  k0 = exp(log(kmin) + L/2.0);
  q = 0.0;
  
  pk = (double*)malloc(sizeof(double)*N);
  assert(pk != NULL);
  xi = (double*)malloc(sizeof(double)*N);
  assert(xi != NULL);
  
  for(i=0;i<N;++i)
    {
      k = exp(i*dlnk + log(kmin));
      pk[i] = linear_powspec_exact(k,a);
    }
  
  r0 = 1.0;
  compute_discrete_spherical_fft(pk,N,k0,L,q,xi,&r0);
  
  if(0)
    {
      for(i=0;i<N;++i)
	xi[i] = xi[i]/pow(2.0*M_PI,3.0);
      compute_discrete_spherical_fft(xi,N,r0,L,q,pk,&k0);
      
      kmin = exp(log(k0)-L/2.0);
      //fprintf(stderr,"k0 = %e, kmin = %e, L = %f, N = %ld\n",k0,kmin,L,N);
      fp = fopen("linxitest.txt","w");
      for(i=0;i<N;++i)
	{
	  k = exp(i*dlnk + log(kmin));
	  //fprintf(stderr,"k = %e, pk in,out = %e|%e\n",k,linear_powspec(k,a),pk[i]);
	  fprintf(fp,"%.20e \t %.20e \t %.20e \n",k,pk[i],linear_powspec_exact(k,a));
	}
      fclose(fp);
      
      return;
    }
  
  rmin = exp(log(r0)-L/2.0);
  dlnr = dlnk;
  fp = fopen("linxitest.txt","w");
  for(i=0;i<N;++i)
    {
      r = exp(i*dlnr + log(rmin));
      xi[i] = xi[i]/pow(2.0*M_PI,3.0);
      //fprintf(stderr,"%.20e \t %.20e \t %.20e \n",r,xi[i],linear_corrfunc_exact(r,a));
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",r,xi[i],linear_corrfunc_exact(r,a));
    }
  fclose(fp);
}

void test_linxi(void)
{
  double rmin=1e-9;
  double rmax=1e5;
  long i,Ntest = 1000;
  double dlnr = log(rmax/rmin)/Ntest;
  FILE *fp;
  double r;
  double a = 1.0;
  
  fp = fopen("linxitest.txt","w");
  for(i=0;i<Ntest;++i)
    {
      r = exp(i*dlnr + log(rmin));
      //fprintf(stderr,"%.20e \t %.20e \t %.20e \n",r,linear_corrfunc_exact(r,a),linear_corrfunc_exact(r,a));
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",r,linear_corrfunc(r,a),linear_corrfunc_exact(r,a));
    }
  fclose(fp);
}

void test_biasfunc(void)
{
  FILE *fp;
  double m,lgm;
  double mmin = 6.0;
  double mmax = 16.0;
  double dm = 0.01;
  double a = 1.0;
  
  fp = fopen("btest.txt","w");
  lgm = mmax;
  while(lgm >= mmin)
    {
      m = pow(10.0,lgm);
      fprintf(fp,"%.20e \t %.20e \n",1.686/sigmaMtophat(m,a),bias_function(m,a));
      lgm -= dm;
    }
  fclose(fp);

}

void test_massfunc(void)
{
  FILE *fp;
  double m,lgm;
  double mmin = 10.0;
  double mmax = 16.0;
  double dm = 0.01;
  double a = 1.0;
  
  fp = fopen("mftest.txt","w");
  lgm = mmax;
  while(lgm >= mmin)
    {
      m = pow(10.0,lgm);
      fprintf(fp,"%.20e \t %.20e \n",m,mass_function(m,a)*m*m/RHO_CRIT/cosmoData.OmegaM);
      lgm -= dm;
    }
  fclose(fp);
}

void test_peakheight(void)
{
  FILE *fp;
  double r;
  double rmin = 1e-40;
  double rmax = 1e3;
  double dr = 0.99;
  double a = 1.0;
  double numin = 0.15;
  double numax = 6.0;
  double dnu = 0.99;
  double nu;
  
  double kmin=1e-10;
  double kmax=1e10;
  long i,Ntest = 1000;
  double dlnk = log(kmax/kmin)/Ntest;
  double k,p = 1e-6;
  
  fp = fopen("phtest_integ.txt","w");
  for(i=0;i<Ntest;++i)
    {
      k = (i*dlnk + log(kmin));
      fprintf(fp,"%.20e \t %.20e \n",exp(k),tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(k,&p));
    }
  fclose(fp);
  
  fp = fopen("phtest.txt","w");
  r = rmax;
  while(r >= rmin)
    {
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",r,sigmaRtophat(r,a),sigmaRtophat_exact(r,a));
      r *= dr;
    }
  fclose(fp);

  fp = fopen("phtestrev.txt","w");
  nu = numax;
  while(nu >= numin)
    {
      fprintf(fp,"%.20e \t %.20e \n",inverse_nuRtophat(nu,a),DELTAC/nu);
      nu *= dnu;
    }
  fclose(fp);
}

void test_gf(void)
{
  FILE *fp;
  double a;
  double amin = 0.1;
  double amax = 1.0;
  double da = 0.0001;
  
  fp = fopen("gftest.txt","w");
  a = amax;
  while(a >= amin)
    {
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",a,growth_function(a),growth_function_exact(a));
      a -= da;
    }
  fclose(fp);
}

void test_linpk(void)
{
  double kmin=1e-8;
  double kmax=1e40;
  long i,Ntest = 10000;
  double dlnk = log(kmax/kmin)/Ntest;
  FILE *fp;
  double k;
  double a = 1.0;
  
  fp = fopen("pktest.txt","w");
  for(i=0;i<Ntest;++i)
    {
      k = exp(i*dlnk + log(kmin));
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",k,linear_powspec(k,a),linear_powspec_exact(k,a));
    }
  fclose(fp);

}

void test_transfunct(void)
{
  double kmin=1e-8;
  double kmax=1e30;
  long i,Ntest = 10000;
  double dlnk = log(kmax/kmin)/Ntest;
  FILE *fp;
  double k;
  
  fp = fopen("tftest.txt","w");
  for(i=0;i<Ntest;++i)
    {
      k = exp(i*dlnk + log(kmin));
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",k,transfer_function(k),transfunct_eh98(k));
    }
  fclose(fp);
  
  cosmoData.OmegaM = 0.2;
  cosmoData.h = 0.5;
  cosmoData.OmegaB = cosmoData.OmegaM*0.5;
  cosmoData.cosmoNum = 2;
  
  fp = fopen("tftest_highb.txt","w");
  for(i=0;i<Ntest;++i)
    {
      k = exp(i*dlnk + log(kmin));
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",k,transfer_function(k),transfunct_eh98(k));
    }
  fclose(fp);
}

void test_distances(void)
{
  FILE *fp;
  
  double amin = 0.1;
  double amax = 1.0;
  double da = 0.0001;
  double dmin = 0.0;
  double dmax = 3000.0;
  double dd = 1.0;
  double a,d;
    
  fp = fopen("disttest_a2d.txt","w");
  a = amax;
  while(a >= amin)
    {
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",a,comvdist(a),comvdist_exact(a));
      a -= da;
    }
  fclose(fp);
  
  fp = fopen("disttest_d2a.txt","w");
  d = dmin;
  while(d <= dmax)
    {
      a = acomvdist(d);
      fprintf(fp,"%.20e \t %.20e \t %.20e \n",d,a,comvdist_exact(a));
      d += dd;
    }
  fclose(fp);
}
