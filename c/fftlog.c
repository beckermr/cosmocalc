#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_result.h>
#include <fftw3.h>

#include "cosmocalc.h"

void compute_discrete_spherical_fft(double *data, int N, double r0, double L, double q, double *result, double *k0)
{
  double mu,rn,kn,tmp[2],um[2],k0guess;
  int i,m;
  static int initFlag = 1;
  static double *rin;
  static fftw_complex *out; 
  static fftw_plan pf,pb;
  
  if(initFlag)
    {
      initFlag = 0;
      
      rin = (double*)fftw_malloc(sizeof(double)*N); 
      out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1));
      pf = fftw_plan_dft_r2c_1d(N,rin,out,FFTW_MEASURE);
      pb = fftw_plan_dft_c2r_1d(N,out,rin,FFTW_MEASURE);
    }
  
  for(i=0;i<N;++i)
    {
      rn = r0*exp(L/((double) N)*((double) (i-N/2)));
      rin[i] = data[i]*4.0*M_PI*pow(rn,1.5)*sqrt(M_PI/2.0);
    }
  fftw_execute(pf);
  
  q = 0.0;
  mu = 0.5;
  k0guess = *k0;
  *k0 = get_k0_fftlog(N,mu,q,r0,L,k0guess);
  for(i=0;i<N/2+1;++i)
    {
      if(i >= N/2)
	m = i-N;
      else
	m = i;
      
      um_fftlog(m,mu,q,*k0,r0,L,&(um[0]),&(um[1]));
      if(m == N/2 || m == -N/2)
	um[1] = 0.0;
      
      tmp[0] = out[i][0]*um[0] - out[i][1]*um[1];
      tmp[1] = out[i][1]*um[0] + out[i][0]*um[1];
      
      out[i][0] = tmp[0];
      out[i][1] = tmp[1];
    }
  
  fftw_execute(pb);
  
  for(i=0;i<N;++i)
    {
      kn = (*k0)*exp(L/((double) N)*((double) (N/2-i)));
      result[N-i-1] = rin[i]/((double) N)/pow(kn,1.5)*exp(L/((double) N)*3.0/2.0);
    }
  
  /* not needed 
     fftw_destroy_plan(pf);
     fftw_destroy_plan(pb);
     
     fftw_free(rin);
     fftw_free(out);
  */
}

double get_k0_fftlog(int N, double mu, double q, double r0, double L, double k0guess)
{
  double zr,zi,arg;
  gsl_sf_result lnr_top,lnr_bot;
  gsl_sf_result arg_top,arg_bot;
  double twox[2];
  double grat[2];
  double tmp[2];
  double lnkr;
  int intval;
  
  //do top part of gamma function ratio
  zr = (mu + 1.0 + q)/2.0;
  zi = M_PI*((double) N)/L/2.0;
  gsl_sf_lngamma_complex_e(zr,zi,&lnr_top,&arg_top);
  
  //do bottom
  zr = (mu + 1.0 - q)/2.0;
  zi = -1.0*M_PI*((double) N)/L/2.0;
  gsl_sf_lngamma_complex_e(zr,zi,&lnr_bot,&arg_bot);
  
  //get real and imag parts of 2^(q+i*2*pi*m/L)
  zr = pow(2.0,q);
  twox[0] = zr*cos(M_PI*((double) N)/L*log(2.0));
  twox[1] = zr*sin(M_PI*((double) N)/L*log(2.0));
  
  //now do computation of real and imag parts of result
  zr = exp(lnr_top.val-lnr_bot.val);
  grat[0] = zr*cos(arg_top.val-arg_bot.val);
  grat[1] = zr*sin(arg_top.val-arg_bot.val);
  
  //do product
  tmp[0] = grat[0]*twox[0] - grat[1]*twox[1];
  tmp[1] = grat[0]*twox[1] + grat[1]*twox[0];
  
  //find k0
  lnkr = log(k0guess*r0);
  arg = atan2(tmp[1],tmp[0])/M_PI*L/((double) N);
  intval = (int) ((lnkr-arg)/(L/((double) N)));
    
  return exp(arg + ((double) intval)*L/((double) N))/r0;
}

void um_fftlog(int m, double mu, double q, double k0, double r0, double L, double *realpart, double *imagpart)
{
  double zr,zi;
  gsl_sf_result lnr_top,lnr_bot;
  gsl_sf_result arg_top,arg_bot;
  double kr[2];
  double twox[2];
  double grat[2];
  double tmp[2];
  
  //do top part of gamma function ratio
  zr = (mu + 1.0 + q)/2.0;
  zi = M_PI*((double) m)/L;
  gsl_sf_lngamma_complex_e(zr,zi,&lnr_top,&arg_top);
    
  //do bottom
  zr = (mu + 1.0 - q)/2.0;
  zi = -1.0*M_PI*((double) m)/L;
  gsl_sf_lngamma_complex_e(zr,zi,&lnr_bot,&arg_bot);
    
  //get real and imag parts of (k0*r0)^(i*2*pi*m/L)
  zr = -1.0*log(k0*r0)*2.0*M_PI*((double) m)/L;
  kr[0] = cos(zr);
  kr[1] = sin(zr);
    
  //get real and imag parts of 2^(q+i*2*pi*m/L)
  zr = pow(2.0,q);
  twox[0] = zr*cos(2.0*M_PI*m/L*log(2.0));
  twox[1] = zr*sin(2.0*M_PI*m/L*log(2.0));
    
  //now do computation of real and imag parts of result
  zr = exp(lnr_top.val-lnr_bot.val);
  grat[0] = zr*cos(arg_top.val-arg_bot.val);
  grat[1] = zr*sin(arg_top.val-arg_bot.val);
    
  tmp[0] = kr[0]*twox[0] - kr[1]*twox[1];
  tmp[1] = kr[0]*twox[1] + kr[1]*twox[0];
    
  *realpart = tmp[0]*grat[0] - tmp[1]*grat[1];
  *imagpart = tmp[0]*grat[1] + tmp[1]*grat[0];
}
