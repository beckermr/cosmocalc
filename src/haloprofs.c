#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>

#include "cosmocalc.h"
#include "haloprofs.h"

double fourierTransformNFW(double k, double m, double rvir, double c)
{
  double rs,rhos;
  double ukm;
  
  rs = rvir/c;
  rhos = m/4.0/M_PI/rs/rs/rs/(log(1.0+c) - c/(1.0+c));
  
  if(k > 0)
    ukm = 4.0*M_PI*rhos*rs*rs*rs/m
      *(sin(k*rs)*(gsl_sf_Si((1.0+c)*k*rs) - gsl_sf_Si(k*rs)) - sin(c*k*rs)/(1.0+c)/k/rs + cos(k*rs)*(gsl_sf_Ci((1.0+c)*k*rs) - gsl_sf_Ci(k*rs)));
  else
    ukm = 1.0;
  
  return ukm;
}

double NFWprof(double r, double m, double rvir, double c)
{
  double rs,rhos;
  
  rs = rvir/c;
  rhos = m/4.0/M_PI/rs/rs/rs/(log(1.0+c) - c/(1.0+c));
  
  return rhos/(r/rs)/(1+r/rs)/(1+r/rs);
}

double NFWprof_menc(double r, double m, double rvir, double c)
{
  double mc = log(1+c) - c/(1.0+c);
  double mr = c*r/rvir;
  mr = log(1.0 + mr) - mr/(1.0+mr);
  
  return mr/mc*m;
}

inline double _duffy2008_concNFW(double m, double a)
{
  return 6.71/pow(1.0/a,0.44)*pow(m/2e12,-0.091);
}

double concNFW(double m, double a)
{
  return _duffy2008_concNFW(m,a);
}

double duffy2008_concNFW(double m, double a)
{
  return _duffy2008_concNFW(m,a);
}

/*
double bh2012_concNFW(double m, double a)
{
  double nu = DELTAC/sigmaMtophat(m,a);
  

}

double prada2011_concNFW(double m, double a)
{


}
*/


