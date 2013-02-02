#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>

#include "cosmocalc.h"
#include "haloprofs.h"

double fourierTransformNFW(double k, double m, double rvir, double c, double a)
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

double NFWprof(double r, double m, double rvir, double c, double a)
{
  double rs,rhos;
  
  rs = rvir/c;
  rhos = m/4.0/M_PI/rs/rs/rs/(log(1.0+c) - c/(1.0+c));
  
  return rhos/(r/rs)/(1+r/rs)/(1+r/rs);
}

double concNFW(double m, double a)
{
  static int init = 1;
  static int currCosmoNum;
  static double mstar;
  
  if(init || currCosmoNum != cosmoData.cosmoNum)
    {
      init = 0;
      currCosmoNum = cosmoData.cosmoNum;
      
      mstar = get_linear_tophatnorm_scale(1.0);
      mstar = 4.0*M_PI/3.0*mstar*mstar*mstar*cosmoData.OmegaM*RHO_CRIT;
    }
  
  return 10.0*pow(m/mstar,-0.2)*a;
}

