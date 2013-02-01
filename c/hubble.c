#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

double weff(double a) 
{
  if(a != 1.0)
    return cosmoData.w0 + cosmoData.wa - cosmoData.wa*(a - 1.0)/log(a);
  else
    return cosmoData.w0;
}

double hubble_noscale(double a)
{
  return sqrt(cosmoData.OmegaM/a/a/a + cosmoData.OmegaK/a/a + cosmoData.OmegaL*exp(3.0*(1.0 + weff(a))));
}
