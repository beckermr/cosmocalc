#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

double hubble_noscale(double a)
{
  double WA = 0.0;
  double W0 =1.0;
  return sqrt(cosmoData.OmegaM/a/a/a + (1.0-cosmoData.OmegaM)*exp(3.0*(WA*(a-1) - log(a)*(1.0 + W0 + WA))));
}
