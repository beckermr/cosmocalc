#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

double hubble_noscale(double a)
{
  return sqrt(cosmoData.OmegaM/a/a/a + (1.0-cosmoData.OmegaM)*exp(3.0*(cosmoData.wa*(a-1) - log(a)*(1.0 + cosmoData.w0 + cosmoData.wa))));
}
