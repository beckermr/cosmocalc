#include "cosmocalc_assert.h"
#include "cosmocalc.h"

#ifndef _HALOMODEL_
#define _HALOMODEL_

//standard physical constants
#define DELTAC 1.686

class MassFunction {
  virtual double operator()(double m, double a) {
    cosmocalc_assert(false,"In operator() method, MassFunction base class cannot be used alone!");
  };
};

class HaloBias {
  virtual double operator()(double m, double a) {
    cosmocalc_assert(false,"In operator() method, HaloBias base class cannot be used alone!");
  };  
};

class HaloProfile {
  virtual double operator()(double r, double m, double a) {
    cosmocalc_assert(false,"In operator() method, HaloProfile base class cannot be used alone!");
  };
};

class HaloConc {
  virtual double operator()(double m, double a) {
    cosmocalc_assert(false,"In operator() method, HaloConc base class cannot be used alone!");
  };
};

#endif /* _HALOMODEL_ */
