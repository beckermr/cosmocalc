#include "cosmocalc_assert.h"
#include "cosmocalc.h"

#ifndef _HALOMODEL_
#define _HALOMODEL_

//standard physical constants
#define DELTAC 1.686

class MassFunction {
  void double operator()(double m, double a) {
    cosmocalc_assert(false,"In operator() method, MassFunction base class cannot be used alone!");
  };
};

class HaloBias {
  void double operator()(double m, double a) {
    cosmocalc_assert(false,"In operator() method, HaloBias base class cannot be used alone!");
  };  
};

class PeakHeight {
  virtual double sigmaRtophat(double topHatRad, double a) {
    cosmocalc_assert(false,"In sigmaRtophat method, PeakHeight base class cannot be used alone!");
  };
  virtual double Rsigmatophat(double sigmaR, double a) {
    cosmocalc_assert(false,"In Rsigmatophat method, PeakHeight base class cannot be used alone!");
  };
}

#endif /* _HALOMODEL_ */
