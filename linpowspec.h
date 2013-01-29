#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

#ifndef _LINPOWSPEC_
#define _LINPOWSPEC_

class LinearPowerSpectrum : public PowerSpectrum {
 protected:
  double _ns;
  double _s8;
  class GrowthFunction *_gf;
  class TransferFunction *_tf;
  
  int LINEAR_POWSPEC_TABLE_LENGTH; //number of spline points in k for linear power spec
  double PL_K_MIN;  //min k value for linear powspec in h/Mpc
  double PL_K_MAX;  //max k value for linear powspec in h/Mpc
  int LINEAR_POWSPEC_FIT_LENGTH; //number of points right before PL_K_MAX to use to fir pl(k) for power-law extrapolation
  double tophatradnorm_linear_powspec_exact_nonorm(double topHatRad); //returns value of unormalized linear powspec filtered by top hat of radius topHatRad
  double _linear_powspec_norm;
  gsl_spline *linear_powspec_spline;
  gsl_interp_accel *linear_powspec_acc;
  double _linear_powspec_c0,_linear_powspec_c1; //parameters of power-law extrapolation
  void init_linear_powspec_table(void);
  
 public:
  LinearPowerSpectrum() {
    LINEAR_POWSPEC_TABLE_LENGTH = 1000;
    PL_K_MIN = 1e-6;
    PL_K_MAX = 1e10; 
    LINEAR_POWSPEC_FIT_LENGTH = 20;
    _linear_powspec_norm = 1.0;
    linear_powspec_spline = NULL;
    linear_powspec_acc = NULL;
  };

  ~LinearPowerSpectrum() {
    if(linear_powspec_spline != NULL)
      gsl_spline_free(linear_powspec_spline);
    if(linear_powspec_acc != NULL)
      gsl_interp_accel_free(linear_powspec_acc);
  };

  void init(double sigma8, double spectral_index, class GrowthFunction& gf, class TransferFunction& tf) {
    _s8 = sigma8;
    _ns = spectral_index;
    _gf = &gf;
    _tf = &tf;
    init_linear_powspec_table();
  };
  
  double operator()(double k, double a) {
    double gf = (*_gf)(1.0,a);
    if(k < PL_K_MIN)
      {
        double tf = (*_tf)(k);
        return tf*tf*pow(k,_ns)*_linear_powspec_norm*gf*gf;
      }
    else if(k < PL_K_MAX)
      return exp(gsl_spline_eval(linear_powspec_spline,log(k),linear_powspec_acc))*gf*gf*_linear_powspec_norm;
    else
      return exp(_linear_powspec_c0+_linear_powspec_c1*log(k))*gf*gf*_linear_powspec_norm;
  };
};

#endif /* LINPOWSPEC */
