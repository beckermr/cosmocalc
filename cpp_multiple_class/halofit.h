#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#include "cosmocalc.h"

#ifndef _HALOFITPOWSPEC_
#define _HALOFITPOWSPEC_

class HaloFitPowerSpectrum : public PowerSpectrum {
 protected:
  double _om,_ol;
  class GrowthFunction *_gf;
  class PowerSpectrum *_pklin;
  class Hubble *_ha;
  class CosmoData *_cd;
  
  int NONLINEAR_POWSPEC_TABLE_LENGTH; //number of spline points in scale factor for nonlinear powspec gaussn norm table
  double PNL_A_MIN; //min scale factor for nonlinear powspec gaussn norm table 
  double PNL_A_MAX; //max scale factor for nonlinear powspec gaussn norm table 
  double get_nonlinear_gaussnorm_scale(double a);
  void init_nonlinear_powspec_table(void);
  int NUM_SPLINE_NONLINEAR_POWSPEC;
  gsl_spline *nonlinear_powspec_spline[5];
  gsl_interp_accel *nonlinear_powspec_acc[5];
  double PNL_RGAUSS_MIN;
  double PNL_RGAUSS_MAX;
    
 public:
  HaloFitPowerSpectrum() {
    NONLINEAR_POWSPEC_TABLE_LENGTH = 100;
    NUM_SPLINE_NONLINEAR_POWSPEC = 5;
    PNL_A_MIN = 0.2;
    PNL_A_MAX = 1.0;
    PNL_RGAUSS_MIN = 0.0001;
    PNL_RGAUSS_MAX = 100.0;
    for(int i=0;i<NUM_SPLINE_NONLINEAR_POWSPEC;++i)
      {
	nonlinear_powspec_spline[i] = NULL;
	nonlinear_powspec_acc[i] = NULL;
      }
  };

  ~HaloFitPowerSpectrum() {
    for(int i=0;i<NUM_SPLINE_NONLINEAR_POWSPEC;++i)
      {
	if(nonlinear_powspec_spline[i] != NULL)
	  gsl_spline_free(nonlinear_powspec_spline[i]);
	if(nonlinear_powspec_acc[i] != NULL)
	  gsl_interp_accel_free(nonlinear_powspec_acc[i]);
      }
  };

  void init(double omegam, double omegal, class CosmoData& cd, class PowerSpectrum& pklin, class Hubble& ha, class GrowthFunction& gf) {
    _om = omegam;
    _ol = omegal;
    _pklin = &pklin;
    _ha = &ha;
    _cd = &cd;
    _gf = &gf;
    init_nonlinear_powspec_table();
  };
  
  double gaussiannorm_linear_powspec_exact(double gaussRad);
  double gaussiannorm_linear_powspec(double gaussRad) {
    return exp(gsl_spline_eval(nonlinear_powspec_spline[2],log(gaussRad),nonlinear_powspec_acc[2]));
  };
  double linear_powspec(double k, double a) {
    return (*_pklin)(k,a);
  };
  
  double operator()(double k, double a);
};

#endif /* _HALOFITPOWSPEC_ */
