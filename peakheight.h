#include "cosmocalc_assert.h"
#include "cosmocalc.h"

#ifndef _PEAKHEIGHT_
#define _PEAKHEIGHT_

//standard physical constants
#define DELTAC 1.686

class PeakHeight {
 private:
  int PEAKHEIGHT_TABLE_LENGTH; //number of spline points in radius for getting peak height
  double PH_R_MIN; //min scale for peak height table
  double PH_R_MAX; //max scale for peak height table
  gsl_spline *R2sigma_spline;
  gsl_interp_accel *R2sigma_acc;
  gsl_spline *sigma2R_spline;
  gsl_interp_accel *sigma2R_acc;
  void init_peakheight_table(void);
  class GrowthFunction *_gf;
  double _om;
  class PowerSpectrum *_pkl;
  
 public:
  PeakHeight() {
    PEAKHEIGHT_TABLE_LENGTH = 50;
    PH_R_MIN = 1e-2;
    PH_R_MAX = 50.0;
    R2sigma_spline = NULL;
    R2sigma_acc = NULL;
    sigma2R_spline = NULL;
    sigma2R_acc = NULL;
  };
  ~PeakHeight() {
    if(R2sigma_spline != NULL)
      gsl_spline_free(R2sigma_spline);
    if(R2sigma_acc != NULL)
      gsl_interp_accel_free(R2sigma_acc);
    if(sigma2R_spline != NULL)
      gsl_spline_free(sigma2R_spline);
    if(sigma2R_acc != NULL)
      gsl_interp_accel_free(sigma2R_acc);
  };
  
  void init(double omegam, class GrowthFunction& gf, class PowerSpectrum& pkl) {
    _om = omegam;
    _gf = &gf;
    _pkl = &pkl;
  };
  
  double sigmaRtophat(double topHatRad, double a)
  {
    return exp(gsl_spline_eval(R2sigma_spline,log(topHatRad),R2sigma_acc))*(*_gf)(a);
  };
  double sigmaMtophat(double m, double a)
  {
    return sigmaRtophat(pow(m/(4.0/3.0*M_PI*RHO_CRIT*_om),1.0/3.0),a);
  };
  double Rsigmatophat(double sigmaR, double a)
  {
    return exp(gsl_spline_eval(sigma2R_spline,log(sigmaR/(*_gf)(a)),sigma2R_acc));
  };
  double Msigmatophat(double sigmaR, double a)
  {
    return (4.0/3.0*M_PI*RHO_CRIT*_om)*pow(Rsigmatophat(sigmaR,a),3.0);
  };
  double Mnutophat(double nu, double a)
  {
    return (4.0/3.0*M_PI*RHO_CRIT*_om)*pow(Rsigmatophat(DELTAC/nu,a),3.0);
  };
  double Rnutophat(double nu, double a)
  {
    return Rsigmatophat(DELTAC/nu,a);
  };
};

#endif /* _PEAKHEIGHT_ */
