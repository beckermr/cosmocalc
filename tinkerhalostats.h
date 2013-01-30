#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include "halomodel.h"
#include "peakheight.h"

#ifndef _TINKERHALOSTATS_
#define _TINKERHALOSTATS_

class Tinker08MassFunction : public MassFunction {
 private:
  gsl_spline *ATm_spline,*aTm_spline,*bTm_spline,*cTm_spline;
  gsl_interp_accel *ATm_acc,*aTm_acc,*bTm_acc,*cTm_acc;
  class PeakHeight *_ph;
  double _om,_delta;
  double A,am,b,c,alpha;
  
 public:
  Tinker08MassFunction() {
    double deltaTm[9] = {200.0, 300.0, 400.0, 600.0, 800.0, 1200.0, 1600.0, 2400.0, 3200.0};
    double ATm[9] = {0.186, 0.2, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260};
    double aTm[9] = {1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66};
    double bTm[9] = {2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41};
    double cTm[9] = {1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44};

    ATm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
    aTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
    bTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
    cTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
    
    gsl_spline_init(ATm_spline,deltaTm,ATm,9);
    gsl_spline_init(aTm_spline,deltaTm,aTm,9);
    gsl_spline_init(bTm_spline,deltaTm,bTm,9);
    gsl_spline_init(cTm_spline,deltaTm,cTm,9);
    
    ATm_acc = gsl_interp_accel_alloc();
    aTm_acc = gsl_interp_accel_alloc();
    bTm_acc = gsl_interp_accel_alloc();
    cTm_acc = gsl_interp_accel_alloc();
  };
  
  ~Tinker08MassFunction() {
    gsl_spline_free(ATm_spline);
    gsl_spline_free(aTm_spline);
    gsl_spline_free(bTm_spline);
    gsl_spline_free(cTm_spline);
    
    gsl_interp_accel_free(ATm_acc);
    gsl_interp_accel_free(aTm_acc);
    gsl_interp_accel_free(bTm_acc);
    gsl_interp_accel_free(cTm_acc);
  };
  
  void init(double omegam, double delta, class PeakHeight& ph) {
    _om = omegam;
    _ph = &ph;
    _delta = delta;
    
    A = gsl_spline_eval(ATm_spline,delta,ATm_acc);
    am = gsl_spline_eval(aTm_spline,delta,aTm_acc);
    b = gsl_spline_eval(bTm_spline,delta,bTm_acc);
    c = gsl_spline_eval(cTm_spline,delta,cTm_acc);
    alpha = pow(10.0,-1.0*pow(0.75/log10(_delta/75.0),1.2));
  };
  
  double operator()(double m, double a) {
    double Aa = A*pow(1.0/a,-0.14);
    double ama = am*pow(1.0/a,-0.06);
    double ba = b*pow(1.0/a,-1.0*alpha);
  
    //fprintf(stderr,"A = %g, a = %g, b = %g, c = %g\n",A,am,b,c);
  
    double sigma = _ph->sigmaMtophat(m,a);
    double fsigma = Aa*(pow(sigma/ba,-1.0*ama) + 1.0)*exp(-1.0*c/sigma/sigma);
    double dm = 1e-6*m;
    double dlnsiginvdm = log(_ph->sigmaMtophat(m-dm/2.0,a)/_ph->sigmaMtophat(m+dm/2.0,a))/dm;
  
    return fsigma*RHO_CRIT*_om/m*dlnsiginvdm;
  };
};

class TinkerHaloBias : public HaloBias {
 private:
  double _delta;
  double y;
  double A;
  double a;
  double B;
  double b;
  double C;
  double c;
  class PeakHeight *_ph;
  
 public:
  TinkerHaloBias() {};
  void init(double delta, class PeakHeight& ph) {
    _ph = &ph;
    _delta = delta;
    y = log10(delta);
    A = 1.0 + 0.24*y*exp(-1.0*pow(4.0/y,4.0));
    a = 0.44*y - 0.88;
    B = 0.183;
    b = 1.5;
    C = 0.019 + 0.107*y + 0.19*exp(-1.0*pow(4.0/y,4.0));
    c = 2.4;
  };
  
  double operator()(double m, double a) {
    double nu = DELTAC/_ph->sigmaMtophat(m,a);
    return 1.0 - A*pow(nu,a)/(pow(nu,a) + pow(1.686,a)) + B*pow(nu,b) + C*pow(nu,c);
  };
};

#endif /* _TINKERHALOSTATS_ */
