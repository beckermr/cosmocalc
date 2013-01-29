#include "cosmocalc.h"

#ifndef _EH98_TRANSFUNCT_
#define _EH98_TRANSFUNCT_

class EH98_TransferFunction : public TransferFunction {
 protected:
  double _om,_ob,_h;
  double eh98_transfer_function(double k);
  
 public:
  EH98_TransferFunction() { };
  void init(double omegam, double omegab, double h) {
    _om = omegam; //total matter density in units of critical at z = 0
    _ob = omegab; //baryon "                                          "
    _h = h;       // little h to convert k in h/Mpc to 1/Mpc units
  };
  
  double operator()(double k) {
    return eh98_transfer_function(k);
  };
};

class EH98Smooth_TransferFunction : public TransferFunction {
 protected:
  double _om,_ob,_h;
  double eh98smooth_transfer_function(double k);
  
 public:
  EH98Smooth_TransferFunction() { };
  void init(double omegam, double omegab, double h) {
    _om = omegam; //total matter density in units of critical at z = 0
    _ob = omegab; //baryon "                                          "
    _h = h;       // little h to convert k in h/Mpc to 1/Mpc units
  };
  
  double operator()(double k) {
    return eh98smooth_transfer_function(k);
  };
};

#endif /* EH98_TRANSFUNCT */
