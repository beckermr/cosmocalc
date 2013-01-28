#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <cstdlib>

#include "cosmocalc.h"

#ifndef _W0WACOSMO_
#define _W0WACOSMO_

//cosmological parameters
class w0wa_CosmoData : public CosmoData {
 public:
  double om;   // total matter density in units of critical at z = 0
  double ol;   // dark energy density in units of critical at z = 0
  double ob;   // baryon density in units of critical at z = 0
  double onu;  // neutrino density in units of critical at z = 0 - NOT CURRENTLY USED
  double ok;   // curvature density in units of critical at z = 0
  double h;    // Hubble constant define as H_0 = h*100 km/s/Mpc
  double s8;   // sigma8 - power spectrum normalization in real space top-hat 8 Mpc/h spheres at z = 0
  double ns;   // spectral index - ns = 1 is scale-invar.
  double w0;   // w0 in w(a) = w0 + (1-a)*wa
  double wa;   // wa in w(a) = w0 + (1-a)*wa
  
  //not required, but allows one to compare cosmo data returned by various objects below using overloaded == operator
  bool operator==(const w0wa_CosmoData& rhs)
  {
    return (om == rhs.om &&
	    ol == rhs.ol &&
	    ob == rhs.ob &&
	    onu == rhs.onu &&
	    ok == rhs.ok &&
	    h == rhs.h &&
	    s8 == rhs.s8 &&
	    ns == rhs.ns &&
	    w0 == rhs.w0 &&
	    wa == rhs.wa);
  };
  
  //required methods
  void print_header(FILE *fp) {
    fprintf(fp,"#omegm omegal omegab omeganu omegak h sigma8 ns w0 wa\n");
  };
  void print_cosmo(FILE *fp) {
    fprintf(fp,"%g %g %g %g %g %g %g %g %g %g\n",om,ol,ob,onu,ok,h,s8,ns,w0,wa);
  };
  double *pack(void) {
    double *data = (double*)malloc(sizeof(double)*10);
    
    data[0] = om;
    data[1] = ol;
    data[2] = ob;
    data[3] = onu;
    data[4] = ok;
    data[5] = h;
    data[6] = s8;
    data[7] = ns;
    data[8] = w0;
    data[9] = wa;
    
    return data;
  }
  void unpack(double *data) {
    om = data[0];
    ol = data[1];
    ob = data[2];
    onu = data[3];
    ok = data[4];
    h = data[5];
    s8 = data[6];
    ns = data[7];
    w0 = data[8];
    wa = data[9];
  };
};

class w0wa_Hubble : public Hubble {
 private:
  class w0wa_CosmoData cd;
  
 public:
  w0wa_Hubble() {};
  void init(class w0wa_CosmoData& _cd) {cd = _cd;};
  
  class w0wa_CosmoData cosmology(void) {return cd;};
  
  double operator()(double a)
  {
    return sqrt(cd.om/a/a/a + cd.ok/a/a + cd.ol*exp(3.0*(cd.wa*(a-1) - log(a)*(1.0 + cd.w0 + cd.wa))));
  };
};

class w0wa_GrowthFunction : public GrowthFunction {
 protected:
  class w0wa_CosmoData cd;
  int GROWTH_FUNCTION_TABLE_LENGTH; //number of spline points for growth function
  double AEXPN_MAX; 
  double AEXPN_MIN; //needs to be in the matter dominated era of your cosmology!
  double LOG_AEXPN_MIN;
  double _growth_function_norm;
  gsl_spline *growth_function_spline;
  gsl_interp_accel *growth_function_acc;
  void init_growth_function_table(class Hubble& h);
  
  void init_table_vars(void)
  {
    //constants
    GROWTH_FUNCTION_TABLE_LENGTH = 50;
    AEXPN_MAX = 1.0;
    AEXPN_MIN = 1.0/31.0;
    LOG_AEXPN_MIN = log(AEXPN_MIN);
    _growth_function_norm = -1.0;
    growth_function_spline = NULL;
    growth_function_acc = NULL;
  }
  
 public:
  w0wa_GrowthFunction() {
    init_table_vars();
  };
  ~w0wa_GrowthFunction() {
    if(growth_function_spline != NULL)
      gsl_spline_free(growth_function_spline);
    if(growth_function_acc != NULL)
      gsl_interp_accel_free(growth_function_acc);
  };
  void init(class w0wa_CosmoData& _cd,  class Hubble& h) {
    init_table_vars();
    cd = _cd;
    init_growth_function_table(h);
  };
  class w0wa_CosmoData cosmology(void) {return cd;};
  
  double operator()(double k, double a) {
    return gsl_spline_eval(growth_function_spline,a,growth_function_acc);
  };
  double growth_function_exact(double k, double a, class Hubble& h);
  double growth_function_norm(void) {return _growth_function_norm;};
};

#endif _W0WACOSMO_ /* _W0WACOSMO_ */
