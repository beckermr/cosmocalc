#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "cosmocalc.h"

#ifndef _W0WACOSMO_
#define _W0WACOSMO_

//cosmological parameters
class w0wa_CosmoData {
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
};

class w0wa_Hubble : public Hubble {
 private:
  class w0wa_CosmoData cd;
  
 public:
  w0wa_Hubble() {};
  w0wa_Hubble(class w0wa_CosmoData& _cd) {cd = _cd;};
  void init(class w0wa_CosmoData& _cd) {cd = _cd;};
  
  class w0wa_CosmoData cosmology(void) {return cd;};
  
  double operator()(double a)
  {
    return sqrt(cd.om/a/a/a + cd.ok/a/a + cd.ol*exp(3.0*(cd.wa*(a-1) - log(a)*(1.0 + cd.w0 + cd.wa))));
  };
};

class w0wa_Distances : public Distances {
 private:
  //cosmology and constants
  double AEXPN_MIN; // minimum expansion factor for standard splines
  double AEXPN_MAX; //maximum expansion factor for standard splines
  class w0wa_CosmoData cd;  
  double DH_sqrtok;
  
  //interpolation table
  int COMVDIST_TABLE_LENGTH; //number of spline points in a for distances
  gsl_spline *aexpn2comvdist_spline;
  gsl_interp_accel *aexpn2comvdist_acc;
  gsl_spline *comvdist2aexpn_spline;
  gsl_interp_accel *comvdist2aexpn_acc;
  void init_comvdist_table(class Hubble& h);
  
  void init_table_vars(void)
  {
    //constants
    AEXPN_MIN = 0.001; // minimum expansion factor for standard splines
    AEXPN_MAX = 1.0; //maximum expansion factor for standard splines
    
    //distances
    COMVDIST_TABLE_LENGTH = 1000;
    aexpn2comvdist_spline = NULL;
    aexpn2comvdist_acc = NULL;
    comvdist2aexpn_spline = NULL;
    comvdist2aexpn_acc = NULL;
  };
  
 public:
  w0wa_Distances() {
    init_table_vars();
  };
  ~w0wa_Distances() {
    if(aexpn2comvdist_spline != NULL)
      gsl_spline_free(aexpn2comvdist_spline);
    if(aexpn2comvdist_acc != NULL)
      gsl_interp_accel_free(aexpn2comvdist_acc);
    if(comvdist2aexpn_spline != NULL)
      gsl_spline_free(comvdist2aexpn_spline);
    if(comvdist2aexpn_acc != NULL)
      gsl_interp_accel_free(comvdist2aexpn_acc);
  };
  w0wa_Distances(class w0wa_CosmoData& _cd, class Hubble& h) {
    init_table_vars();
    cd = _cd;
    DH_sqrtok = DH/sqrt(fabs(cd.ok));
    init_comvdist_table(h);    
  };
  void init(class w0wa_CosmoData& _cd,  class Hubble& h) {
    init_table_vars();
    cd = _cd;
    DH_sqrtok = DH/sqrt(fabs(cd.ok));
    init_comvdist_table(h);
  };
  class w0wa_CosmoData cosmology(void) {return cd;};
  
  //distances
  double comvdist_exact(double a, class Hubble& h);
  double acomvdist(double dist)
  {
    return gsl_spline_eval(comvdist2aexpn_spline,dist,comvdist2aexpn_acc);
  };
  double comvdist(double a)
  {
    return gsl_spline_eval(aexpn2comvdist_spline,a,aexpn2comvdist_acc);
  };
  double angdist(double a)
  {
    if(cd.ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)*a;
    else if(cd.ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)*a;
    else
      return comvdist(a)*a;
  };
  double lumdist(double a)
  {
    if(cd.ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)/a;
    else if(cd.ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)/a;
    else
      return comvdist(a)/a;
  };
  double angdistdiff(double a1, double a2)
  {
    if(a1 < a2)
      {
        if(cd.ok > 0.0)
          return DH_sqrtok*sinh((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
        else if(cd.ok < 0)
          return DH_sqrtok*sin((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
        else
          return (comvdist(a1)-comvdist(a2))*a1;
      }
    else
      {
        if(cd.ok > 0.0)
          return DH_sqrtok*sinh((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
        else if(cd.ok < 0)
          return DH_sqrtok*sin((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
        else
          return (comvdist(a2)-comvdist(a1))*a2;
      }
  };
};


#endif _W0WACOSMO_ /* _W0WACOSMO_ */
