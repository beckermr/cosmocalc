#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>

#include "cosmocalc.h"

#ifndef _FLRW_DISTANCES_
#define _FLRW_DISTANCES_

class FLRWDistances : public Distances {
 private:
  //cosmology and constants
  double AEXPN_MIN; // minimum expansion factor for standard splines
  double AEXPN_MAX; //maximum expansion factor for standard splines
  double ok;
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
  FLRWDistances() {
    init_table_vars();
  };
  ~FLRWDistances() {
    if(aexpn2comvdist_spline != NULL)
      gsl_spline_free(aexpn2comvdist_spline);
    if(aexpn2comvdist_acc != NULL)
      gsl_interp_accel_free(aexpn2comvdist_acc);
    if(comvdist2aexpn_spline != NULL)
      gsl_spline_free(comvdist2aexpn_spline);
    if(comvdist2aexpn_acc != NULL)
      gsl_interp_accel_free(comvdist2aexpn_acc);
  };
  void init(double _ok,  class Hubble& h) {
    init_table_vars();
    ok = _ok;
    DH_sqrtok = DH/sqrt(fabs(ok));
    init_comvdist_table(h);
  };
    
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
    if(ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)*a;
    else if(ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)*a;
    else
      return comvdist(a)*a;
  };
  double lumdist(double a)
  {
    if(ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)/a;
    else if(ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)/a;
    else
      return comvdist(a)/a;
  };
  double angdistdiff(double a1, double a2)
  {
    if(a1 < a2)
      {
        if(ok > 0.0)
          return DH_sqrtok*sinh((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
        else if(ok < 0)
          return DH_sqrtok*sin((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
        else
          return (comvdist(a1)-comvdist(a2))*a1;
      }
    else
      {
        if(ok > 0.0)
          return DH_sqrtok*sinh((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
        else if(ok < 0)
          return DH_sqrtok*sin((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
        else
          return (comvdist(a2)-comvdist(a1))*a2;
      }
  };
};

#endif _FLRW_DISTANCES_ /* _FLRW_DISTANCES_ */
