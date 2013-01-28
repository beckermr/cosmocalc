#include <cstdio>
#include "cosmocalc_assert.h"

#ifndef _COSMOCALC_
#define _COSMOCALC_

//standard physical constants
#define RHO_CRIT 2.77519737e11 // Critial mass density  in h^2 M_sun/Mpc^3
#define CSOL 299792.458 // velocity of light in km/s
#define DH 2997.92458 // Hubble Distance assuming H0 = 100 km/s

/*CosmoData base class
  This base class just defines the interface. You must define your own class 
  to actually use it. See the w0wacosmo.h header file for details.
*/
class CosmoData {
 public:
  //public methods to output parameters to a file
  virtual void print_header(FILE *fp) {
    cosmocalc_assert(false,"In print_header method, CosmoData base class cannot be used alone!");
  };
  virtual void print_cosmo(FILE *fp) {
    cosmocalc_assert(false,"In print_cosmo method, CosmoData base class cannot be used alone!");
  };
  
  //methods to pack and unpack parameters for sending over MPI, binary I/O, etc.
  virtual double *pack(void) {
    cosmocalc_assert(false,"In pack method, CosmoData base class cannot be used alone!");
  };
  virtual void unpack(double *data) {
    cosmocalc_assert(false,"In unpack method, CosmoData base class cannot be used alone!");
  };    
};

/*Hubble factor base class 
  This base class just defines the interface. You must define your own class 
  to actually use it. See the w0wacosmo.h header file for details.
*/
class Hubble {
 public:
  virtual double operator()(double a) {
    cosmocalc_assert(false,"In operator() method, Hubble base class cannot be used alone!");
  };
};

/*Distances base class 
  This base class just defines the interface. You must define your own class 
  to actually use it. See the w0wacosmo.h header file for details.
*/
class Distances {
 public:
  virtual double comvdist_exact(double a, class Hubble& h) {
    cosmocalc_assert(false,"In comvdist_exact method, Distances base class cannot be used alone!");
  };
  virtual double acomvdist(double dist) {
    cosmocalc_assert(false,"In acomvdist method, Distances base class cannot be used alone!");
  };
  virtual double comvdist(double a) {
    cosmocalc_assert(false,"In comvdist method, Distances base class cannot be used alone!");
  };
  virtual double angdist(double a) {
    cosmocalc_assert(false,"In angdist method, Distances base class cannot be used alone!");
  };
  virtual double lumdist(double a)
  {
    cosmocalc_assert(false,"In lumdist method, Distances base class cannot be used alone!");
  };
  virtual double angdistdiff(double a1, double a2) {
    cosmocalc_assert(false,"In angdistdiff method, Distances base class cannot be used alone!");
  };
};

/*GrowthFunction base class.
  This base class just defines the interface. You must define your own class
  to actually use it. See the w0wacosmo.h header file for details.   
*/
class GrowthFunction {
 public:
  virtual double operator()(double k, double a) {
    cosmocalc_assert(false,"In operator() method, GrowthFunction base class cannot be used alone!");
  };
};

#endif /* _COSMOCALC_ */
