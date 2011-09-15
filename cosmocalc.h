#include <stdio.h>

#ifndef _COSMOCALC_
#define _COSMOCALC_

/* Some debugging macros
   undef DEBUG for no debugging
   DEBUG_LEVEL = 0 is for basic debugging
   DEBUG_LEVEL = 1 is for messages printed by a single task but not as critical
   DEBUG_LEVEL = 2 or above is used for messages that every task will print 
   
   define DEBUG_IO some output as follows
*/

#ifdef NDEBUG
#undef DEBUG
#define DEBUG_LEVEL -1
#undef DEBUG_IO
#endif

#ifdef DEBUG
#ifndef DEBUG_LEVEL
#define DEBUG_LEVEL 0
#endif
#endif

#define MAX_FILENAME 1024

#define RHO_CRIT 2.77519737e11 /* Critial mass density  in h^2 M_sun/Mpc^3 */
#define CSOL 299792.458 /* velocity of light in km/s */
#define DELTAC 1.686 /* peak height factor */

typedef struct {
  int cosmoNum;
  double OmegaM;
  double h;
  double Sigma8;
  double SpectralIndex;
  double OmegaB;
  double delta;
} cosmocalcData;

extern cosmocalcData cosmoData;

#ifdef TEST_CODE
void test_nonlinear_corrfunc(void);
void test_nonlinear_powspec(void);
void test_linxi_fftlog(void);
void test_linxi(void);
void test_biasfunc(void);
void test_massfunc(void);
void test_peakheight(void);
void test_gf(void);
void test_distances(void);
void test_transfunct(void);
void test_linpk(void);
#endif

/* in utils.c */
double wtime(void);
void gauleg(double x1, double x2, double x[], double w[], int n);

/* in fftlog.c */
void compute_discrete_spherical_fft(double *data, int N, double r0, double L, double q, double *result, double *k0);
double get_k0_fftlog(int N, double mu, double q, double r0, double L, double k0guess);
void um_fftlog(int m, double mu, double q, double k0, double r0, double L, double *realpart, double *imagpart);

/* distances.c - computes distances - assumes flat lambda */
void init_cosmocalc_comvdist_table(void);
double comvdist_integ_funct(double a, void *p);
double angdist(double a);
double comvdist(double z);
double angdistdiff(double amin, double amax);
double acomvdist(double dist);
double comvdist_exact(double a);

/* in transfer_function.c */
double transfer_function(double k);
double transfunct_eh98(double kin);

/* in linear_powspec.c */
double linear_powspec(double k, double a);
double linear_powspec_exact(double k, double a);
double tophatradnorm_linear_powspec_exact_nonorm(double topHatRad);
double tophatradnorm_linear_powspec_exact_nonorm_k_integ_funct_I0(double k, void *p);
double tophatradnorm_linear_powspec_exact_nonorm_lnk_integ_funct_I0(double lnk, void *p);
double fourierTransformTopHat(double y);

/* linear_corrfunc.c */
double linear_corrfunc_integ_funct(double k, void *p);
double linear_corrfunc_exact(double r, double a);
double linear_corrfunc(double r, double a);

/* in nonlinear_powspec.c */
double gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
double gaussiannorm_linear_powspec_exact(double gaussRad);
double onederiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
double onederiv_gaussiannorm_linear_powspec_exact(double gaussRad);
double twoderiv_gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p);
double twoderiv_gaussiannorm_linear_powspec_exact(double gaussRad);
double nonlinear_gaussnorm_scale_funct(double gaussR, void *p);
double get_nonlinear_gaussnorm_scale(double a);
double nonlinear_powspec_exact(double k, double a);
double nonlinear_powspec(double k, double a);

/* nonlinear_corrfunc.c */
double nonlinear_corrfunc_integ_funct(double k, void *p);
double nonlinear_corrfunc_exact(double r, double a);
double nonlinear_corrfunc(double r, double a);

/* in growth_function.c */
double growth_function_integ_funct(double a, void *p);
double growth_function_exact(double a);
double growth_function(double a);

/* in hubble.c */
double hubble_noscale(double a);

/* in peakheight.c */
double sigmaRtophat_exact(double topHatRad, double a);
double sigmaRtophat(double topHatRad, double a);
double sigmaMtophat(double m, double a);
double inverse_sigmaRtophat(double sigmaR, double a);
double inverse_sigmaMtophat(double sigmaR, double a);
double inverse_nuMtophat(double nu, double a);
double inverse_nuRtophat(double nu, double a);

/* in mass_bias_functions.c */
double mass_function(double m, double a);
double bias_function(double m, double a);

#endif /* _COSMOCALC_ */
