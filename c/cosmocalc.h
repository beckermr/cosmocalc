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

//constants
#define RHO_CRIT 2.77519737e11 /* Critial mass density  in h^2 M_sun/Mpc^3 */
#define CSOL 299792.458 /* velocity of light in km/s */
#define DELTAC 1.686 /* peak height factor */
#define DH 2997.92458 /* Hubble distance in Mpc/h */
#define TCMB 2.7255 /* EH98 uses 2.728 */
//#define TCMB 2.728 /* EH98 uses 2.728 */

typedef struct {
  int cosmoNum;
  double OmegaM;
  double OmegaB;
  double OmegaL;
  double OmegaK;
  double OmegaNu;
  double h;
  double Sigma8;
  double As;
  double As_pivot;
  double SpectralIndex;
  double delta;
  double w0;
  double wa;
  int useSmoothTransFunc;
} cosmocalcData;

extern cosmocalcData cosmoData;

//lengths and ranges of spline tables
#define AEXPN_MIN 0.001
#define AEXPN_MAX 1.0
#define AEXPN_MIN_GROWTH (1.0/31.0)
#define K_MIN 1e-9
#define K_MAX 1e20
#define PH_R_MIN 1e-3
#define PH_R_MAX 1e3
#define CF_R_MIN 1e-10
#define CF_R_MAX 3e2

#define COSMOCALC_COMVDIST_TABLE_LENGTH            1000
#define COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH     500
#define COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH   5000
#define COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH     20
#define COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH      5000
#define COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH        20
#define COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH   100
#define COSMOCALC_PEAKHEIGHT_TABLE_LENGTH          100
#define COSMOCALC_NONLINEAR_CORRFUNC_TABLE_LENGTH  500
#define COSMOCALC_LINEAR_CORRFUNC_TABLE_LENGTH     500
#define COSMOCALC_LINEAR_POWSPEC_NORM_TABLE_LENGTH 100

#define GSL_SPLINE_TYPE gsl_interp_cspline

/* in utils.c */
double wtime(void);
void gauleg(double x1, double x2, double x[], double w[], int n);

#ifdef FFTLOG
/* in fftlog.c */
void compute_discrete_spherical_fft(double *data, int N, double r0, double L, double q, double *result, double *k0);
double get_k0_fftlog(int N, double mu, double q, double r0, double L, double k0guess);
void um_fftlog(int m, double mu, double q, double k0, double r0, double L, double *realpart, double *imagpart);
#endif

/* distances.c - computes distances - assumes flat lambda */
void init_cosmocalc_comvdist_table(void);
double angdist(double a);
double lumdist(double a);
double comvdist(double a);
double angdistdiff(double amin, double amax);
double acomvdist(double dist);
double comvdist_exact(double a);

/* in transfer_function.c */
double transfer_function(double k);
double transfunct_eh98(double kin);
double transfunct_eh98_smooth(double kin);

/* in linear_powspec.c */
double linear_powspec(double k, double a);
double linear_powspec_exact(double k, double a);
double tophatradnorm_linear_powspec_exact_nonorm(double topHatRad);
double tophatnorm_linear_powspec(double topHatRad);
double get_linear_tophatnorm_scale(double a);
double convert_cmbnorm2sigma8(void);

/* linear_corrfunc.c */
double linear_corrfunc_exact(double r, double a);
double linear_corrfunc(double r, double a);

/* in nonlinear_powspec.c */
double nonlinear_powspec(double k, double a);
double get_nonlinear_gaussnorm_scale(double a);
double gaussiannorm_linear_powspec(double gaussRad);

/* nonlinear_corrfunc.c */
double nonlinear_corrfunc_integ_funct(double k, void *p);
double nonlinear_corrfunc_exact(double r, double a);
double nonlinear_corrfunc(double r, double a);

/* in growth_function.c */
double growth_function_exact(double a);
double growth_function(double a);
double growth_function_exact_nonorm(double a);

/* in hubble.c */
double weff(double a);
double hubble_noscale(double a);

/* in peakheight.c */
double sigmaRtophat_exact(double topHatRad, double a);
double sigmaRtophat(double topHatRad, double a);
double sigmaMtophat_exact(double m, double a);
double sigmaMtophat(double m, double a);
double inverse_sigmaRtophat(double sigmaR, double a);
double inverse_sigmaMtophat(double sigmaR, double a);
double inverse_nuMtophat(double nu, double a);
double inverse_nuRtophat(double nu, double a);

/* in mass_bias_functions.c */
double mass_function(double m, double a);
double bias_function(double m, double a);
double tinker2010_mass_function(double m, double a, double delta);
double tinker2008_mass_function(double m, double a, double delta);

#endif /* _COSMOCALC_ */
