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

/* distances.c - computes distances - assumes flat lambda */
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

/* linear_corrfunc.c */
double linear_corrfunc_exact(double r, double a);
double linear_corrfunc(double r, double a);

/* in nonlinear_powspec.c */
double nonlinear_powspec(double k, double a);

/* nonlinear_corrfunc.c */
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

/* in mass_bias_functions.c */
double tinker2008_mass_function(double m, double a, double delta);
double tinker2010_bias(double m, double a, double delta);
