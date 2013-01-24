#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>

#ifndef _COSMOCALC_
#define _COSMOCALC_

#ifdef _OPENMP
#include <omp.h>
#endif

class cosmoCalc {
 private:
  //constants
  double RHO_CRIT;  // Critial mass density  in h^2 M_sun/Mpc^3
  double CSOL;      // velocity of light in km/s
  double DH;        // Hubble Distance assuming H0 = 100 km/s
  double DELTAC;    // peak height factor
  double AEXPN_MIN; // minimum expansion factor for standard splines
  double AEXPN_MAX; //maximum expansion factor for standard splines
  
  //cosmological parameters
  double _om;   // total matter density in units of critical at z = 0
  double _ol;   // dark energy density in units of critical at z = 0
  double _ob;   // baryon density in units of critical at z = 0
  double _onu;  // neutrino density in units of critical at z = 0 - NOT CURRENTLY USED
  double _ok;   // curvature density in units of critical at z = 0
  double _h;    // Hubble constant define as H_0 = h*100 km/s/Mpc
  double _s8;   // sigma8 - power spectrum normalization in real space top-hat 8 Mpc/h spheres at z = 0
  double _ns;   // spectral index - ns = 1 is scale-invar. 
  double _w0;   // w0 in w(a) = w0 + (1-a)*wa
  double _wa;   // wa in w(a) = w0 + (1-a)*wa
  
  double _delta; // overdensity for mass function
  double DH_sqrtok; //Hubble distance / sqrt(Omega_k)
  
  long _cosmo_num; //used to track states of various interpolation tables - if table states (i.e. _dist_cosmo_num) differ
                   // from the global state (_cosmo_num), then tables are remade
  
  int _num_threads; //sets number of threads if wanted
  
  //distances
  int COSMOCALC_COMVDIST_TABLE_LENGTH; //number of spline points in a for distances
  gsl_spline *cosmocalc_aexpn2comvdist_spline;
  gsl_interp_accel *cosmocalc_aexpn2comvdist_acc;
  gsl_spline *cosmocalc_comvdist2aexpn_spline;
  gsl_interp_accel *cosmocalc_comvdist2aexpn_acc;
  void init_cosmocalc_comvdist_table(void);
  long _dist_cosmo_num;
  
  //linear growth function
  int COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH; //number of spline points for growth function
  double AEXPN_MIN_GROWTH; //needs to be in the matter dominated era of your cosmology!
  double LOG_AEXPN_MIN_GROWTH;
  double _growth_function_norm;
  gsl_spline *cosmocalc_growth_function_spline;
  gsl_interp_accel *cosmocalc_growth_function_acc;
  void init_cosmocalc_growth_function_table(void);
  long _gf_cosmo_num;
  
  //transfer_function
  int COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH; //number of spline points in k for transfer function
  double TF_K_MIN; //min k value for transfer function in h/Mpc
  double TF_K_MAX; //max k value for transfer function in h/Mpc
  int COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH; //number of points right before TF_K_MAX to use to fit tf(k) for power-law extrapolation
  gsl_spline *cosmocalc_transfer_function_spline;
  gsl_interp_accel *cosmocalc_transfer_function_acc;
  double _transfer_function_c0,_transfer_function_c1; //parameters of power law extrapolation
  void init_cosmocalc_transfer_function_table(void);
  int _transfer_function_type; //type of transfer function, define through macros below
#define COSMOCALC_TRANS_FUNC_EH98        0  
#define COSMOCALC_TRANS_FUNC_EH98_SMOOTH 1  
  double (cosmoCalc::*_transfer_function_function)(double);
  long _tf_cosmo_num;
  
  //linear powspec
  int COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH; //number of spline points in k for linear power spec
  double PL_K_MIN;  //min k value for linear powspec in h/Mpc
  double PL_K_MAX;  //max k value for linear powspec in h/Mpc
  int COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH; //number of points right before PL_K_MAX to use to fir pl(k) for power-law extrapolation
  double tophatradnorm_linear_powspec_exact_nonorm(double topHatRad); //returns value of unormalized linear powspec filtered by top hat of radius topHatRad
  double _linear_powspec_norm;
  gsl_spline *cosmocalc_linear_powspec_spline;
  gsl_interp_accel *cosmocalc_linear_powspec_acc;
  double _linear_powspec_c0,_linear_powspec_c1; //parameters of power-law extrapolation
  void init_cosmocalc_linear_powspec_table(void);
  long _pl_cosmo_num;
  
  //nonlinear powspec
  int COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH; //number of spline points in scale factor for nonlinear powspec gaussn norm table
  double PNL_A_MIN; //min scale factor for nonlinear powspec gaussn norm table 
  double PNL_A_MAX; //max scale factor for nonlinear powspec gaussn norm table 
  double get_nonlinear_gaussnorm_scale(double a);
  void init_cosmocalc_nonlinear_powspec_table(void);
  gsl_spline *cosmocalc_nonlinear_powspec_spline[3];
  gsl_interp_accel *cosmocalc_nonlinear_powspec_acc[3];
  long _pnl_cosmo_num;
  double PNL_RGAUSS_MIN;
  double PNL_RGAUSS_MAX;
  
  //peak heights
  int COSMOCALC_PEAKHEIGHT_TABLE_LENGTH; //number of spline points in radius for getting peak height
  double PH_R_MIN; //min scale for peak height table
  double PH_R_MAX; //max scale for peak height table
  int _ph_cosmo_num;
  gsl_spline *cosmocalc_R2sigma_spline;
  gsl_interp_accel *cosmocalc_R2sigma_acc;
  gsl_spline *cosmocalc_sigma2R_spline;
  gsl_interp_accel *cosmocalc_sigma2R_acc;
  void init_cosmocalc_peakheight_table(void);

 public:
  //con- and de-structors
  cosmoCalc ()
    {
      long i;
      
      //constants
      RHO_CRIT = 2.77519737e11; // Critial mass density  in h^2 M_sun/Mpc^3
      CSOL =  299792.458; // velocity of light in km/s
      DH =  2997.92458; // Hubble Distance assuming H0 = 100 km/s
      DELTAC = 1.686; // peak height factor
      AEXPN_MIN = 0.001; // minimum expansion factor for standard splines
      AEXPN_MAX = 1.0; //maximum expansion factor for standard splines
      
      _cosmo_num = -1;
      _num_threads = -1;
      
      //distances
      COSMOCALC_COMVDIST_TABLE_LENGTH = 1000;
      cosmocalc_aexpn2comvdist_spline = NULL;
      cosmocalc_aexpn2comvdist_acc = NULL;
      cosmocalc_comvdist2aexpn_spline = NULL;
      cosmocalc_comvdist2aexpn_acc = NULL;
      _dist_cosmo_num = -1;
      
      //linear growth function
      COSMOCALC_GROWTH_FUNCTION_TABLE_LENGTH = 50;
      AEXPN_MIN_GROWTH = 1.0/31.0;
      LOG_AEXPN_MIN_GROWTH = log(AEXPN_MIN_GROWTH);
      _growth_function_norm = -1.0;
      cosmocalc_growth_function_spline = NULL;
      cosmocalc_growth_function_acc = NULL;
      _gf_cosmo_num = -1;
      
      //transfer function 
      COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH = 5000;
      TF_K_MIN = 1e-7;
      TF_K_MAX = 1e20;
      COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH = 20;
      cosmocalc_transfer_function_spline = NULL;
      cosmocalc_transfer_function_acc = NULL;
      _transfer_function_type = -1;
      _transfer_function_function = NULL;
      _tf_cosmo_num = -1;
      
      //linear powspec
      COSMOCALC_LINEAR_POWSPEC_TABLE_LENGTH = 1000;
      PL_K_MIN = 1e-9;
      PL_K_MAX = 1e20; 
      COSMOCALC_LINEAR_POWSPEC_FIT_LENGTH = 20;
      _linear_powspec_norm = -1;
      cosmocalc_linear_powspec_spline = NULL;
      cosmocalc_linear_powspec_acc = NULL;
      _pl_cosmo_num = -1;
      
      //nonlinear powspec
      COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH = 100;
      PNL_A_MIN = 0.2;
      PNL_A_MAX = 1.0;
      PNL_RGAUSS_MIN = 0.0001;
      PNL_RGAUSS_MAX = 100.0;
      for(i=0;i<3;++i)
	{
	  cosmocalc_nonlinear_powspec_spline[i] = NULL;
	  cosmocalc_nonlinear_powspec_acc[i] = NULL;
	}
      _pnl_cosmo_num = -1;
      
      //peak heights
      COSMOCALC_PEAKHEIGHT_TABLE_LENGTH = 50;
      PH_R_MIN = 1e-2;
      PH_R_MAX = 50.0;
      cosmocalc_R2sigma_spline = NULL;
      cosmocalc_R2sigma_acc = NULL;
      cosmocalc_sigma2R_spline = NULL;
      cosmocalc_sigma2R_acc = NULL;
      _ph_cosmo_num = -1;
    };
  
  ~cosmoCalc () 
    {
      long i;
      
      //distances
      if(cosmocalc_aexpn2comvdist_spline != NULL)
	gsl_spline_free(cosmocalc_aexpn2comvdist_spline);
      if(cosmocalc_aexpn2comvdist_acc != NULL)
	gsl_interp_accel_free(cosmocalc_aexpn2comvdist_acc);
      if(cosmocalc_comvdist2aexpn_spline != NULL)
	gsl_spline_free(cosmocalc_comvdist2aexpn_spline);
      if(cosmocalc_comvdist2aexpn_acc != NULL)
	gsl_interp_accel_free(cosmocalc_comvdist2aexpn_acc);
      
      //linear growth function
      if(cosmocalc_growth_function_spline != NULL)
        gsl_spline_free(cosmocalc_growth_function_spline);
      if(cosmocalc_growth_function_acc != NULL)
        gsl_interp_accel_free(cosmocalc_growth_function_acc);
      
      //transfer function 
      if(cosmocalc_transfer_function_spline != NULL)
	gsl_spline_free(cosmocalc_transfer_function_spline);
      if(cosmocalc_transfer_function_acc != NULL)
	gsl_interp_accel_free(cosmocalc_transfer_function_acc);
      
      //linear powspec
      if(cosmocalc_linear_powspec_spline != NULL)
	gsl_spline_free(cosmocalc_linear_powspec_spline);
      if(cosmocalc_linear_powspec_acc != NULL)
	gsl_interp_accel_free(cosmocalc_linear_powspec_acc);
      
      //nonlinear powspec
      for(i=0;i<3;++i)
	{
	  if(cosmocalc_nonlinear_powspec_spline[i] != NULL)
	    gsl_spline_free(cosmocalc_nonlinear_powspec_spline[i]);
	  if(cosmocalc_nonlinear_powspec_acc[i] != NULL)
	    gsl_interp_accel_free(cosmocalc_nonlinear_powspec_acc[i]);
	}
      
      //peak heights
      if(cosmocalc_R2sigma_spline != NULL)
	gsl_spline_free(cosmocalc_R2sigma_spline);
      if(cosmocalc_R2sigma_acc != NULL)
	gsl_interp_accel_free(cosmocalc_R2sigma_acc);
      if(cosmocalc_sigma2R_spline != NULL)
	gsl_spline_free(cosmocalc_sigma2R_spline);
      if(cosmocalc_sigma2R_acc != NULL)
	gsl_interp_accel_free(cosmocalc_sigma2R_acc);
    };

  //deal with threads
  void num_threads(int num_threads) {_num_threads = num_threads;};
  
  //function to init params
  void init_cosmology(double omegam, double omegal, double omegab, double omeganu, 
		      double hubble, double sigma8, double spectral_index, 
		      double w0, double wa, int transfer_function_type);
  
  //functions to init interpolation tables - make sure to catch dependencies
  void init_distances(void);
  void init_growth_function(void);
  void init_transfer_function(void);
  void init_linear_powspec(void);
  void init_nonlinear_powspec(void);
  void init_peakheight(void);
  void init_all(void);
  
  //parameters
  double omegam(void) {return _om;};
  double omegal(void) {return _ol;};
  double omegab(void) {return _ob;};
  double omeganu(void) {return _onu;};
  double omegak(void) {return _ok;};
  double h(void) {return _h;};
  double sigma8(void) {return _s8;};
  double ns(void) {return _ns;};
  double w0(void) {return _w0;};
  double wa(void) {return _wa;};
  
  //distances
  double comvdist_exact(double a);
  double acomvdist(double dist)
  {
    return gsl_spline_eval(cosmocalc_comvdist2aexpn_spline,dist,cosmocalc_comvdist2aexpn_acc);
  };
  double comvdist(double a)
  {
    return gsl_spline_eval(cosmocalc_aexpn2comvdist_spline,a,cosmocalc_aexpn2comvdist_acc);
  };
  double angdist(double a)
  {
    if(_ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)*a;
    else if(_ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)*a;
    else
      return comvdist(a)*a;
  };
  double lumdist(double a)
  {
    if(_ok > 0.0)
      return DH_sqrtok*sinh(comvdist(a)/DH_sqrtok)/a;
    else if(_ok < 0)
      return DH_sqrtok*sin(comvdist(a)/DH_sqrtok)/a;
    else
      return comvdist(a)/a;
  };
  double angdistdiff(double a1, double a2)
  {
    if(a1 < a2)
      {
	if(_ok > 0.0)
	  return DH_sqrtok*sinh((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
	else if(_ok < 0)
	  return DH_sqrtok*sin((comvdist(a1)-comvdist(a2))/DH_sqrtok)/a1;
	else
	  return (comvdist(a1)-comvdist(a2))*a1;
      }
    else
      {
	if(_ok > 0.0)
	  return DH_sqrtok*sinh((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
	else if(_ok < 0)
	  return DH_sqrtok*sin((comvdist(a2)-comvdist(a1))/DH_sqrtok)/a2;
	else
	  return (comvdist(a2)-comvdist(a1))*a2;
      }
  };
  
  //hubble
  double hubble_noscale(double a)
  {
    return sqrt(_om/a/a/a + _ok/a/a + _ol*exp(3.0*(_wa*(a-1) - log(a)*(1.0 + _w0 + _wa))));
  };
  
  //linear growth function
  double growth_function_exact(double a);
  double growth_function(double a)
  {
    return gsl_spline_eval(cosmocalc_growth_function_spline,a,cosmocalc_growth_function_acc)/_growth_function_norm;
  };
  double growth_function_norm(void) {return _growth_function_norm;};

  //transfer_function 
  double transfer_function(double k)
  {
    if(k < TF_K_MIN)
      return 1.0;
    else if(k < TF_K_MAX)
      return exp(gsl_spline_eval(cosmocalc_transfer_function_spline,log(k),cosmocalc_transfer_function_acc));
    else
      return exp(_transfer_function_c0+_transfer_function_c1*log(k));
  };
  double transfer_function_eh98(double kin);
  double transfer_function_eh98_smooth(double kin);
  
  //linear_powspec
  double linear_powspec_exact(double k, double a)
  {
    double gf = growth_function(a);
    double tf = transfer_function(k);
  
    return tf*tf*pow(k,_ns)*gf*gf*_s8*_s8/tophatradnorm_linear_powspec_exact_nonorm(8.0);
  }
  double linear_powspec(double k, double a)
  {
    double gf = growth_function(a);
    double tf;

    if(k < PL_K_MIN)
      {
	tf = transfer_function(k);
	return tf*tf*pow(k,_ns)*_linear_powspec_norm*gf*gf;
      }
    else if(k < PL_K_MAX)
      return exp(gsl_spline_eval(cosmocalc_linear_powspec_spline,log(k),cosmocalc_linear_powspec_acc))*gf*gf;
    else
      return exp(_linear_powspec_c0+_linear_powspec_c1*log(k))*gf*gf;
  };

  //nonlinear powspec
  double gaussiannorm_linear_powspec_exact(double gaussRad);
  double gaussiannorm_linear_powspec(double gaussRad)
  {
    return exp(gsl_spline_eval(cosmocalc_nonlinear_powspec_spline[2],log(gaussRad),cosmocalc_nonlinear_powspec_acc[2]));
  };
  double nonlinear_powspec(double k, double a);
  
  //peak heights
  double sigmaRtophat_exact(double topHatRad, double a)
  {
    return sqrt(_linear_powspec_norm*tophatradnorm_linear_powspec_exact_nonorm(topHatRad))*growth_function(a);
  };
  double sigmaRtophat(double topHatRad, double a)
  {
    return exp(gsl_spline_eval(cosmocalc_R2sigma_spline,log(topHatRad),cosmocalc_R2sigma_acc))*growth_function(a);
  };
  double sigmaMtophat(double m, double a)
  {
    return sigmaRtophat(pow(m/(4.0/3.0*M_PI*RHO_CRIT*_om),1.0/3.0),a);
  };
  double Rsigmatophat(double sigmaR, double a)
  {
    return exp(gsl_spline_eval(cosmocalc_sigma2R_spline,log(sigmaR/growth_function(a)),cosmocalc_sigma2R_acc));
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

/*
  
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

// in utils.c 
double wtime(void);
void gauleg(double x1, double x2, double x[], double w[], int n);

#ifdef FFTLOG
// in fftlog.c 
void compute_discrete_spherical_fft(double *data, int N, double r0, double L, double q, double *result, double *k0);
double get_k0_fftlog(int N, double mu, double q, double r0, double L, double k0guess);
void um_fftlog(int m, double mu, double q, double k0, double r0, double L, double *realpart, double *imagpart);
#endif

// linear_corrfunc.c 
double linear_corrfunc_exact(double r, double a);
double linear_corrfunc(double r, double a);

// in nonlinear_powspec.c 
double nonlinear_powspec_exact(double k, double a);


// nonlinear_corrfunc.c 
double nonlinear_corrfunc_integ_funct(double k, void *p);
double nonlinear_corrfunc_exact(double r, double a);
double nonlinear_corrfunc(double r, double a);

// in mass_bias_functions.c 
double mass_function(double m, double a);
double bias_function(double m, double a);
*/

#endif /* _COSMOCALC_ */
