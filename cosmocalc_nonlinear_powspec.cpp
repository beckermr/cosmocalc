#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "cosmocalc.h"
#include "cosmocalc_assert.h"

struct nonlinear_powspec_data {
  double param;
  cosmoCalc *cd;
};

static double gaussiannorm_linear_powspec_exact_lnk_integ_funct(double lnk, void *p)
{
  struct nonlinear_powspec_data *dat = (struct nonlinear_powspec_data*)p;
  double gaussRad = dat->param;
  double k = exp(lnk);
  return dat->cd->linear_powspec(k,1.0)*k*k*k/2.0/M_PI/M_PI*exp(-1.0*k*k*gaussRad*gaussRad);
}

void cosmoCalc::init_cosmocalc_nonlinear_powspec_table(void)
{
  long i;
  double *xtab = (double*)malloc(sizeof(double)*COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH);
  double *ytab = (double*)malloc(sizeof(double)*COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH);
  
  cosmocalc_assert(xtab != NULL,"out of memory for nonlinear powspec table!");
  cosmocalc_assert(ytab != NULL,"out of memory for nonlinear powspec table!");
  
  for(i=0;i<3;++i)
    {
      if(cosmocalc_nonlinear_powspec_spline[i] != NULL)
	gsl_spline_free(cosmocalc_nonlinear_powspec_spline[i]);
      cosmocalc_nonlinear_powspec_spline[i] = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
      cosmocalc_assert(cosmocalc_nonlinear_powspec_spline[i] != NULL,"could not alloc spline %ld for nonlinear powspec!",i);
            
      if(cosmocalc_nonlinear_powspec_acc[i] != NULL)
	gsl_interp_accel_reset(cosmocalc_nonlinear_powspec_acc[i]);
      else
	{
	  cosmocalc_nonlinear_powspec_acc[i] = gsl_interp_accel_alloc();
	  cosmocalc_assert(cosmocalc_nonlinear_powspec_acc[i] != NULL,"could not alloc accel %ld for nonlinear powspec table!",i);
	}
    }
  
  double dlnr = log(PNL_RGAUSS_MAX/PNL_RGAUSS_MIN)/(COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH-1.0);
  double lnrmin = log(PNL_RGAUSS_MIN);
  double lnr;
  
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
  struct nonlinear_powspec_data dat;
  double gaussRad;
  
  dat.cd = this;
  F.params = &dat;
  F.function = &gaussiannorm_linear_powspec_exact_lnk_integ_funct;
  
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 0.0
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  cosmocalc_assert(workspace != NULL,"could not alloc workspace for GSL integration in gauss norm for nonlinear powspec!");
  
  for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
    {
      lnr = dlnr*i + lnrmin;
      xtab[i] = lnr;
      gaussRad = exp(lnr);
      dat.param = gaussRad;
      
      gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
      gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
            
      ytab[i] = log(I0 + I1);
    }

  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  gsl_spline_init(cosmocalc_nonlinear_powspec_spline[2],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
    
  for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
    {
      xtab[i] = i*(PNL_A_MAX-PNL_A_MIN)/(COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH-1.0) + PNL_A_MIN;
      ytab[i] = get_nonlinear_gaussnorm_scale(xtab[i]);
    }
  gsl_spline_init(cosmocalc_nonlinear_powspec_spline[0],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
  
  for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
    xtab[i] = ytab[i];
  for(i=0;i<COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH;++i)
    ytab[i] = gaussiannorm_linear_powspec(xtab[i]);
  gsl_spline_init(cosmocalc_nonlinear_powspec_spline[1],xtab,ytab,(size_t) (COSMOCALC_NONLINEAR_POWSPEC_TABLE_LENGTH));
  
  free(xtab);
  free(ytab);
}

double cosmoCalc::nonlinear_powspec(double k, double a)
{
  double Rsigma,C,neff,ksigma,sigma2;
  double an,bn,cn,alphan,gamman,betan,mun,nun;
  double f1,f2,f3;
  double DeltakNL,dsigma2dR,d2sigma2d2R,PkNL,PkL;
  double y,DeltakL,fy,DeltakQ,DeltakHprime,DeltakH;
  
  Rsigma = gsl_spline_eval(cosmocalc_nonlinear_powspec_spline[0],
			   a,cosmocalc_nonlinear_powspec_acc[0]); //get_nonlinear_gaussnorm_scale(a);
  sigma2 = gsl_spline_eval(cosmocalc_nonlinear_powspec_spline[1],
			   Rsigma,cosmocalc_nonlinear_powspec_acc[1]); //gaussiannorm_linear_powspec_exact(Rsigma);
  dsigma2dR = gsl_spline_eval_deriv(cosmocalc_nonlinear_powspec_spline[1],
				    Rsigma,cosmocalc_nonlinear_powspec_acc[1]); //onederiv_gaussiannorm_linear_powspec_exact(Rsigma);
  d2sigma2d2R = gsl_spline_eval_deriv2(cosmocalc_nonlinear_powspec_spline[1],
				       Rsigma,cosmocalc_nonlinear_powspec_acc[1]); //twoderiv_gaussiannorm_linear_powspec_exact(Rsigma);
  
  ksigma = 1.0/Rsigma;
  neff = -1.0*Rsigma/sigma2*dsigma2dR - 3.0;
  C = -1.0*(d2sigma2d2R*Rsigma*Rsigma/sigma2 + dsigma2dR*Rsigma/sigma2 - dsigma2dR*dsigma2dR*Rsigma*Rsigma/sigma2/sigma2);

#ifdef SMITH03
  an = pow(10.0,1.4861 + 1.8369*neff + 1.6762*neff*neff + 0.7940*neff*neff*neff + 0.1670*neff*neff*neff*neff - 0.6206*C);
  bn = pow(10.0,0.9463 + 0.9466*neff + 0.3084*neff*neff - 0.9400*C);
  cn = pow(10.0,-0.2807 + 0.6669*neff + 0.3214*neff*neff - 0.0793*C);
  gamman = 0.8649 + 0.2989*neff + 0.1631*C;
  alphan = 1.3884 + 0.3700*neff - 0.1452*neff*neff;
  betan = 0.8291 + 0.9854*neff + 0.3401*neff*neff;
  mun = pow(10.0,-3.5442 + 0.1908*neff);
  nun = pow(10.0,0.9589 + 1.2857*neff);
  
  f1 = pow(cosmoData.OmegaM,-0.0307);
  f2 = pow(cosmoData.OmegaM,-0.0585);
  f3 = pow(cosmoData.OmegaM,0.0743);
#else
  double w0,wa,ha,weffa,omegaMz,omegaDEwz;
  
  w0 = _w0;
  wa = _wa;
  ha = hubble_noscale(a);
  if(a != 1.0)
    weffa = w0 + wa - wa*(a - 1.0)/log(a);
  else
    weffa = w0;
  omegaMz = _om/a/a/a/ha/ha;
  omegaDEwz = _ol/ha/ha/pow(a,3.0*(1.0 + weffa));
      
  an = pow(10.0,1.5222 + 2.8553*neff + 2.3706*neff*neff + 0.9903*neff*neff*neff + 0.2250*neff*neff*neff*neff - 0.6038*C + 0.1749*omegaDEwz*(1.0 + weffa));
  bn = pow(10.0,-0.5642 + 0.5864*neff + 0.5716*neff*neff - 1.5474*C + 0.2279*omegaDEwz*(1.0 + weffa));
  cn = pow(10.0,0.3698 + 2.0404*neff + 0.8161*neff*neff + 0.5869*C);
  gamman = 0.1971 - 0.0843*neff + 0.8460*C;
  alphan = fabs(6.0835 + 1.3373*neff - 0.1959*neff*neff - 5.5274*C);
  betan = 2.0379 - 0.7354*neff + 0.3157*neff*neff + 1.2490*neff*neff*neff + 0.3980*neff*neff*neff*neff - 0.1682*C;
  mun = 0.0;
  nun = pow(10.0,5.2105 + 3.6902*neff);
      
  f1 = pow(omegaMz,-0.0307);
  f2 = pow(omegaMz,-0.0585);
  f3 = pow(omegaMz,0.0743);
#endif
  
  PkL = linear_powspec(k,a);
  y = k/ksigma;
  fy = y/4.0 + y*y/8.0;
  DeltakL = PkL*k*k*k/2.0/M_PI/M_PI;
  
  DeltakQ = DeltakL*pow(1.0 + DeltakL,betan)/(1.0 + alphan*DeltakL)*exp(-1.0*fy);
  
  DeltakHprime = an*pow(y,3.0*f1)/(1.0 + bn*pow(y,f2) + pow(cn*f3*y,3.0 - gamman));
  DeltakH = DeltakHprime/(1.0 + mun/y + nun/y/y);
  
  DeltakNL = DeltakQ + DeltakH;
  PkNL = DeltakNL/(k*k*k/2.0/M_PI/M_PI);
  
  return PkNL;
}

double cosmoCalc::gaussiannorm_linear_powspec_exact(double gaussRad)
{
  double I0,I1;
  double abserr;
  gsl_integration_workspace *workspace;
  gsl_function F;
  struct nonlinear_powspec_data dat;  
  
  dat.param = gaussRad;
  dat.cd = this;
  
#define WORKSPACE_NUM 10000000
#define ABSERR 1e-6
#define RELERR 1e-2
  workspace = gsl_integration_workspace_alloc((size_t) WORKSPACE_NUM);
  cosmocalc_assert(workspace != NULL,"could not alloc workspace for GSL integration in gauss norm for nonlinear powspec!");
  
  F.params = &dat;
  F.function = &gaussiannorm_linear_powspec_exact_lnk_integ_funct;
  gsl_integration_qags(&F,log(1e-4),log(2.0*M_PI/gaussRad),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I0,&abserr);
  gsl_integration_qags(&F,log(2.0*M_PI/gaussRad),log(1e3),ABSERR,RELERR,(size_t) WORKSPACE_NUM,workspace,&I1,&abserr);
    
  gsl_integration_workspace_free(workspace);
#undef ABSERR
#undef RELERR
#undef WORKSPACE_NUM
  
  return I0 + I1;
}

static double nonlinear_gaussnorm_scale_funct(double gaussR, void *p)
{
  struct nonlinear_powspec_data *dat = (struct nonlinear_powspec_data*)p;
  double gf = dat->param;
  
  return dat->cd->gaussiannorm_linear_powspec(gaussR)*gf*gf-1.0;
}

inline double cosmoCalc::get_nonlinear_gaussnorm_scale(double a)
{
  double gf = growth_function(a);
  double Rsigma,Rlow=PNL_RGAUSS_MIN,Rhigh=PNL_RGAUSS_MAX;
  int itr,maxItr=1000,status;
  struct nonlinear_powspec_data dat;  
  dat.param = gf;
  dat.cd = this;
  
#define ABSERR 1e-6
#define RELERR 1e-6
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  gsl_function F;
       
  F.function = &nonlinear_gaussnorm_scale_funct;
  F.params = &dat;
  
  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  cosmocalc_assert(s != NULL,"could not alloc GSL root solver for nonlinear powspec non-lin k computation!");
  
  gsl_root_fsolver_set(s,&F,Rlow,Rhigh);
  itr = 0;
  
  do
    {
      itr++;
      status = gsl_root_fsolver_iterate(s);
      Rsigma = gsl_root_fsolver_root(s);
      Rlow = gsl_root_fsolver_x_lower(s);
      Rhigh = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(Rlow,Rhigh,ABSERR,RELERR);
    }
  while(status == GSL_CONTINUE && itr < maxItr);
  
#undef ABSERR
#undef RELERR

  gsl_root_fsolver_free(s);
  
  return Rsigma;
}
