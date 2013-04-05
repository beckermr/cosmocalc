#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>

#include "cosmocalc.h"

static double deltaTm[9] = {200.0, 300.0, 400.0, 600.0, 800.0, 1200.0, 1600.0, 2400.0, 3200.0};
static double ATm[9] = {0.186, 0.2, 0.212, 0.218, 0.248, 0.255, 0.260, 0.260, 0.260};
static double aTm[9] = {1.47, 1.52, 1.56, 1.61, 1.87, 2.13, 2.30, 2.53, 2.66};
static double bTm[9] = {2.57, 2.25, 2.05, 1.87, 1.59, 1.51, 1.46, 1.44, 1.41};
static double cTm[9] = {1.19, 1.27, 1.34, 1.45, 1.58, 1.80, 1.97, 2.24, 2.44};

static double tinker2010_bias(double nu, double delta);
static double tinker2010_mf_norm_integ(double nu, void *p);
  
double mass_function(double m, double a)
{
  return tinker2008_mass_function(m,a,cosmoData.delta);
}

double tinker2008_mass_function(double m, double a, double delta)
{
  double A,am,b,c,alpha;
  
  static int init = 1;
  static gsl_spline *ATm_spline,*aTm_spline,*bTm_spline,*cTm_spline;
  static gsl_interp_accel *ATm_acc,*aTm_acc,*bTm_acc,*cTm_acc;
  
  if(init)
    {
#define TINKERCODEPARMS
#ifdef TINKERCODEPARMS
      //first parameter
      ATm[0] = 1.858659e-01 ;
      ATm[1] = 1.995973e-01 ;
      ATm[2] = 2.115659e-01 ;
      ATm[3] = 2.184113e-01 ;
      ATm[4] = 2.480968e-01 ;
      ATm[5] = 2.546053e-01 ;
      ATm[6] = 2.600000e-01 ;
      ATm[7] = 2.600000e-01 ;
      ATm[8] = 2.600000e-01 ;
      
      //second parameter
      aTm[0] = 1.466904e+00 ;
      aTm[1] = 1.521782e+00 ;
      aTm[2] = 1.559186e+00 ;
      aTm[3] = 1.614585e+00 ;
      aTm[4] = 1.869936e+00 ;
      aTm[5] = 2.128056e+00 ;
      aTm[6] = 2.301275e+00 ;
      aTm[7] = 2.529241e+00 ;
      aTm[8] = 2.661983e+00 ;
      
      //third parameter
      bTm[0] = 2.571104e+00 ;
      bTm[1] = 2.254217e+00 ;
      bTm[2] = 2.048674e+00 ;
      bTm[3] = 1.869559e+00 ;
      bTm[4] = 1.588649e+00 ;
      bTm[5] = 1.507134e+00 ;
      bTm[6] = 1.464374e+00 ;
      bTm[7] = 1.436827e+00 ;
      bTm[8] = 1.405210e+00 ;

      //fourth parameter
      cTm[0] = 1.193958e+00;
      cTm[1] = 1.270316e+00;
      cTm[2] = 1.335191e+00;
      cTm[3] = 1.446266e+00;
      cTm[4] = 1.581345e+00;
      cTm[5] = 1.795050e+00;
      cTm[6] = 1.965613e+00;
      cTm[7] = 2.237466e+00;
      cTm[8] = 2.439729e+00;
#endif
      
      ATm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      aTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      bTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      cTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      
      gsl_spline_init(ATm_spline,deltaTm,ATm,9);
      gsl_spline_init(aTm_spline,deltaTm,aTm,9);
      gsl_spline_init(bTm_spline,deltaTm,bTm,9);
      gsl_spline_init(cTm_spline,deltaTm,cTm,9);
      
      ATm_acc = gsl_interp_accel_alloc();
      aTm_acc = gsl_interp_accel_alloc();
      bTm_acc = gsl_interp_accel_alloc();
      cTm_acc = gsl_interp_accel_alloc();
      
      init = 0;
    }
  
  A = gsl_spline_eval(ATm_spline,delta,ATm_acc);
#ifdef TINKERCODEPARMS
  if(delta > 1600.0)
    A = 0.26;
#endif
  am = gsl_spline_eval(aTm_spline,delta,aTm_acc);
  b = gsl_spline_eval(bTm_spline,delta,bTm_acc);
  c = gsl_spline_eval(cTm_spline,delta,cTm_acc);
  
  double onepz = 1.0/a;
#ifdef TINKERCODEPARMS
  if(onepz > 4.0)
    onepz = 4.0;
#endif
  
  A *= pow(onepz,-0.14);
  am *= pow(onepz,-0.06);
  alpha = pow(10.0,-1.0*pow(0.75/log10(delta/75.0),1.2));
  b *= pow(onepz,-1.0*alpha);
  
  //fprintf(stderr,"A = %f, a = %f, b = %f, c = %f\n",A,am,b,c);
  
  double sigma = sigmaMtophat(m,a);
  double fsigma = A*(pow(sigma/b,-1.0*am) + 1.0)*exp(-1.0*c/sigma/sigma);
  double dm = 1e-6*m;
  double dlnsiginvdm = log(sigmaMtophat(m-dm/2.0,a)/sigmaMtophat(m+dm/2.0,a))/dm;
  
  return fsigma*RHO_CRIT*cosmoData.OmegaM/m*dlnsiginvdm;
}

static double alphaTm[9] = {0.368, 0.363, 0.385, 0.389, 0.393, 0.365, 0.379, 0.355, 0.327};
static double betaTm[9] = {0.589, 0.585, 0.544, 0.543, 0.564, 0.623, 0.637, 0.673, 0.702};
static double gammaTm[9] = {0.864, 0.922, 0.987, 1.09, 1.20, 1.34, 1.50, 1.68, 1.81};
static double phiTm[9] = {-0.729, -0.789, -0.910, -1.05, -1.20, -1.26, -1.45, -1.50, -1.49};
static double etaTm[9] = {-0.243, -0.261, -0.261, -0.273, -0.278, -0.301, -0.301, -0.319, -0.336};

double tinker2010_mass_function(double m, double a, double delta)
{
  double z,gsigma,nu,dlnsiginvdm;
  double dm = 1e-6*m;
  double alpha,beta,gamma,phi,eta;
  
  static int init = 1;
  static gsl_spline *alphaTm_spline,*betaTm_spline,*gammaTm_spline,*phiTm_spline,*etaTm_spline;
  static gsl_interp_accel *alphaTm_acc,*betaTm_acc,*gammaTm_acc,*phiTm_acc,*etaTm_acc;
  
  if(init)
    {
      alphaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      betaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      gammaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      phiTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      etaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      
      gsl_spline_init(alphaTm_spline,deltaTm,alphaTm,9);
      gsl_spline_init(betaTm_spline,deltaTm,betaTm,9);
      gsl_spline_init(gammaTm_spline,deltaTm,gammaTm,9);
      gsl_spline_init(phiTm_spline,deltaTm,phiTm,9);
      gsl_spline_init(etaTm_spline,deltaTm,etaTm,9);
      
      alphaTm_acc = gsl_interp_accel_alloc();
      betaTm_acc = gsl_interp_accel_alloc();
      gammaTm_acc = gsl_interp_accel_alloc();
      phiTm_acc = gsl_interp_accel_alloc();
      etaTm_acc = gsl_interp_accel_alloc();
      
      init = 0;
    }
  
  alpha = gsl_spline_eval(alphaTm_spline,delta,alphaTm_acc);
  beta = gsl_spline_eval(betaTm_spline,delta,betaTm_acc);
  gamma = gsl_spline_eval(gammaTm_spline,delta,gammaTm_acc);
  phi = gsl_spline_eval(phiTm_spline,delta,phiTm_acc);
  eta = gsl_spline_eval(etaTm_spline,delta,etaTm_acc);
  
  nu = 1.686/sigmaMtophat(m,a);
  z = 1.0/a - 1.0;
  if(z > 3.0)
    z = 3.0;
  
  beta = beta*pow(1.0+z,0.20);
  phi = phi*pow(1.0+z,-0.08);
  eta = eta*pow(1.0+z,0.27);
  gamma = gamma*pow(1.0+z,-0.01);
  
  gsigma = alpha*(1.0 + pow(beta*nu,-2.0*phi))*pow(nu,2.0*eta)*exp(-1.0*gamma*nu*nu/2.0)*nu;
  dlnsiginvdm = log(sigmaMtophat(m-dm/2.0,a)/sigmaMtophat(m+dm/2.0,a))/dm;
  
  return gsigma*RHO_CRIT*cosmoData.OmegaM/m*dlnsiginvdm;
}

/*double tinker2010_mass_function(double m, double a, double delta)
{
  double z,gsigma,nu,dlnsiginvdm;
  double dm = 1e-6*m;
  double alpha,beta,gamma,phi,eta;
  double params[5];
  gsl_function f;
  size_t limit = 1000;
  double epsabs,epsrel;
  double abserr;
  double numin,numax,dnu;
  int i,Nnu;
  
  static int init = 1;
  static gsl_spline *betaTm_spline,*gammaTm_spline,*phiTm_spline,*etaTm_spline;
  static gsl_interp_accel *betaTm_acc,*gammaTm_acc,*phiTm_acc,*etaTm_acc;
  static gsl_integration_glfixed_table *gltab;
  static gsl_integration_workspace *w;
  
  if(init)
    {
      betaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      gammaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      phiTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      etaTm_spline = gsl_spline_alloc(gsl_interp_cspline,9);
      
      gsl_spline_init(betaTm_spline,deltaTm,betaTm,9);
      gsl_spline_init(gammaTm_spline,deltaTm,gammaTm,9);
      gsl_spline_init(phiTm_spline,deltaTm,phiTm,9);
      gsl_spline_init(etaTm_spline,deltaTm,etaTm,9);
      
      betaTm_acc = gsl_interp_accel_alloc();
      gammaTm_acc = gsl_interp_accel_alloc();
      phiTm_acc = gsl_interp_accel_alloc();
      etaTm_acc = gsl_interp_accel_alloc();
      
      gltab = gsl_integration_glfixed_table_alloc((size_t) 10);
      w = gsl_integration_workspace_alloc(limit);
      
      init = 0;
    }
  
  beta = gsl_spline_eval(betaTm_spline,delta,betaTm_acc);
  gamma = gsl_spline_eval(gammaTm_spline,delta,gammaTm_acc);
  phi = gsl_spline_eval(phiTm_spline,delta,phiTm_acc);
  eta = gsl_spline_eval(etaTm_spline,delta,etaTm_acc);
  
  nu = 1.686/sigmaMtophat(m,a);
  z = 1.0/a - 1.0;
  if(z > 3.0)
    z = 3.0;
  
  beta = beta*pow(1.0+z,0.20);
  phi = phi*pow(1.0+z,-0.08);
  eta = eta*pow(1.0+z,0.27);
  gamma = gamma*pow(1.0+z,-0.01);
  
  //determine alpha by normalization condition
  params[0] = beta;
  params[1] = phi;
  params[2] = eta;
  params[3] = gamma;
  params[4] = delta;
  f.function = &tinker2010_mf_norm_integ;
  f.params = params;
  numin = 0.0;
  numax = 20.0;
  epsabs = 0.0;
  epsrel = 0.05;
  
  if(1)
    gsl_integration_qag(&f,numin,numax,epsabs,epsrel,limit,GSL_INTEG_GAUSS61,w,&alpha,&abserr);
  else if(0)
    alpha = gsl_integration_glfixed(&f,numin,numax,gltab);
  else if(0)
    {
      alpha = 0.0;
      Nnu = 100;
      
      dnu = (numax-numin)/Nnu;
      for(i=0;i<Nnu;++i)
	{
	  alpha += tinker2010_mf_norm_integ(numin+dnu*i,f.params)*dnu;
	  //fprintf(stderr,"nu = %lf, integ = %lf\n",numin+dnu*i,tinker2010_mf_norm_integ(numin+dnu*i,f.params)*dnu);
	}
    }
  else if(0)
    gsl_integration_qng(&f,numin,numax,epsabs,epsrel,&alpha,&abserr,&limit);
  
  alpha = 1.0/alpha;
  //fprintf(stderr,"alpha = %lf\n",alpha);
  
  gsigma = alpha*(1.0 + pow(beta*nu,-2.0*phi))*pow(nu,2.0*eta)*exp(-1.0*gamma*nu*nu/2.0)*nu;
  dlnsiginvdm = log(sigmaMtophat(m-dm/2.0,a)/sigmaMtophat(m+dm/2.0,a))/dm;
  
  return gsigma*RHO_CRIT*cosmoData.OmegaM/m*dlnsiginvdm;
}

static double tinker2010_mf_norm_integ(double nu, void *p)
{
  double *params = (double*) p;
  double fnu,bnu;
  
  if(nu == 0.0)
    return 0.0;
  else
    {
      fnu = (1.0 + pow(params[0]*nu,-2.0*params[1]))*pow(nu,2.0*params[2])*exp(-1.0*params[3]*nu*nu/2.0);
      bnu = tinker2010_bias(nu,params[4]);
      return bnu*fnu;
    }
}
*/

double bias_function(double m, double a)
{
  double nu;
  nu =  1.686/sigmaMtophat(m,a);
  
  return tinker2010_bias(nu,cosmoData.delta);
}

static double tinker2010_bias(double nu, double delta)
{
  double y = log10(delta);
  double A = 1.0 + 0.24*y*exp(-1.0*pow(4.0/y,4.0));
  double a = 0.44*y - 0.88;
  double B = 0.183;
  double b = 1.5;
  double C = 0.019 + 0.107*y + 0.19*exp(-1.0*pow(4.0/y,4.0));
  double c = 2.4;
    
  return 1.0 - A*pow(nu,a)/(pow(nu,a) + pow(1.686,a)) + B*pow(nu,b) + C*pow(nu,c);
}
