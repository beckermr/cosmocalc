#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fit.h>

#include "cosmocalc.h"
#include "cosmocalc_assert.h"

void cosmoCalc::init_cosmocalc_transfer_function_table(void)
{
  double *ln_transfer_function_table = (double*)malloc(sizeof(double)*COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH);
  double *lnk_table = (double*)malloc(sizeof(double)*COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH);
  long i;
  double cov00,cov01,cov11,sumsq;
  double dlnk,lnkmin;
  
  cosmocalc_assert(ln_transfer_function_table != NULL,"out of memory for transfer function table!");
  cosmocalc_assert(lnk_table != NULL,"out of memory for transfer function table!");
  
  dlnk = log(TF_K_MAX/TF_K_MIN)/(COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH-1.0);
  lnkmin = log(TF_K_MIN);

#ifdef _OPENMP
  if(_num_threads > 0) omp_set_num_threads(_num_threads);
#endif

#pragma omp parallel for default(none) private(i) shared(dlnk,lnkmin,lnk_table,ln_transfer_function_table)
  for(i=0;i<COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH;++i)
    {
      lnk_table[i] = dlnk*i + lnkmin;
      ln_transfer_function_table[i] = log((this->*_transfer_function_function)(exp(lnk_table[i])));
    }
      
  //init the spline and accelerators
  if(cosmocalc_transfer_function_spline != NULL)
    gsl_spline_free(cosmocalc_transfer_function_spline);
  cosmocalc_transfer_function_spline = gsl_spline_alloc(gsl_interp_cspline,(size_t) (COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH));
  cosmocalc_assert(cosmocalc_transfer_function_spline != NULL,"could not alloc spline for transfer function table!");
  gsl_spline_init(cosmocalc_transfer_function_spline,lnk_table,ln_transfer_function_table,(size_t) (COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH));
  if(cosmocalc_transfer_function_acc != NULL)
    gsl_interp_accel_reset(cosmocalc_transfer_function_acc);
  else
    {
      cosmocalc_transfer_function_acc = gsl_interp_accel_alloc();
      cosmocalc_assert(cosmocalc_transfer_function_acc != NULL,"could not alloc accel for transfer function table!");
    }
  
  gsl_fit_linear(lnk_table+COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH-COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH,(size_t) 1,
		 ln_transfer_function_table+COSMOCALC_TRANSFER_FUNCTION_TABLE_LENGTH-COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH,(size_t) 1,
		 (size_t) COSMOCALC_TRANSFER_FUNCTION_FIT_LENGTH,&_transfer_function_c0,&_transfer_function_c1,&cov00,&cov01,&cov11,&sumsq);
  
  //free mem
  free(ln_transfer_function_table);
  free(lnk_table);
}

double cosmoCalc::transfer_function_eh98(double kin)
{
  //vars
  double k;
  double Tk,Tb,Tc,T0t,T0t11,T0t1bc;
  double omb,om0,omc,h;
  double theta2p7,f,ac,bc,s,q;
  double C,C11,C1bc,st,bb,ab,ksilk,bnode,Gy,y,zeq;
  double a1,a2,b1,b2,keq,Rd,zd,Req,b1d,b2d;

  //get cosmoparms
  omb = _ob;
  om0 = _om;
  omc = _om-_ob;
  h = _h;

  //convert k from hMpc^-1 to Mpc^-1
  k = kin*h;
  
  //-----------
  //input parms 
  //-----------
  theta2p7 = 2.728/2.7;
  
  //eqn 2
  zeq = 2.50*1e4*om0*h*h/(theta2p7*theta2p7*theta2p7*theta2p7);
  
  //eqn 3
  keq = 7.46e-2*om0*h*h/(theta2p7*theta2p7); //Mpc^-{1} (NOT h/Mpc)
  
  //eqn 4
  b1d = 0.313*pow(om0*h*h,-0.419)*(1.0 + 0.607*pow(om0*h*h,0.674));
  b2d = 0.238*pow(om0*h*h,0.223);
  zd = 1291.0*pow(om0*h*h,0.251)/(1.0 + 0.659*pow(om0*h*h,0.828))
    *(1.0 + b1d*pow(omb*h*h,b2d));
  
  //eqn 5
  Rd = 31.5*omb*h*h/(theta2p7*theta2p7*theta2p7*theta2p7)/(zd/1e3);
  Req = 31.5*omb*h*h/(theta2p7*theta2p7*theta2p7*theta2p7)/(zeq/1e3);
  
  //eqn 6
  s = 2.0/3.0/keq*sqrt(6.0/Req)*log((sqrt(1.0 + Rd) + sqrt(Rd + Req))/(1.0 + sqrt(Req)));
  
  //eqn 7
  ksilk = 1.6*pow(omb*h*h,0.52)*pow(om0*h*h,0.73)*(1.0 + pow(10.4*om0*h*h,-0.95));
  
  //eqn 10
  q = k/13.41/keq;
  
  //eqn 11
  a1 = pow(46.9*om0*h*h,0.670)*(1.0 + pow(32.1*om0*h*h,-0.532));
  a2 = pow(12.0*om0*h*h,0.424)*(1.0 + pow(45.0*om0*h*h,-0.582));
  ac = pow(a1,-1.0*omb/om0)*pow(a2,-1.0*(omb/om0)*(omb/om0)*(omb/om0));
  
  //eqn 12
  b1 = 0.944/(1.0 + pow(458.0*om0*h*h,-0.708));
  b2 = pow(0.395*om0*h*h,-0.0266);
  bc = 1.0/(1.0 + b1*(pow(omc/om0,b2) - 1.0));
  
  //eqn 15
  y = (1.0 + zeq)/(1.0 + zd);
  Gy = y*(-6.0*sqrt(1.0 + y) + (2.0 + 3.0*y)*log((sqrt(1.0 + y) + 1.0)/(sqrt(1.0 + y) - 1.0)));
  
  //eqn 14
  ab = 2.07*keq*s*pow(1.0 + Rd,-3.0/4.0)*Gy;
  
  //----------------------------------
  // Get CDM part of transfer function
  //----------------------------------
  
  //eqn 18
  f = 1.0/(1.0 + (k*s/5.4)*(k*s/5.4)*(k*s/5.4)*(k*s/5.4));
  
  //eqn 20
  C = 14.2/ac + 386.0/(1.0 + 69.9*pow(q,1.08));           
  
  //eqn 19
  T0t = log(M_E + 1.8*bc*q)/(log(M_E + 1.8*bc*q) + C*q*q);    
  
  //eqn 17
  C1bc = 14.2 + 386.0/(1.0 + 69.9*pow(q,1.08));           
  T0t1bc = log(M_E + 1.8*bc*q)/(log(M_E + 1.8*bc*q) + C1bc*q*q);    
  Tc = f*T0t1bc + (1.0 - f)*T0t;
  
  //-------------------------------------
  // Get baryon part of transfer function
  //-------------------------------------
  
  //eqn 24
  bb = 0.5 + omb/om0 + (3.0 - 2.0*omb/om0)*sqrt((17.2*om0*h*h)*(17.2*om0*h*h) + 1.0);
  
  //eqn 23
  bnode = 8.42*pow(om0*h*h,0.435);
  
  //eqn 22
  st = s/pow(1.0 + (bnode/k/s)*(bnode/k/s)*(bnode/k/s),1.0/3.0);
    
  //eqn 21
  C11 = 14.2 + 386.0/(1.0 + 69.9*pow(q,1.08));
  T0t11 = log(M_E + 1.8*q)/(log(M_E + 1.8*q) + C11*q*q);    
  Tb = (T0t11/(1.0 + (k*s/5.2)*(k*s/5.2)) 
	+ ab/(1.0 + (bb/k/s)*(bb/k/s)*(bb/k/s))/exp(pow(k/ksilk,1.4)))*sin(k*st)/(k*st);

  //------------------------
  // total transfer function
  //------------------------
  Tk = omb/om0*Tb + omc/om0*Tc;
  
  return Tk;
}	

double cosmoCalc::transfer_function_eh98_smooth(double kin)
{
  //vars
  double k,Tk;
  double omb,om0,omc,h;
  double theta2p7,s,q;
  double Gamma,alphaGamma,L0,C0;
  
  //get cosmoparms
  omb = _ob;
  om0 = _om;
  omc = _om - _ob;
  h = _h;
  
  //convert k from hMpc^-1 to Mpc^-1
  k = kin*h;
  
  //-----------
  //input parms 
  //-----------
  theta2p7 = 2.728/2.7;
  
  //eqn 26
  s = 44.5*log(9.83/om0/h/h)/sqrt(1.0 + 10.0*pow(omb*h*h,0.75));
  
  //eqn 31
  alphaGamma = 1.0 - 0.328*log(431.0*om0*h*h)*omb/om0 + 0.38*log(22.3*om0*h*h)*(omb/om0)*(omb/om0);
  
  //eqn 30
  Gamma = om0*h*(alphaGamma + (1.0 - alphaGamma)/(1.0 + pow(0.43*k*s,4.0)));
  
  //eqn 28
  q = kin*theta2p7*theta2p7/Gamma;
  
  //eqns 29
  C0 = 14.2 + 731.0/(1.0 + 62.5*q);
  L0 = log(2.0*exp(1.0) + 1.8*q);
  Tk = L0/(L0 + C0*q*q);
  
  return Tk;
}	
