#include <cstdlib>
#include <cstdio>

#include "flrw_distances.h"
#include "w0wacosmo.h"
#include "eh98_transfunct.h"
#include "linpowspec.h"
#include "halofit.h"
#include "peakheight.h"
#include "tinkerhalostats.h"

#ifndef _OPENMP
#include <sys/time.h>
static double wtime(void)
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ((double) (tp.tv_sec)) + ((double) (tp.tv_usec))/1e6;
}
#else
#include <omp.h>
#define wtime omp_get_wtime
#endif

void printcosmo(double a, class CosmoData& cd, class Hubble& h, class Distances& d, class GrowthFunction& gf,
		class TransferFunction& tf, class TransferFunction& tfs,
		class PowerSpectrum& lp, class PowerSpectrum& pknl, class PeakHeight& ph)
{
  fprintf(stderr,"weff(%f) = %f\n",a,cd(a));
  fprintf(stderr,"H(%f)/H100 = %f\n",a,h(a));
  fprintf(stderr,"comvdist(%f) = %f\n",a,d.comvdist(a));
  fprintf(stderr,"exact comvdist(%f) = %f\n",a,d.comvdist_exact(a,h));
  fprintf(stderr,"angdist(%f) = %f\n",a,d.angdist(a));
  fprintf(stderr,"lumdist(%f) = %f\n",a,d.lumdist(a));
  fprintf(stderr,"growth function(%f) = %f\n",a,gf(1.0,a));
  fprintf(stderr,"sigma8 = %f\n",ph.sigmaRtophat(8.0,1.0));
  
  double kmin = 1e-6;
  double kmax = 1e4;
  long Nk = 2000;
  double k,dlnk = log(kmax/kmin)/Nk;
  long i;
  
  FILE *fp;
  fp = fopen("test.dat","w");
  for(i=0;i<Nk;++i)
    {
      k = exp(dlnk*i)*kmin;
      fprintf(fp,"%e %e %e %e %e\n",k,tf(k),tfs(k),lp(k,a),pknl(k,a));
    }
  fclose(fp);
}

void initcosmo(class w0wa_CosmoData& cd, class w0wa_Hubble& h, class FLRWDistances& d, class w0wa_GrowthFunction& gf,
	       class EH98_TransferFunction& tf, class EH98Smooth_TransferFunction& tfs,
	       class LinearPowerSpectrum& lp, class HaloFitPowerSpectrum& pknl, class PeakHeight& ph)
{
  double t;
  double t0 = -wtime();
  t = -wtime();
  h.init(cd);
  t += wtime();
  fprintf(stderr,"\nh init took %g seconds.\n",t);
  t = -wtime();
  d.init(cd.ok,h);
  t += wtime();
  fprintf(stderr,"d init took %g seconds.\n",t);
  t = -wtime();
  gf.init(cd,h);
  t += wtime();
  fprintf(stderr,"gf init took %g seconds.\n",t);
  t = -wtime();
  tf.init(cd.om,cd.ob,cd.h);
  t += wtime();
  fprintf(stderr,"tf init took %g seconds.\n",t);
  t = -wtime();
  tfs.init(cd.om,cd.ob,cd.h);
  t += wtime();
  fprintf(stderr,"tfs init took %g seconds.\n",t);
  t = -wtime();
  lp.init(cd.s8,cd.ns,gf,tf);
  t += wtime();
  fprintf(stderr,"pkl init took %g seconds.\n",t);
  
  t = -wtime();
  pknl.init(cd.om,cd.ol,cd,lp,h,gf);
  t += wtime();
  fprintf(stderr,"pknl init took %g seconds.\n",t);
  
  t = -wtime();
  ph.init(cd.om,gf,lp);
  t += wtime();
  fprintf(stderr,"ph init took %g seconds.\n",t);
  
  t0 += wtime();
  fprintf(stderr,"total init took %g seconds.\n\n",t0);
}

int main(int argc, char **argv)
{
  w0wa_CosmoData cd;
  w0wa_Hubble h;
  FLRWDistances d;
  w0wa_GrowthFunction gf;
  EH98_TransferFunction tf;
  EH98Smooth_TransferFunction tfs;
  LinearPowerSpectrum lp;
  HaloFitPowerSpectrum pknl;
  PeakHeight ph;
  Tinker08MassFunction mf;
  double a = atof(argv[1]);
  
  fprintf(stderr,"Testing the cosmology routines...\n");
    
  //set params
  cd.om = 0.25;
  cd.ol = 0.75;
  cd.ob = 0.045;
  cd.onu = 0.0;
  cd.ok = 0.0;
  cd.h = 0.7;
  cd.s8 = 0.8;
  cd.ns = 1.0;
  cd.w0 = -1.0;
  cd.wa = 0.0;
  
  initcosmo(cd,h,d,gf,tf,tfs,lp,pknl,ph);
  printcosmo(a,cd,h,d,gf,tf,tfs,lp,pknl,ph);
  fprintf(stderr,"exact growth function(%f) = %f, norm = %f\n",a,gf.growth_function_exact(1.0,a,h),gf.growth_function_norm());
  
  mf.init(cd.om,200.0,ph);
  fprintf(stderr,"n(%g,%f) = %g\n",1e13,1.0,mf(1e13,1.0));
  
  initcosmo(cd,h,d,gf,tf,tfs,lp,pknl,ph);
  printcosmo(a,cd,h,d,gf,tf,tfs,lp,pknl,ph);
  fprintf(stderr,"exact growth function(%f) = %f, norm = %f\n",a,gf.growth_function_exact(1.0,a,h),gf.growth_function_norm());
  
  //set params
  cd.om = 0.25;
  cd.ol = 0.75;
  cd.ob = 0.045;
  cd.onu = 0.0;
  cd.ok = 0.0;
  cd.h = 0.7;
  cd.s8 = 0.85;
  cd.ns = 1.0;
  cd.w0 = -1.0;
  cd.wa = 0.0;
  
  initcosmo(cd,h,d,gf,tf,tfs,lp,pknl,ph);
  printcosmo(a,cd,h,d,gf,tf,tfs,lp,pknl,ph);
  fprintf(stderr,"exact growth function(%f) = %f, norm = %f\n",a,gf.growth_function_exact(1.0,a,h),gf.growth_function_norm());
    
  cosmocalc_assert(cd == h.cosmology(),"cosmology in h object is not the same as that use to init!");
  
  /*
  cd.init_cosmology(0.25,0.75,0.145,0.0,0.7,0.8,1.0,-1.0,0.0,COSMOCALC_TRANS_FUNC_EH98);
  cd.init_all();
  if(argc > 3)
    {
      fprintf(stderr,"setting library to use %d threads!\n",atoi(argv[3]));
      cd.num_threads(atoi(argv[3]));
    }
    
  

  t = -wtime();
  cd.init_cosmology(0.20,0.80,0.045,0.0,0.7,0.8,1.0,-1.0,0.0,COSMOCALC_TRANS_FUNC_EH98);
  cd.init_all();
  t += wtime();
  fprintf(stderr,"second init took %g seconds.\n\n",t);
  
  fprintf(stderr,"omegam = %f\n",cd.omegam());
  fprintf(stderr,"omegal = %f\n",cd.omegal());
  fprintf(stderr,"omegab = %f\n",cd.omegab());
  fprintf(stderr,"omeganu = %f\n",cd.omeganu());
  fprintf(stderr,"omegak = %f\n",cd.omegak());
  fprintf(stderr,"h = %f\n",cd.h());
  fprintf(stderr,"sigma8 = %f\n",cd.sigma8());
  fprintf(stderr,"ns = %f\n",cd.ns());
  fprintf(stderr,"w0,wa = %f|%f\n",cd.w0(),cd.wa());
  
  fprintf(stderr,"\n");
  
  double a = atof(argv[1]);
  fprintf(stderr,"comvdist(%f) = %f\n",a,cd.comvdist(a));
  //fprintf(stderr,"comvdist(%f) = %f\n",a,cd.comvdist_exact(a));
  fprintf(stderr,"angdist(%f) = %f\n",a,cd.angdist(a));
  fprintf(stderr,"lumdist(%f) = %f\n",a,cd.lumdist(a));
  
  //fprintf(stderr,"growth function(%f) = %f|%f\n",a,cd.growth_function(a),cd.growth_function_exact(a));
  fprintf(stderr,"growth function(%f) = %f, norm = %f\n",a,cd.growth_function(a),cd.growth_function_norm());
  
  fprintf(stderr,"\n");
  
  double k = atof(argv[2]);
  fprintf(stderr,"transfer function(%f) = %e\n",k,cd.transfer_function(k));
  
  fprintf(stderr,"\n");
  fprintf(stderr,"sigma(R=8 Mpc/h) = %e\n",cd.sigmaRtophat(8.0,1.0));
  
  FILE *fp;
  double kmin = 1e-6;
  double kmax = 1e4;
  long Nk = 2000;
  double dlnk = log(kmax/kmin)/Nk;
  long i;
  fp = fopen("test.dat","w");
  t = -wtime();
  for(i=0;i<Nk;++i)
    {
      k = exp(dlnk*i)*kmin;
      fprintf(fp,"%e %e %e %e\n",k,cd.transfer_function(k),cd.linear_powspec(k,a),cd.nonlinear_powspec(k,a));
    }
  t += wtime();
  fprintf(stderr,"trans func loop took %g seconds.\n",t);
  fclose(fp);
  */
  return 0;
}
