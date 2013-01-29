#include <cstdlib>
#include <cstdio>
#include "cosmocalc.h"
#include <omp.h>

#include <sys/time.h>
static double wtime(void)
{
  struct timeval tp;
  gettimeofday(&tp,NULL);
  return ((double) (tp.tv_sec)) + ((double) (tp.tv_usec))/1e6;
}

int main(int argc, char **argv)
{
  cosmoCalc cd;
  
  fprintf(stderr,"Testing the cosmology routines...\n\n");
    
  double t = -wtime();
  cd.init_cosmology(0.25,0.75,0.145,0.0,0.7,0.8,1.0,-1.0,0.0,COSMOCALC_TRANS_FUNC_EH98);
  cd.init_all();
  if(argc > 3)
    {
      fprintf(stderr,"setting library to use %d threads!\n",atoi(argv[3]));
      cd.num_threads(atoi(argv[3]));
    }
    
  t += wtime();
  fprintf(stderr,"first init took %g seconds.\n",t);

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
      fprintf(fp,"%e %e %e %e %e\n",k,cd.transfer_function(k),cd.linear_powspec(k,a),cd.linear_powspec(k,a),cd.nonlinear_powspec(k,a));
    }
  t += wtime();
  fprintf(stderr,"trans func loop took %g seconds.\n",t);
  fclose(fp);
  
  return 0;
}
