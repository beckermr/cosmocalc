#ifndef _HALOPROFS_
#define _HALOPROFS_

double fourierTransformNFW(double k, double m, double rvir, double c);
double NFWprof(double r, double m, double rvir, double c);
double NFWprof_menc(double r, double m, double rvir, double c);
double concNFW(double m, double a);
double duffy2008_concNFW(double m, double a);

#endif /* _HALOPROFS_ */
