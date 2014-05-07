#ifndef _TF_FIT_EH98_
#define _TF_FIT_EH98_

void TFset_parameters(double omega0hh, double f_baryon, double Tcmb);
double TFfit_onek(double k, double *tf_baryon, double *tf_cdm);

void TFfit_hmpc(double omega0, double f_baryon, double hubble, double Tcmb,
		int numk, double *k, double *tf_full, double *tf_baryon, double *tf_cdm);

double TFsound_horizon_fit(double omega0, double f_baryon, double hubble);
double TFk_peak(double omega0, double f_baryon, double hubble);
double TFnowiggles(double omega0, double f_baryon, double hubble,
		  double Tcmb, double k_hmpc);
double TFzerobaryon(double omega0, double hubble, double Tcmb, double k_hmpc);

#endif /* _TF_FIT_EH98_ */
