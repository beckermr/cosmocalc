#include "cosmocalc.h"
#include "cosmocalc_assert.h"

void cosmoCalc::init_cosmology(double omegam, double omegal, double omegab, double omeganu,
			       double hubble, double sigma8, double spectral_index,
			       double w0, double wa, int transfer_function_type)
{
  //set parameters                                                                                                                                                                                       
  _om = omegam;
  _ol = omegal;
  _ob = omegab;
  _onu = 0.0;//omeganu;                                                                                                                                                                                  
  _h = hubble;
  _s8 = sigma8;
  _ns = spectral_index;
  _w0 = w0;
  _wa = wa;
  _ok = 1.0 - _om - _ol - _onu; //curvature density is a derived parameter                                                                                                                               
  DH_sqrtok = DH/sqrt(fabs(_ok));
  
  _growth_function_norm = -1.0;
  _linear_powspec_norm = -1.0;
  
  ++_cosmo_num;

  if(transfer_function_type == COSMOCALC_TRANS_FUNC_EH98)
    {
      _transfer_function_type = COSMOCALC_TRANS_FUNC_EH98;
      _transfer_function_function = &cosmoCalc::transfer_function_eh98;
    }
  else if(transfer_function_type == COSMOCALC_TRANS_FUNC_EH98_SMOOTH)
    {
      _transfer_function_type = COSMOCALC_TRANS_FUNC_EH98_SMOOTH;
      _transfer_function_function = &cosmoCalc::transfer_function_eh98_smooth;
    }
  else
    {
      cosmocalc_assert(_transfer_function_type == transfer_function_type,"invalid transfer function type specified! transfer_function_type = %d\n",transfer_function_type);
    }
};

//functions to init interpolation tables - make sure to catch dependencies
void cosmoCalc::init_distances(void)
{
  if(_cosmo_num != _dist_cosmo_num)
    {
      init_cosmocalc_comvdist_table();
      _dist_cosmo_num = _cosmo_num;
    }
};

void cosmoCalc::init_growth_function(void)
{
  if(_cosmo_num != _gf_cosmo_num)
    {
      init_cosmocalc_growth_function_table();
      _gf_cosmo_num = _cosmo_num;
    }
};

void cosmoCalc::init_transfer_function(void)
{
  if(_cosmo_num != _tf_cosmo_num)
    {
      init_cosmocalc_transfer_function_table();
      _tf_cosmo_num = _cosmo_num;
    }
};

void cosmoCalc::init_linear_powspec(void)
{
  init_growth_function();
  init_transfer_function();
  if(_cosmo_num != _pl_cosmo_num)
    {
      init_cosmocalc_linear_powspec_table();
      _pl_cosmo_num = _cosmo_num;
    }
};

void cosmoCalc::init_nonlinear_powspec(void)
{
  init_growth_function();
  init_transfer_function();
  init_linear_powspec();
  if(_cosmo_num != _pnl_cosmo_num)
    {
      init_cosmocalc_nonlinear_powspec_table();
      _pnl_cosmo_num = _cosmo_num;
    }
};

void cosmoCalc::init_peakheight(void)
{
  init_growth_function();
  init_linear_powspec();
  if(_cosmo_num != _ph_cosmo_num)
    {
      init_cosmocalc_peakheight_table();
      _ph_cosmo_num = _cosmo_num;
    }
}

void cosmoCalc::init_all(void)
{
  init_distances();
  init_growth_function();
  init_transfer_function();
  init_linear_powspec();
  init_nonlinear_powspec();
  init_peakheight();
};
