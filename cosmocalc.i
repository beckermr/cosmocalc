/* File : cosmocalc.i */
%module(docstring="Cosmology and LSS calculations - python interface to (fast) C code") cosmocalc
%{
/* Put headers and other declarations here */
#define SWIG_FILE_WITH_INIT
#include "cosmocalc.h"
#include "haloprofs.h"
#include "weaklens.h"
%}

/* in global.c */
%feature("autodoc", "turn off GSL error handling") turn_off_gsl_errs;

/* in distances.c */
%feature("autodoc", "comoving distance in Mpc/h - comvdist(scale factor)") comvdist;
%feature("autodoc", "comoving distance in Mpc/h - comvdist_exact(scale factor) [does integration as opposed to using spline]") comvdist_exact;
%feature("autodoc", "angular diameter distance in Mpc/h - angdist(scale factor)") angdist;
%feature("autodoc", "luminosity distance in Mpc/h - lumdist(scale factor)") lumdist;
%feature("autodoc", "scale factor for a given comoving distance in Mpc/h - acomvdist(comv. dist in Mpc/h)") acomvdist;
%feature("autodoc", "angular diatameter distance difference in Mpc/h for lensing - angdistdiff(amin,amax) [amin <= amax]") angdistdiff;

/* ages.c */
%feature("autodoc", "age in yr - age(scale factor)") age;
%feature("autodoc", "lookback time in yr - lookback(scale factor)") lookback;

/* in transfer_function.c */
%feature("autodoc", "EH98 transfer function - transfer_function(k in h/Mpc) [uses spline]") transfer_function;
%feature("autodoc", "EH98 transfer function - transfunct_eh98(k in h/Mpc) [evaluates fitting formula - slower than transfer_function]") transfunct_eh98;
%feature("autodoc", "EH98 transfer function w/ no BAO wiggles - transfunct_eh98_smooth(k in h/Mpc)") transfunct_eh98_smooth;

/* in linear_powspec.c */
%feature("autodoc", "linear power spectrum P(k) - linear_powspec(k in h/Mpc, scale factor)") linear_powspec;
%feature("autodoc", "linear power spectrum P(k) - linear_powspec_exact(k in h/Mpc, scale factor) [does integration as opposed to using spline]") linear_powspec_exact;
%feature("autodoc", "convert CMB power spectrum amplitude to sigma8 - convert_cmbnorm2sigma8()") convert_cmbnorm2sigma8;

/* linear_corrfunc.c */
%feature("autodoc", "linear corr. function xi(r) - linear_corrfunc(r in Mpc/h, scale factor)") linear_corrfunc;
%feature("autodoc", "linear corr. function xi(r) - linear_corrfunc(r in Mpc/h, scale factor) [does integration as opposed to using spline]") linear_corrfunc_exact;

/* in nonlinear_powspec.c */
%feature("autodoc", "Takahashi+12 HaloFit nonlinear power spectrum - nonlinear_powspec(k in h/Mpc, scale factor)") nonlinear_powspec;

/* nonlinear_corrfunc.c */
%feature("autodoc", "Takahashi+12 HaloFit nonlinear corr. function xi(r) - nonlinear_corrfunc(r in Mpc/h, scale factor)") nonlinear_corrfunc;
%feature("autodoc", "Takahashi+12 HaloFit nonlinear corr. function xi(r) - nonlinear_corrfunc(r in Mpc/h, scale factor) [does integration as opposed to using spline]") nonlinear_corrfunc_exact;

/* in growth_function.c */
%feature("autodoc", "growth function ODE integration - growth_function(scale factor)") growth_function;
%feature("autodoc", "growth function ODE integration - growth_function_exact(scale factor) [does ODE as opposed to using spline]") growth_function_exact;
%feature("autodoc", "growth function ODE integration w/ no normalization - growth_function_exact_nonorm(scale factor) [does ODE as opposed to using spline]") growth_function_exact_nonorm;

/* in hubble.c */
%feature("autodoc", "hubble function H(a)/H0 - hubble_noscale(scale factor)") hubble_noscale;
%feature("autodoc", "effective w for DE so that rho_{DE} = Omega_{DE}(a=1)/a^(3*(1+weff(a))) - weff(scale factor)") weff;

/* in peakheight.c */
%feature("autodoc", "RMS variance in top hat sphere of radius R of linear power spectrum - sigmaRtophat(radius in Mpc/h, scale factor)") sigmaRtophat;
%feature("autodoc", "RMS variance in top hat sphere of radius R of linear power spectrum - sigmaRtophat(radius in Mpc/h, scale factor) [does integration as opposed to using spline]") sigmaRtophat_exact;
%feature("autodoc", "RMS variance for mass M of linear power spectrum - sigmaMtophat(mass in M_{sun}/h, scale factor)") sigmaMtophat;
%feature("autodoc", "RMS variance for mass M of linear power spectrum - sigmaMtophat(mass in M_{sun}/h, scale factor) [does integration as opposed to using spline]") sigmaMtophat_exact;
%feature("autodoc", "sphere of radius R given RMS variance in linear power spectrum - inverse_sigmaRtophat(RMS variance sigma(R,scale factor), scale factor)") inverse_sigmaRtophat;
%feature("autodoc", "mass M given RMS variance in linear power spectrum - inverse_sigmaMtophat(RMS variance sigma(M,scale factor), scale factor)") inverse_sigmaMtophat;

/* in mass_bias_functions.c */
%feature("autodoc", "Tinker+08 mass function - tinker2008_mass_function(mass, scale factor, delta w/ mean density)") tinker2008_mass_function;
%feature("autodoc", "Tinker+10 halos bias - tinker2010_bias(mass, scale factor, delta w/ mean density)") tinker2010_bias;
%feature("autodoc", "Tinker+10 mass function - tinker2010_mass_function(mass, scale factor, delta w/ mean density)") tinker2010_mass_function;

%include "cosmocalc_swig.h"

%pythoncode %{
_cosmocalc.cvar.cosmoData.cosmoNum = 1
_cosmocalc.turn_off_gsl_errs()
_cosmocalc.cvar.cosmoData.useSmoothTransFunc = 0
_cosmocalc.cvar.cosmoData.useNoTransFunc = 0

def useEH98SmoothTransFunc():
    _cosmocalc.cvar.cosmoData.useSmoothTransFunc = 1
    _cosmocalc.cvar.cosmoData.useNoTransFunc = 0
    _cosmocalc.cvar.cosmoData.cosmoNum += 1

def useNoTransFunc():
    _cosmocalc.cvar.cosmoData.useSmoothTransFunc = 0
    _cosmocalc.cvar.cosmoData.useNoTransFunc = 1
    _cosmocalc.cvar.cosmoData.cosmoNum += 1

def useEH98TransFunc():
    _cosmocalc.cvar.cosmoData.useSmoothTransFunc = 0
    _cosmocalc.cvar.cosmoData.useNoTransFunc = 0
    _cosmocalc.cvar.cosmoData.cosmoNum += 1

def _init(cd):
    _cosmocalc.cvar.cosmoData.cosmoNum += 1
    _cosmocalc.cvar.cosmoData.OmegaM  = _cosmodict_resolve(cd,'om')
    _cosmocalc.cvar.cosmoData.OmegaL  = _cosmodict_resolve(cd,'ol')
    _cosmocalc.cvar.cosmoData.OmegaB  = _cosmodict_resolve(cd,'ob')
    okval = _cosmodict_resolve(cd,'ok')
    if okval is None:
        _cosmocalc.cvar.cosmoData.OmegaK = 1.0 - _cosmocalc.cvar.cosmoData.OmegaM - _cosmocalc.cvar.cosmoData.OmegaL
    else:
        _cosmocalc.cvar.cosmoData.OmegaK = okval
    _cosmocalc.cvar.cosmoData.h = _cosmodict_resolve(cd,'h')
    _cosmocalc.cvar.cosmoData.SpectralIndex = _cosmodict_resolve(cd,'ns')

    _cosmocalc.cvar.cosmoData.w0 = -1.0
    _cosmocalc.cvar.cosmoData.wa = 0.0
    if _cosmodict_resolve(cd,'w') is not None:
        _cosmocalc.cvar.cosmoData.w0 = _cosmodict_resolve(cd,'w')
        _cosmocalc.cvar.cosmoData.wa = 0.0
    elif _cosmodict_resolve(cd,'w0') is not None and _cosmodict_resolve(cd,'wa') is not None:
        _cosmocalc.cvar.cosmoData.w0 = _cosmodict_resolve(cd,'w0')
        _cosmocalc.cvar.cosmoData.wa = _cosmodict_resolve(cd,'wa')

    asval = _cosmodict_resolve(cd,'as')
    as_pivot_val = _cosmodict_resolve(cd,'as_pivot')
    s8val = _cosmodict_resolve(cd,'s8')
    if asval is not None and as_pivot_val is not None and s8val is None:
        _cosmocalc.cvar.cosmoData.As = asval
        _cosmocalc.cvar.cosmoData.As_pivot = as_pivot_val
	_cosmocalc.cvar.cosmoData.Sigma8 = _cosmocalc.convert_cmbnorm2sigma8()
    elif s8val is not None:
        _cosmocalc.cvar.cosmoData.Sigma8 = s8val

def _cosmodict_resolve(cd,skey,keynames=None):
    _keynames = [
        ['omega_m','om','OmegaM','Omega_M','omegam'],
        ['omega_l','ol','OmegaL','Omega_L','omegal',
         'omega_de','ode','OmegaDE','Omega_DE','omegade'],
        ['omega_k','ok','OmegaK','Omega_K','omegak'],
        ['omega_r','or','OmegaR','Omega_R','omegar'],
        ['omega_nu','onu','OmegaNu','Omega_Nu','omeganu',
         'OmegaNU','Omega_NU'],
        ['omega_b','ob','OmegaB','Omega_B','omegab'],
        ['H0'],
        ['h','hubble'],
        ['s8','sigma8','Sigma8','sigma_8','Sigma_8'],
        ['as','As'],
	['as_pivot','As_pivot','AsPivot','as_Pivot','As_Pivot'],
        ['ns','SpectralIndex','spectral_index','spectralindex'],
        ['w','W'],
        ['w0','W0'],
        ['wa','wA','WA']
        ]
    if keynames is None:
        keynames = _keynames
    for keylist in keynames:
        if skey in keylist:
            for tkey in keylist:
                if tkey in cd:
                    return cd[tkey]
    return None

def set_cosmology(cd):
    """Set the cosmology using a dictionary like this

           import cosmocalc
           cd = {
               "OmegaM":0.3,
               "OmegaB":0.045,
               "OmegaDE":0.7,
               "OmegaK":0.0,
               "h":0.7,
               "Sigma8":0.8,
               "SpectralIndex":0.95,
               "w0":-1.0,
               "wa":0.0
               "As":2.1e-9,
               "As_pivot":0.05}
           cosmocalc.set_cosmology(cd)
       
       Note that the cosmology is set *globally*, so you can only use 1 at a time!
       """
    _init(cd)

def get_cosmology():
    "Return the cosmology as a dictionary. See cosmocalc.set_cosmology() for an example."
    return {'OmegaM':_cosmocalc.cvar.cosmoData.OmegaM,
            'OmegaDE':_cosmocalc.cvar.cosmoData.OmegaL,
            'OmegaB':_cosmocalc.cvar.cosmoData.OmegaB,
            'OmegaK':_cosmocalc.cvar.cosmoData.OmegaK,
            'h':_cosmocalc.cvar.cosmoData.h,
            'Sigma8':_cosmocalc.cvar.cosmoData.Sigma8,
            'SpectralIndex':_cosmocalc.cvar.cosmoData.SpectralIndex,
            'w0':_cosmocalc.cvar.cosmoData.w0,
            'wa':_cosmocalc.cvar.cosmoData.wa}

def lcdm():
    "Return a random cosmology."
    return {'OmegaM':0.3,'OmegaDE':0.7,'OmegaB':0.045,'OmegaK':0.0,'h':0.7,'Sigma8':0.8,'SpectralIndex':0.95,'w0':-1.0,'wa':0.0}

%}


