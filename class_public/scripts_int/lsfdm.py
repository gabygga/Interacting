import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math

k_out = [0.1] # 1/Mpc
#

# Cosmological parameters and other CLASS parameters, updated to Planck 2018, suma de masas 0.06eV

common_settings = {# we need to set the output field to something although
                   # the really releveant outpout here will be set with 'k_output_values'
                   'output':'mPk,tCl,pCl,lCl',
                   'P_k_max_1/Mpc':3.0,
                   # value of k we want to polot in [1/Mpc]
                   'k_output_values': str(k_out).strip('[]'),
                   # LambdaCDM parameters
                   'h':0.6732117,
                   'omega_b':0.02238280,
                   'Omega_cdm':0.264,
                   'A_s':2.100549e-09 ,
                   'n_s':0.9660499,
                   'tau_reio':0.05430842,
                   'N_ur':0.00641,
                   'T_cmb':2.726,
                   'N_ncdm':1,
                   'deg_ncdm':3,
                   'm_ncdm': 0.02,
                # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe': 0.2454006,
                   #'Omega_Lambda': ,
                   'Omega_scf' :0.,
                   'Omega_fld':0,
                   #'cs2_fld' : 1,
               # other options and settings
                   'compute damping scale':'yes', # needed to output the time of damping scale crossing
                   'lensing':'yes'
                    }
##############
#
# call CLASS
#
M = {}
MQ={}
MQ2={}
MQ3={}
all_k = {}
one_k = {}
background= {}
all_kQ = {}
one_kQ = {}
backgroundQ= {}
all_kQ2 = {}
one_kQ2 = {}
backgroundQ2= {}
all_kQ3 = {}
one_kQ3 = {}
backgroundQ3= {}


    # call CLASS
    #
M = Class()
M.set(common_settings)
MQ=Class()
MQ.set(common_settings)
MQ.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.25,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-22., 0.0, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2})
M.compute()
MQ.compute()
#load perturbations
all_k=M.get_perturbations()
one_k=all_k['scalar']
background = M.get_background()
    ##
all_kQ=MQ.get_perturbations()
one_kQ=all_kQ['scalar']
backgroundQ = MQ.get_background()
 ###

###
#CMB lensed spectra    
cl_lens0= M.lensed_cl(2500)
cl_lens1= MQ.lensed_cl(2500)
# plotting
#################
#
fig, (D_ell, DeltaD, DeltaD2) = plt.subplots(3,sharex=True,figsize=(8,12))
fig.subplots_adjust(hspace=0)
plt.xlim([2,2510])
#
ell_l=cl_lens0['ell']
#Factor for units
factor_l = (2.726e6)*(2.726e6)*ell_l**2*ell_l*(ell_l+1)/2/np.pi
# D_ell suma de masas 0.06
D_ell.plot(ell_l,1e-9*factor_l*cl_lens0['tt'],lw ='2.', color='r',label=r'$\Lambda CDM$')
D_ell.plot(ell_l,1e-9*factor_l*cl_lens1['tt'],lw ='2.',color='b',label=r'$m=-22$')
### Delta D_ell, suma de masas 0.06
DeltaD.plot(ell_l,1e-8*factor_l*(cl_lens1['tt']-cl_lens0['tt']),lw ='2.', color='b')

## Delta D_ell/Delta_LCDM: LCDM es considerado con suma de masas 0.06
DeltaD2.plot(ell_l,(cl_lens1['tt']-cl_lens0['tt'])/cl_lens0['tt'],lw ='2.',color='b')

### Label options
D_ell.legend(title='sfdm', fontsize='10',ncol=2)
#D_ell.set_xticks([ 500, 1000, 1500, 2000],["500", "1000", "1500", "2000"])
D_ell.tick_params(direction='in', top=True, right=True, labeltop=False,length=7,width=1.5)
D_ell.set_ylabel(r"$\ell^2 D_{\ell} [\mu K^2] \,\,\, (\times 10^{-9})$", fontsize=16)
D_ell.set_ylim([-0.495, 2.2])
D_ell.ticklabel_format(useMathText=True)
DeltaD.set_ylabel(r"$\Delta (\ell^2 D_{\ell}) [\mu K^2] \,\,\,(\times 10^{-8}) $", fontsize=16)
DeltaD.tick_params(direction='in', top=True, right=True, labeltop=False,length=7,width=1.5)
DeltaD.ticklabel_format(useMathText=True)
DeltaD2.set_ylabel(r"$\Delta  D_{\ell}/ D_{\ell}^{\Lambda\rm CDM}  \,\,\, $", fontsize=16)
DeltaD2.set_xlabel(r"$\ell  \,\,\, $", fontsize=16)
DeltaD2.tick_params(direction='in', top=True, right=True, labeltop=False,length=7,width=1.5)
DeltaD2.legend()
plt.savefig('lsfdm.pdf',bbox_inches='tight')

    
