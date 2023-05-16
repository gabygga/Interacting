import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
import scienceplots

plt.style.use(['science','ieee'])

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
                   'a_ini_over_a_today_default': 1.e-14,
                   'h':0.6732117,
                   'omega_b':0.02238280,
                   'omega_cdm':0.1201075,
                   'A_s':2.100549e-09 ,
                   'n_s':0.9660499,
                   'tau_reio':0.05430842,
                   'N_ur':0.00641,
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
MQ.set({'Omega_Lambda':1e-7,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '9, 2.0, 0.01, 30.2203, 19, 0.0',
                   'scf_tuning_index':0,})

MQ2=Class()
MQ2.set(common_settings)
MQ2.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '7.81, 2.0, 0.01, 34.8, 22, 0.0',
                   'scf_tuning_index':0,})
MQ3=Class()
MQ3.set(common_settings)
MQ3.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '10, 2.0, 0.008,27.2, 17.2, 0.0',
                   'scf_tuning_index':0,})
MQ4=Class()
MQ4.set(common_settings)
MQ4.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '12, 2.0, 0.004,22.66, 14.2, 0.0',
                   'scf_tuning_index':0,})

M.compute()
MQ.compute()
MQ2.compute()
MQ3.compute()
MQ4.compute()

    #load perturbations
all_k=M.get_perturbations()
one_k=all_k['scalar']
background = M.get_background()
    ##
all_kQ=MQ.get_perturbations()
one_kQ=all_kQ['scalar']
backgroundQ = MQ.get_background()
 ###
all_kQ2=MQ2.get_perturbations()
one_kQ2=all_kQ2['scalar']
backgroundQ2 = MQ2.get_background()
  ###
all_kQ3=MQ3.get_perturbations()
one_kQ3=all_kQ3['scalar']
backgroundQ3 = MQ3.get_background()
  ###
all_kQ4=MQ4.get_perturbations()
one_kQ4=all_kQ4['scalar']
backgroundQ4 = MQ4.get_background()
     ###
#CMB lensed spectra    
cl_lens0= M.lensed_cl(2500)
cl_lens1= MQ.lensed_cl(2500)
cl_lens2= MQ2.lensed_cl(2500)
cl_lens3= MQ3.lensed_cl(2500)
cl_lens4= MQ4.lensed_cl(2500)
ell_l=cl_lens0['ell']

# plotting
#################
#
plt.figure(figsize=(3.3,2.5)) 
#plt.loglog(kvec/h, np.array(pkMQ),'b',linestyle='-',label='AS') 
plt.plot(ell_l,(cl_lens2['tt']-cl_lens0['tt'])/cl_lens0['tt'],  color='m',alpha=0.6, linestyle=':',label='$\lambda=8, A=0.01, B=34.80$')
plt.plot(ell_l,(cl_lens3['tt']-cl_lens0['tt'])/cl_lens0['tt'],  color='m',alpha=0.75,linestyle='--',label='$\lambda=10, A=0.008, B=27.20$')
plt.plot(ell_l,(cl_lens4['tt']-cl_lens0['tt'])/cl_lens0['tt'], color='m',linestyle='-', label='$\lambda=12, A=0.004, B=22.66$')

#plt.legend(title='Albrecht Skordis')
plt.xlim([1.5,2560])
plt.ylim([-0.04,0.215 ])
plt.legend()

plt.xlabel(r'$\ell$')
plt.ylabel(r'$\Delta C_{\ell}/C_{\ell}^{\Lambda \rm CDM}$')
#plt.tick_params(which='minor',axis='both',direction='in',right=True,top=True)
#plt.tick_params(which='major',axis='both',direction='in',right=True,top=True)

plt.savefig('scripts/Plots/AS_TT3.pdf', bbox_inches='tight', pad_inches=0.04)
    

