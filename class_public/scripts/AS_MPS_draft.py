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
  
#Defining an array of k-values
kvec = np.logspace(-4,np.log10(3),1000)

pkM = []
pkMQ= []
pkMQ2 = []
pkMQ3 = []
pkMQ4 = []
for k in kvec:
    pkM.append(M.pk(k,0.))
    pkMQ.append(MQ.pk(k,0.))
    pkMQ2.append(MQ2.pk(k,0.))
    pkMQ3.append(MQ3.pk(k,0.))
    pkMQ4.append(MQ4.pk(k,0.))
h = M.h()
# plotting
#################
#
#plt.loglog(kvec/h,np.array(pkM),'r',linestyle='-', label='$\Lambda$CDM', lw='2.5') 

#plt.loglog(kvec/h, np.array(pkMQ),'b',linestyle='-',label='AS') 

plt.figure(figsize=(3.3,2.5)) 
plt.semilogx(kvec/h,  np.array(pkMQ2)/np.array(pkM),color='m',alpha=0.6,linestyle=':', label='$\lambda=8, A=0.01, B=34.80$')
plt.semilogx(kvec/h, np.array(pkMQ3)/np.array(pkM),color='m',alpha=0.8,linestyle='--', label='$\lambda=10, A=0.008, B=27.20$')
plt.semilogx(kvec/h, np.array(pkMQ4)/np.array(pkM),color='m',linestyle='-', label='$\lambda=12, A=0.004, B=22.66$')

plt.rcParams["figure.autolayout"] = True

#
#plt.legend(title='Albrecht Skordis')
plt.xlim([0.00012,3])
plt.legend()
plt.ylim([0.717,1.08])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$', labelpad = 1)
plt.ylabel(r'$P^{\rm AS}/P_{\Lambda\rm CDM}$')
#plt.tick_params(which='minor',axis='both',direction='in',right=True,top=True)
#plt.tick_params(which='major',axis='both',direction='in',right=True,top=True)

plt.savefig('scripts/Plots/AS_MPS_ratio.pdf', bbox_inches='tight', pad_inches=0.04)
    

