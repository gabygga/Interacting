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
                   'beta':0,
                   # LambdaCDM parameters
                   'h':0.6732117,
                   'omega_b':0.02238280,
                   #'omega_cdm':0.1201075,
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
MQ1={}
MQ2={}
MQ3={}
all_k = {}
one_k = {}
background= {}
all_kQ = {}
one_kQ = {}
backgroundQ= {}
all_kQ1 = {}
one_kQ1 = {}
backgroundQ1= {}
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
##
MQ=Class()
MQ.set(common_settings)
MQ.set({'a_ini_over_a_today_default':1.e-14, 'beta':0,'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ2=Class()
MQ2.set(common_settings)
MQ2.set({'a_ini_over_a_today_default':1.e-14,'beta':1.e-6, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ3=Class()
MQ3.set(common_settings)

MQ3.set({'a_ini_over_a_today_default':1.e-14, 'beta':5.e-6, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ4=Class()
MQ4.set(common_settings)
MQ4.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,})

M.compute()
MQ3.compute()

MQ2.compute()
MQ.compute()
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

#Defining an array of k-values
kvec = np.logspace(-4,np.log10(3),1000)

pkM   = []
pkMQ  = []
pkMQ3 = []
pkMQ2 = []
pkMQ4 = []

for k in kvec:
    pkM.append(M.pk(k,0.))
    pkMQ.append(MQ.pk(k,0.))
    pkMQ3.append(MQ3.pk(k,0.))
    pkMQ2.append(MQ2.pk(k,0.))
    pkMQ4.append(MQ4.pk(k,0.))
    h = M.h()
# plotting
#################
#
plt.plot(kvec/h, np.array(pkMQ4)/np.array(pkM),'c',  linestyle='-',label='$\phi, \lambda_{\phi}>0$')

plt.plot(kvec/h, np.array(pkMQ)/np.array(pkM),'b' ,alpha=0.6, linestyle='--', label='$\\beta=0$')

plt.plot(kvec/h, np.array(pkMQ2)/np.array(pkM),'deeppink',alpha=0.9,linestyle=':',label='$\\beta= 10^{-6}$')

#plt.plot(kvec/h, np.array(pkMQ3)/np.array(pkM),'r', alpha=0.6, linestyle='-.',label='$\\beta=5\\times 10^{-6}$')


#plt.plot(kvec/h, (kvec/h)**3*np.array(pkM),'k' ,alpha=0.6, linestyle='--', label='$\Lambda$')

#plt.plot(kvec/h, (kvec/h)**3*np.array(pkMQ),'b' ,alpha=0.6, linestyle='--', label='$\\beta=0$')

#plt.plot(kvec/h, (kvec/h)**3*np.array(pkMQ2),'lightseagreen',alpha=0.6,linestyle='-.',label='$\\beta= 10^{-6}$')

#plt.plot(kvec/h, (kvec/h)**3*np.array(pkMQ3),'lightseagreen', alpha=0.6, linestyle=':',label='$\\beta=5\\times 10^{-6}$')


#
#plt.legend(title='Albrecht Skordis')
plt.xscale('log')
#plt.yscale('log')
plt.legend()
#plt.ylim([2e-7,2e5])
plt.xlim([2.e-4,6])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P^{\phi,\psi}/P_{\Lambda\rm CDM}$')

plt.savefig('scripts_int/Plots/MPS_AS_SF_beta_trig_lin.pdf')
    
    
