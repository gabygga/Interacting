import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math
import scienceplots

plt.style.use(['science','ieee'])

k_out = [0.00067, 0.067, 0.743 ] # 1/Mpc
#

# Cosmological parameters and other CLASS parameters, updated to Planck 2018, suma de masas 0.06eV

common_settings = {# we need to set the output field to something although
                   # the really releveant outpout here will be set with 'k_output_values'
                   'output':'mPk,tCl,pCl,lCl,mTk',
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
MQ.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 0., 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ2=Class()
MQ2.set(common_settings)
MQ2.set({'a_ini_over_a_today_default':1.e-14,'beta':1.e-6, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24.,0., 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ3=Class()
MQ3.set(common_settings)

MQ3.set({'a_ini_over_a_today_default':1.e-14, 'beta':5.e-6, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 0., 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})
MQ4=Class()
MQ4.set(common_settings)
MQ4.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 0., 1.e-2, 1.e-16, 1.e-30',
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
Transfers=M.get_transfer()
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
a_LCDM_k0 = one_k[0]['a']
delta_LCDM_k0 = one_k[0]['delta_cdm']
a_z=1/(1+background['z'])
rho_cdm= interp1d(a_z, background['(.)rho_cdm'])

a_beta0_k0 = one_kQ[0]['a']
delta_sfdm_k0_beta0 = one_kQ[0]['delta_sfdm_1']
##
a_beta1_k0 = one_kQ2[0]['a']
delta_sfdm_k0_beta1 = one_kQ2[0]['delta_sfdm_1']
# plotting
#################
fig, (k0, k1, k2) = plt.subplots(3,sharex=True,figsize=(3.3, 3.8))
fig.subplots_adjust(hspace=0)
#plt.figure(figsize=(3.3, 2.5))

k0.plot(a_LCDM_k0, np.abs(delta_LCDM_k0),'k' , linestyle='--', label='CDM')
#k1.plot(one_k[1]['a'], np.abs(one_k[1]['delta_cdm']),'k' , linestyle='--', label='$CDM$')
#k2.plot(one_k[2]['a'], np.abs(one_k[2]['delta_cdm']),'k',linestyle='--', label='CDM')

k0.plot(one_kQ4[0]['a'], np.abs(one_kQ4[0]['delta_sfdm_1']), 'c', label='$\phi, \lambda_{\phi}=0$')
k1.plot(one_kQ4[1]['a'], np.abs(one_kQ4[1]['delta_sfdm_1']), 'c', label='$\phi$')
k2.plot(one_kQ4[2]['a'], np.abs(one_kQ4[2]['delta_sfdm_1']), 'c',label='$\phi$')

k0.plot(one_kQ[0]['a'], np.abs(one_kQ[0]['delta_sfdm_1']),'blue',alpha=0.6,linestyle='--',label='$\\beta=0$')
k1.plot(one_kQ[1]['a'], np.abs(one_kQ[1]['delta_sfdm_1']),'blue',alpha=0.6,linestyle='--',label='$\\beta=0$')
k2.plot(one_kQ[2]['a'], np.abs(one_kQ[2]['delta_sfdm_1']),'blue',alpha=0.6,linestyle='--',label='$\\beta=0$')

k0.plot(one_kQ2[0]['a'], np.abs(one_kQ2[0]['delta_sfdm_1']),'deeppink',alpha=0.9,linestyle=':',label='$\\beta= 10^{-6}$')
k1.plot(one_kQ2[1]['a'], np.abs(one_kQ2[1]['delta_sfdm_1']),'deeppink',alpha=0.9,linestyle=':',label='$\\beta= 10^{-6}$')
k2.plot(one_kQ2[2]['a'], np.abs(one_kQ2[2]['delta_sfdm_1']),'deeppink',alpha=0.9,linestyle=':',label='$\\beta= 10^{-6}$')


k0.plot(one_k[0]['a'], np.abs(one_k[0]['delta_cdm']),'k',linestyle='--',)
k1.plot(one_k[1]['a'], np.abs(one_k[1]['delta_cdm']),'k' , linestyle='--', label='$CDM$')
k2.plot(one_k[2]['a'], np.abs(one_k[2]['delta_cdm']),'k',linestyle='--', label='CDM')

#plt.plot(one_k[1]['a'], np.abs(one_k[1]['delta_cdm']),'k' , linestyle='--', label='$CDM$')
#plt.plot(one_kQ4[1]['a'], np.abs(one_kQ4[1]['delta_sfdm_1']), 'c', label='$\phi$')
#plt.plot(one_kQ[1]['a'], np.abs(one_kQ[1]['delta_sfdm_1']),'blue',alpha=0.6,linestyle='--',label='$\\beta=0$')
#plt.plot(one_kQ2[1]['a'], np.abs(one_kQ2[1]['delta_sfdm_1']),'deeppink',alpha=0.9,linestyle=':',label='$\\beta= 10^{-6}$')
#plt.plot(one_k[1]['a'], np.abs(one_k[1]['delta_cdm']),'k' , linestyle='--', label='$CDM$')
#plt.plot(one_kQ[1]['a'], np.abs(one_kQ[1]['delta_sfdm_1']),'blue',alpha=0.6,linestyle='--',label='$\\beta=0$')

#
#plt.legend(title='Albrecht Skordis')
plt.xscale('log')
plt.yscale('log')
k0.set_yscale('log')
k1.set_yscale('log')
#k2.set_yscale('log')
k0.legend(fontsize=8,loc='lower right')
k0.set_title('$k=10^{-3} h/Mpc $',fontsize=8, x=0.48, y=0.04)
k1.set_title('$k=0.1 h/Mpc$',fontsize=8, x=0.48, y=0.04)
k2.set_title('$k=1.1 h/Mpc$',fontsize=8, x=0.48, y=0.04)

k0.set_ylim([6.e-6,50])
k1.set_ylim([2.e-5,8.e4])
k2.set_ylim([2.e-5,3.3e5])

plt.xlim([4.e-6,1.09])
plt.xlabel(r'$a$')
#plt.legend()
k0.set_ylabel(r'$|\delta|$')
k1.set_ylabel(r'$|\delta|$')
k2.set_ylabel(r'$|\delta|$')
plt.savefig('scripts_int/Plots/Contrast_3k_quad_log.pdf')
    
