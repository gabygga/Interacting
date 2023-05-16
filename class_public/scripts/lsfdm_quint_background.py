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
                   #'h':0.6732117,
                   'h':0.72,
                   'beta':0,
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
MQ.set({'a_ini_over_a_today_default':1.e-15, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.25,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-24., 0.0, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '7.812, 2.0, 0.01, 34.8, 21, 0.0',
                   'scf_tuning_index':0,})


M.compute()
MQ.compute()

all_k=M.get_perturbations()
one_k=all_k['scalar']
background = M.get_background()
    ##
all_kQ=MQ.get_perturbations()
one_kQ=all_kQ['scalar']
backgroundQ = MQ.get_background()

    ###
#  print(background.keys())
#rho_scf=background['rho_scf'] 
#################
#
baH = background['H [1/Mpc]']
baT = background['conf. time [Mpc]']
baa = 1/(1 + background['z'])
#baCC= background['(.)rho_lambda']
baCrit = background['(.)rho_crit']
##
baCritQ = backgroundQ['(.)rho_crit']
rho_sfdmQ=backgroundQ['(.)rho_sfdm_1']
p_sfdmQ=backgroundQ['(.)p_sfdm_1']
rho_gQ =backgroundQ['(.)rho_g']
#rho_lambdaQ=backgroundQ['(.)rho_lambda']
rho_bQ =backgroundQ['(.)rho_b']
##

plt.semilogx(1/(1+backgroundQ['z']), (backgroundQ['(.)rho_g']+backgroundQ['(.)rho_ncdm[0]'])/baCritQ, color='b',linestyle='-', label='$\Omega_{rad}$')
plt.semilogx(1/(1+backgroundQ['z']), (backgroundQ['(.)rho_b']+backgroundQ['(.)rho_cdm'])/baCritQ, color='g',linestyle='--', label='$\Omega_{b}$')
plt.semilogx(1/(1+backgroundQ['z']), (backgroundQ['(.)rho_scf'])/baCritQ, color='c',linestyle='-', label='$\Omega_{quint}$')
plt.semilogx(1/(1+backgroundQ['z']), (backgroundQ['(.)rho_sfdm_1'])/baCritQ, color='m',linestyle='-', label='$\Omega_{sfdm}$')
#plt.semilogx(1/(1+backgroundQ['z']), (backgroundQ['(.)p_scf'])/(backgroundQ['(.)rho_scf']), color='c',linestyle='-', label='$\omega_{quint}$')
#plt.semilogx(1/(1+backgroundQ['z']),  (backgroundQ['(.)p_sfdm_1'])/(backgroundQ['(.)rho_sfdm_1']), color='m',linestyle='-', label='$\omega_{sfdm}$')
#
#
plt.xlim([1e-15, 1.1])
#plt.ylim([0.0, 20])

plt.xlabel(r"$a$")
plt.ylabel(r"$\mathrm{\Omega}$")
plt.tight_layout()
plt.legend(title='$m_{\phi}=1e-24$')

# In[ ]:


plt.savefig('scripts/Plots/Omega_sfdm_quint.pdf')

