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
MQ.set({'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '1.4, 0.0, 10, 0.0, 20, 0.0',
                   'scf_tuning_index':0,})
MQ2=Class()
MQ2.set(common_settings)
MQ2.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '1.4, 0.0, 40, 0.0, 20, 0.0',
                   'scf_tuning_index':0,})
MQ3=Class()
MQ3.set(common_settings)
MQ3.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '1.4, 0.0, 0., 0.0, 20, 0.0',
                   'scf_tuning_index':0,})
MQ4=Class()
MQ4.set(common_settings)
MQ4.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '1.4, 1.0, 0., 0.0, 20, 0.0',
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
#  print(background.keys())
#rho_scf=background['rho_scf'] 
#################
#
baH = background['H [1/Mpc]']
baT = background['conf. time [Mpc]']
baa = 1/(1 + background['z'])
baCC= background['(.)rho_lambda']
baCrit = background['(.)rho_crit']
##
bVQ = backgroundQ['V_scf']
bppQ = backgroundQ["phi'_scf"]
baCritQ = backgroundQ['(.)rho_crit']
rho_scfQ = (bppQ*bppQ/(2*baa*baa) + bVQ)/3.
Omega_scfQ = rho_scfQ/baCritQ
##
bVQ2 = backgroundQ2['V_scf']
bppQ2 = backgroundQ2["phi'_scf"]
baCritQ2 = backgroundQ2['(.)rho_crit']
rho_scfQ2 = (bppQ2*bppQ2/(2*baa*baa) + bVQ2)/3.
Omega_scfQ2 = rho_scfQ2/baCritQ2
##
bVQ3 = backgroundQ3['V_scf']
bppQ3 = backgroundQ3["phi'_scf"]
baCritQ3 = backgroundQ3['(.)rho_crit']
rho_scfQ3 = (bppQ3*bppQ3/(2*baa*baa) + bVQ3)/3.
Omega_scfQ3 = rho_scfQ3/baCritQ3
##
bVQ4 = backgroundQ4['V_scf']
bppQ4 = backgroundQ4["phi'_scf"]
baCritQ4 = backgroundQ4['(.)rho_crit']
rho_scfQ4 = (bppQ4*bppQ4/(2*baa*baa) + bVQ4)/3.
Omega_scfQ4 = rho_scfQ4/baCritQ4


colours = ['g']
#for name in namelist:
#    idx = namelist.index(name)
#    plt.loglog(baLCDM['a'],fLCDM*baLCDM[name],colours[idx]+'-')
#plt.legend(namelist,loc='upper left')
#for name in namelist:
#    idx = namelist.index(name)
#plt.xlim([0.0, 10.])
#plt.ylim([0.0, 1.0])
plt.semilogx(1/(1+background['z']), baCC/baCritQ,color='r', label='LCDM')
plt.semilogx(1/(1+backgroundQ['z']),rho_scfQ/baCritQ, color='b', label='$\lambda=1.4, A=10$')
plt.semilogx(1/(1+backgroundQ2['z']),rho_scfQ2/baCritQ2,color='g', label='$\lambda=1.4, A=40$')
plt.semilogx(1/(1+backgroundQ3['z']),rho_scfQ3/baCritQ3,color='c', label='$\lambda=1.4, \\alpha=0$')
plt.semilogx(1/(1+backgroundQ4['z']),rho_scfQ4/baCritQ4,color='m', label='$\lambda=1.4, \\alpha=1$')

plt.xlim([1e-30, 1.1])
#plt.ylim([0.0, 20])

plt.xlabel(r"$a$")
plt.ylabel(r"$\mathrm{\Omega_{scf}}$")
plt.tight_layout()
plt.legend()

# In[ ]:


plt.savefig('Omega_scf.pdf')

