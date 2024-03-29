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
                   'Omega_cdm':0.264,
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
Mscf={}
MSFDM={}
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
    #
M = Class()
M.set(common_settings)
MQ=Class()
MQ.set(common_settings)
Mscf=Class()
Mscf.set(common_settings)
MSFDM=Class()
MSFDM.set(common_settings)
#####
Mscf.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004,22.656, 14.2, 0.0',
                   'scf_tuning_index':0,})

MSFDM.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-22., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2})
                  
MQ.set({'a_ini_over_a_today_default':1.e-14, 'Omega_cdm':0.0001,'Omega_sfdm_1':0.264,'Omega_sfdm_2':0.0,'attractor_ic_sfdm_1': 'yes',
                   'sfdm_parameters_1': '-22., 1.e4, 1.e-2, 1.e-16, 1.e-30',
                   'sfdm_tuning_index_1':2,'Omega_Lambda':0,'Omega_scf':-0.1,'Omega_fld':0.,'Omega_Lambda':0,'attractor_ic_scf': 'no',
                   'scf_parameters': '11.96, 2.0, 0.004, 22.656, 14.2, 0',
                   'scf_tuning_index':0,})

M.compute()
MQ.compute()
Mscf.compute()
MSFDM.compute()
    #load perturbations
all_k=M.get_perturbations()
one_k=all_k['scalar']
background = M.get_background()
    ##
all_kQ=MQ.get_perturbations()
one_kQ=all_kQ['scalar']
backgroundQ = MQ.get_background()
#
all_kscf=Mscf.get_perturbations()
one_kscf=all_kscf['scalar']
background_scf = Mscf.get_background()
#
all_kSFDM=MSFDM.get_perturbations()
one_kSFDM=all_kSFDM['scalar']
background_SFDM = MSFDM.get_background()
#

#Defining an array of k-values
kvec = np.logspace(-4,np.log10(3),1000)

pkM = []
pkMQ= []
pkMscf = []
pkMSFDM = []

for k in kvec:
    pkM.append(M.pk(k,0.))
    pkMQ.append(MQ.pk(k,0.))
    pkMscf.append(Mscf.pk(k,0.))
    pkMSFDM.append(MSFDM.pk(k,0.))
    h = M.h()
# plotting
#################
#
plt.loglog(kvec/h,np.array(pkM)/np.array(pkM),'r',linestyle='-', label='$\Lambda$CDM') 

plt.loglog(kvec/h, np.array(pkMscf)/np.array(pkM),'m',linestyle='-',label='$\psi$') 

#plt.loglog(kvec/h, np.array(pkMSFDM)/np.array(pkM),'c',linestyle='-', label='$\phi$') 

#plt.loglog(kvec/h, np.array(pkMQ)/np.array(pkM),'k',linestyle='-', label='$\psi+\phi$') 

#
#plt.legend(title='Albrecht Skordis')
plt.xlim([0.00012,10])
plt.legend()
plt.ylim([0.717,1.08])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$')
plt.ylabel(r'$P^{\phi, \psi}/P_{\Lambda\rm CDM}$')
#plt.tick_params(which='minor',axis='both',direction='in',right=True,top=True,
#                    length=4,width=1)
#plt.tick_params(which='major',axis='both',direction='in',right=True,top=True,
#                    length=7,width=1.5)
plt.savefig('scripts/Plots/MPS_AS_sfdm2.pdf')
    
