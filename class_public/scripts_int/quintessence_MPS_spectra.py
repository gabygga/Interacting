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
MQ.set({'Omega_Lambda':1e-5,'Omega_scf':-0.1,'Omega_fld':0.,'attractor_ic_scf': 'no',
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
                   'scf_parameters': '-1.4, 0.0, 10, 0.0, 20, 0.0',
                   'scf_tuning_index':0,})

M.compute()
MQ.compute()
MQ2.compute()
MQ3.compute()
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
   


#Defining an array of k-values
kvec = np.logspace(-4,np.log10(3),1000)

pkM = []
pkMQ= []
pkMQ2 = []
pkMQ3 = []

for k in kvec:
    pkM.append(M.pk(k,0.))
    pkMQ.append(MQ.pk(k,0.))
    pkMQ2.append(MQ2.pk(k,0.))
    pkMQ3.append(MQ3.pk(k,0.))
h = M.h()
# plotting
#################
#
plt.semilogx(kvec/h,np.array(pkM)/np.array(pkM),'r',linestyle='-', label='$\Lambda$CDM', lw='2.5') 

plt.semilogx(kvec/h,(pkM[0]/pkMQ[0])*np.array(pkMQ)/np.array(pkM),'b',linestyle='-', label=r'$A=10, \lambda = 1.4$', lw='2.5') 
plt.semilogx(kvec/h,(pkM[0]/pkMQ2[0])*np.array(pkMQ2)/np.array(pkM),'g',linestyle='-',label=r'$A=40, \lambda = 1.4$')
#plt.semilogx(kvec/h,(pkM[0]/pkMQ3[0])*np.array(pkMQ3)/np.array(pkM),'c',linestyle='-', label='$g_{\\rm eff}=1x10^{-9}, \Sigma m_{\\nu}=0.06\\rm{eV}$', lw='2.5')

#
plt.legend(title='Quintessence: $V(\phi)=(1+A)\exp ^{-\lambda \phi}$')
plt.xlim([0.00012,3])
#plt.ylim([0.88,1.08])
plt.xlabel(r'$k \,\,\,\, [h/\mathrm{Mpc}]$',fontsize='16')
plt.ylabel(r'$P_{\phi}(k)/P_{\Lambda CDM}(k)$',fontsize='16')
plt.tick_params(which='minor',axis='both',direction='in',right=True,top=True,
                    length=4,width=1)
plt.tick_params(which='major',axis='both',direction='in',right=True,top=True,
                    length=7,width=1.5)
plt.savefig('MPS_quintessence.pdf')
    
