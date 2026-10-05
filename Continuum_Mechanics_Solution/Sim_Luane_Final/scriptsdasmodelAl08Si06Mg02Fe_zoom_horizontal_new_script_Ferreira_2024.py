from numpy import trapezoid as trapz, abs, zeros, sqrt, exp, log, pi, cos, sin
import operator as op
import numpy as np
import matplotlib.pyplot as plt
from math  import exp, gamma, log
from scipy.optimize import brentq
import scipy.optimize as so
from midpoint import midpoint

# Alloy Al-0.8Si-0.6Mg-0.20Fe Horizontal

global Error

vc01 = np.zeros(200)
vk01 = np.zeros(200)
vbeta1 = np.zeros(200)
vgama1 = np.zeros(200)
vc02 = np.zeros(200)
vk02 = np.zeros(200)
vbeta2 = np.zeros(200)
vgama2 = np.zeros(200)
vc03 = np.zeros(200)
vk03 = np.zeros(200)
vbeta3 = np.zeros(200)
vgama3 = np.zeros(200)

vFo1 = np.zeros(200)
vFo2 = np.zeros(200)
vFo3 = np.zeros(200)
vtSL = np.zeros(200)


def Diffusion_coefficient(T,D0,Qd,R):
 
    Diff = D0 * exp(-Qd/(R*T))

    return Diff 


xsc = np.zeros(20)
ysc = np.zeros(20)
yscp = np.zeros(20)
yscm = np.zeros(20)
error_ysc = np.zeros(20)

# SDAS Experimental Data


xsc[0] = 4.9592 ##2.1404022424680296
ysc[0] = 11.6712 #7.22470895906067
xsc[1] = 10.19725  ##5.495451957722622
ysc[1] = 13.69645 #9.608248538333466
xsc[2] = 15.54597 ##9.434625617168729
ysc[2] = 15.69645 #11.271402893473208
xsc[3] = 19.87864 ##14.10853323919925
ysc[3] = 18.3601  #13.22502247555545
xsc[4] = 31.96606 ##24.21721922429099
ysc[4] = 20.626 #17.00352690913498
xsc[5] = 54.37656 ##36.220167953428565
ysc[5] = 25.67347  #18.413006868615454
xsc[6] = 77.1587  ##48.288756228643976
ysc[6] = 30.89785  #19.049192733084187
xsc[7] = 100.20633 ##63.62622510989415
ysc[7] = 32.94919 #22.354896434181878
#xsc[8] = #77.36522989642918
#ysc[8] = 23.396834092669057

yscp[0] = 2.3528
yscp[1] = 5.40855
yscp[2] = 6.10955
yscp[3] = 5.7559
yscp[4] = 7.891
yscp[5] = 7.05053
yscp[6] = 5.04215
yscp[7] = 6.28781


yscm[0] = 3.4562
yscm[1] = 4.21245
yscm[2] = 3.99045
yscm[3] = 5.3821
yscm[4] = 5.078
yscm[5] = 5.36347
yscm[6] = 8.15885
yscm[7] = 5.23319

"""

xsc[0] = 48.8 
ysc[0] = 31.718275 
xsc[1] = 85.4 
ysc[1] = 33.6492615 
xsc[2] = 128 
ysc[2] = 36.525645 
xsc[3] = 208
ysc[3] = 42.912604  
xsc[4] = 269 
ysc[4] = 50.38717 
xsc[5] = 328
ysc[5] = 54.3  




yscp[0] = 5.918675
yscp[1] = 2.94236
yscp[2] = 1.96885
yscp[3] = 2.9276
yscp[4] = 7.37857
yscp[5] = 7.04


yscm[0] = 4.65143
yscm[1] = 1.89694
yscm[2] = 1.77025
yscm[3] = 2.6122
yscm[4] = 8.88003
yscm[5] = 4.92

"""

for i in range(12):
    error_ysc[i] = 0.1 * ysc[i]

R = 8.31451

# Input data
TL = 925.32 #923.86 #866.5 #933.15   # Verify after
T_Al = TL
TF = 933.15
Ts = 525.0 + 273.15
DeltaH = 335300 #260000.0 
rhos = 2555.72 #2556.0
rhol = 2378.33 #2386.0
cps  = 1270.0
sigmasl = 0.914
gamma_0 = 0.154
mu0 = 1.6
lambda0 = -5.0
rho0 = 2.56e-6
mu = 2.594e10
lambda_ = 5.034e10
sigma = sigmasl



Delta_Sv_hom = rhos * DeltaH/TL


# Data solute Mg
Dol_Cu = 9.9e-5 #1.06e-7
Qol_Cu = 71600.0 #24000
Dos_Cu = 6.23e-6 #6.5e-5
Qos_Cu = 115000.0 #136100
c01=0.6
cf1= 34.2 #28.816 
DL1= Diffusion_coefficient(TL,Dol_Cu,Qol_Cu,R) 
DS1= Diffusion_coefficient(Ts,Dos_Cu,Qos_Cu,R) 
k01= 0.36 #0.1308  
mL1= -5.1 #-3.64103  


# Data solute Si
Dol_Si = 1.34E-07
Qol_Si = 30000.0 
Dos_Si = 2.48E-04
Qos_Si = 137000.0
c02= 0.8 #4.0
cf2= 12.7 #10.7275 
DL2= Diffusion_coefficient(TL,Dol_Si,Qol_Si,R) 
DS2= Diffusion_coefficient(Ts,Dos_Si,Qos_Si,R) 
k02= 0.11 #0.1075
mL2= -6.2 #-7.198322069


# Data solute Fe
Dol_Fe = 2.34E-07
Qol_Fe = 35000.0
Dos_Fe = 5.30E-03
Qos_Fe = 183400.0
c03= 0.2  #0.11
cf3= 1.81 #1.655
DL3= Diffusion_coefficient(TL,Dol_Fe,Qol_Fe,R) 
DS3= Diffusion_coefficient(Ts,Dos_Fe,Qos_Fe,R) 
k03= 0.03 #0.0162 
mL3= -4.2 #-2.888217523
  


# Nucleation Calculation data





r_hom = 3.9821e-06  #4.84982e-06              #4.86322e-06              # [m]
r_het = 1.93141e-06 #1.87165e-06              #1.86834e-06              # [m]
thet =  0.952871 # 0.783201                 #0.780197                 # [rad]
dTdr = 7902.38 #4267.76                   #4205                     # [K.m^-1]
DT_= 0.79083 #0.52019                     #0.51692                  # [K]
sigma = -2.09008 #-1.66031                 #-1.65425                 # [N.m^-1]
surface_stress = -1.00008 #-0.570309       #-0.564251                # [N.m^-1]
gamsl = 1.45814 #0.879113                 #0.871318                 # [J.m^-2]
dgamsldr = -777538 #-576854               #-574516                  # [J.m^-3]
Gibbs_Thomson_hom = 7.53446e-07 #4.86812e-07  #4.80167e-07              # [m.K]
Gibbs_Thomson_het_1st = 7.63702e-07 #4.86809e-07      #4.80e-07         # [m.K]
thet_2nd = 3.04624 #3.04699               #3.04703                  # [rad]
Gibbs_Thomson_het =  1.57452e-06 #1.26141e-06 #1.25023e-06 ### 2nd-Order ### 3.41374e-07 # [m.K]      # 2nd Homogeneous --> 2.208e-07 
Gibbs_Thomson_non_Eq = Gibbs_Thomson_het
Gibbs_Thomson = Gibbs_Thomson_het
GIBBSTHOMSON = Gibbs_Thomson_het 


# start guess on the solvent concentration in solid alpha phase

csv = 100 - k01 * c01 - k02 * c02 - k03 * c03

# Experimental solidification kinetics
P  = np.arange(1,101, 1)
VL = 1.2*P**(-0.27)     # VL Horizontal
#V  = VL / 1000.0
tSL= 3.58 * P ** (1.09)       # tSL Horizontal 

# Experimental solidification kinetics
P_non_Eq  = np.arange(1,101, 1.0)
VL_non_Eq = 5.0*P**(-0.43)
#V  = VL / 1000.0
tSL_non_Eq= 0.32*P**(1.36)



def diff_length_scale(P, VL, tSL):
   V  = VL / 1000.0
   a = min(V)
   b = max(V)
   # Model Equations
   # Liquid diffusion length scale
   delta1 = DL1/V
   delta2 = DL2/V
   delta3 = DL3/V
   # Mean liquid boundary layer affected diffusion
   a1 = trapz(V,delta1)
   a2 = trapz(V,delta2)
   a3 = trapz(V,delta3)

   minV = min(V)
   maxV = max(V)

   delta1m=abs((1.0/(maxV-minV))*a1)
   delta2m=abs((1.0/(maxV-minV))*a2)
   delta3m=abs((1.0/(maxV-minV))*a3)


   # Effective partition coefficient
   kef1= k01/(k01+(1.0 - k01)*(np.exp((-V*delta1m/DL1))))
   kef2= k02/(k02+(1.0 - k02)*(np.exp((-V*delta2m/DL2))))
   kef3= k03/(k03+(1.0 - k03)*(np.exp((-V*delta3m/DL3))))


   # Solute solid fractions
   w1=(c01/100)/((c01/100)+(c02/100)+(c03/100))
   w2=(c02/100)/((c01/100)+(c02/100)+(c03/100))
   w3=(c03/100)/((c01/100)+(c02/100)+(c03/100))


   return delta1m, delta2m, delta3m, kef1, kef2, kef3, w1, w2, w3


def func_f(f, c01, k01, beta1, c02, k02, beta2, c03, k03, beta3, csv):
    csv = -c01 * k01 * (1 - f * (1 - beta1 * k01)) ** ((k01-1)/(1-beta1 * k01)) - c02 * k02 * (1 - f * (1 - beta2 * k02)) ** ((k02-1)/(1-beta2 * k02)) - c03 * k03 * (1 - f * (1 - beta3 * k03)) ** ((k03-1)/(1-beta3 * k03)) + 100
#    print(f'function_f({f}) = {function_f}')
    return csv

def cs(f, c0, k0, beta):
    cs = c0 * k0 * (1 - f * (1 - beta * k0)) ** ((k0-1)/(1-beta * k0))
    return cs

def T(f, TF, TL, k0, beta):
    T = TF + (TL-TF)* (1 - f * (1 - beta * k0)) ** ((k0-1)/(1-beta * k0))
    return T


def dfdtau(f, Fo, gama, beta):
    df_dtau = Fo * (gama - beta)/(f * gama * beta)
    return df_dtau

#def func_DH(f, rhos, cps, TL, Tf, beta, k):
    funcDH = ( rhos * cps * (TF - TL) * (1 - k) * (1 - f * (1 - beta * k)) ** ((k - 1)/(1 - beta * k)) )/((1 - beta * k)*(1 - f*(1- beta * k)))
    return funcDH
 
def dcs_approx(f, c0, k0, beta):
    dcs_app = c0*k0*(1 - k0)*(-f*(-beta*k0 + 1) + 1)**((k0 - 1)/(-beta*k0 + 1))/((-f*(-beta*k0 + 1) + 1))
    return dcs_app

def dcs(f, c0, k0, beta, dbeta):
    dcs_nonap = c0*k0*(-f*(-k0*beta + 1) + 1)**((k0 - 1)/(-k0*beta + 1)) * (k0*(k0 - 1)*log(-f*(-k0*beta + 1) + 1)*dbeta/(-k0*beta + 1)**2 + (k0 - 1)*(f*k0*dbeta + k0*beta - 1)/((-f*(-k0*beta + 1) + 1)*(-k0*beta + 1)))
    return dcs_nonap

def dbetadf(f, dfdtau, Fo, gama, beta):
    dbeta = -Fo * gama**2 * dfdtau / ((Fo + gama * f * dfdtau)**2)
    return dbeta


SDAS2 = 1.0
beta_CU = 0.0
beta_SI = 0.0
Error=1.0
def calc(G_T, delta1m, delta2m, delta3m, kef1, kef2, kef3, w1, w2, w3):
  GIBBSTHOMSON = G_T 
  SDAS2 = 1.0
  Error=1.0
  while np.all(op.gt(Error, 1.e-10)):
    SDAS = SDAS2
    # Fourier number / SDAS initial guess
    Fo1=DS1*tSL/((SDAS/2)**2.0)
    Fo2=DS2*tSL/((SDAS/2)**2.0)
    Fo3=DS3*tSL/((SDAS/2)**2.0)
    # Macrosegregation  parameters Gama
    gama1=(4.0*Fo1*k01)/(4.0*Fo1*k01+1.0)
    gama2=(4.0*Fo2*k02)/(4.0*Fo2*k02+1.0)
    gama3=(4.0*Fo3*k03)/(4.0*Fo3*k03+1.0)
    # Segregation parameter betha
    beta1=(2.0*gama1*Fo1)/(2.0*Fo1+gama1)
    beta2=(2.0*gama2*Fo2)/(2.0*Fo2+gama2) 
    beta3=(2.0*gama3*Fo2)/(2.0*Fo3+gama3)  
    # Weighted betha
    bethaweighted=(beta1*w1)+(beta2*w2)+(beta3*w3)  
    # Rappaz and Boettinger Modified Model
    p1=(mL1*(1.0 - kef1)*(cf1-c01)/DL1)+(mL2*(1.0-kef2)*(cf2-c02)/DL2)+(mL3*(1.0-kef3)*(cf3-c03)/DL3)
    s2=(mL1*(1.0 - kef1)*( cf1)/ DL1)+(mL2*(1.0-kef2)*(cf2)/DL2)+(mL3*(1.0-kef3)*(cf3)/DL3)
    t3=(mL1*(1.0 - kef1)*( c01)/ DL1)+(mL2*(1.0-kef2)*(c02)/DL2)+(mL3*(1.0-kef3)*(c03)/DL3)
    # M e M^1/3
    mrbm=bethaweighted*((-GIBBSTHOMSON/p1)*np.log(s2/t3))
    mrbm2=5.5*(mrbm**(1/3))
    # Corrected SDAS value
    SDAS2= mrbm2*(tSL)**(1/3)
    Error = np.abs(SDAS-SDAS2)/SDAS
    
  return SDAS, beta1, beta2, beta3

  


def calc_RB(G_T, k01, k02, k03, c01, cf1, c02, cf2,c03,cf3, mL1, mL2, mL3, DL1, DL2, DL3, tSL):
  GIBBSTHOMSON = G_T
  # Rappaz and Boettinger
  p1=(mL1*(1.0 - k01)*(cf1-c01)/DL1)+(mL2*(1.0-k02)*(cf2-c02)/DL2)+(mL3*(1.0-k03)*(cf3-c03)/DL3)
  s2=(mL1*(1.0 - k01)*( cf1)/ DL1)+(mL2*(1.0-k02)*(cf2)/DL2)+(mL3*(1.0-k03)*(cf3)/DL3)
  t3=(mL1*(1.0 - k01)*( c01)/ DL1)+(mL2*(1.0-k02)*(c02)/DL2)+(mL3*(1.0-k03)*(c03)/DL3)
  # M e M^1/3
  mrbm=((-GIBBSTHOMSON/p1)*np.log(s2/t3))
  mrbm2=5.5*(mrbm**(1/3))
  # Calculate SDAS value
  SDAS= mrbm2*(tSL)**(1/3)
  return SDAS


SDAS_RB = calc_RB(Gibbs_Thomson_het, k01, k02, k03, c01, cf1, c02, cf2, c03, cf3, mL1, mL2, mL3, DL1, DL2,DL3, tSL)
print(f'Gibbs-Thomson = {Gibbs_Thomson}\n')
print(f'k01 = {k01}')
print(f'k02 = {k02}')
print(f'k03 = {k03}')
print(f'c01 = {c01}')
print(f'cf1 = {cf1}')
print(f'c02 = {c02}')
print(f'cf2 = {cf2}')
print(f'c03 = {c03}')
print(f'cf3 = {cf3}')
print(f'mL1 = {mL1}')
print(f'mL2 = {mL2}')
print(f'mL3 = {mL3}')
print(f'DL1 = {DL1}')
print(f'DL2 = {DL2}')
print(f'DL3 = {DL3}')
print(f'tSL = {tSL}')

print(f'SDAS_RB = {SDAS_RB} [K.m]')
 

delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq = diff_length_scale(P, VL, tSL)
print(f'delta1m_Eq = {delta1m_Eq}')
SDAS_Eq, beta_CU, beta_SI, beta_FE = calc(Gibbs_Thomson_het, delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)
print(f'SDAS_Eq = {SDAS_Eq} [K.m]')

SDAS_het_1st_Eq, beta_het_1st_CU, beta_het_1st_SI, beta_het_1st_FE = calc(Gibbs_Thomson_het_1st, delta1m_Eq, delta2m_Eq, delta3m_Eq, kef1_Eq, kef2_Eq, kef3_Eq, w1_Eq, w2_Eq, w3_Eq)
print(f'SDAS_Eq = {SDAS_Eq} [K.m]')



delta1m_non_Eq, delta2m_non_Eq, delta3m_non_Eq, kef1_non_Eq, kef2_non_Eq, kef3_non_Eq, w1_non_Eq, w2_non_Eq, w3_non_Eq = diff_length_scale(P_non_Eq, VL_non_Eq, tSL_non_Eq)
print(f'delta1m_non_Eq = {delta1m_non_Eq}')
SDAS_non_Eq, beta_CU_non_Eq, beta_SI_non_Eq, beta_FE_non_Eq = calc(Gibbs_Thomson_het, delta1m_non_Eq, delta2m_non_Eq, delta3m_non_Eq, kef1_non_Eq, kef2_non_Eq, kef3_non_Eq, w1_non_Eq, w2_non_Eq, w3_non_Eq)
print(f'SDAS_non_Eq = {SDAS_non_Eq} [K.m]')

print('### Este arquivo ####')

# plot data
lx_RB = np.log(tSL)
ly_RB = np.log(SDAS_RB*1e06)
z3_RB= np.polyfit(lx_RB,ly_RB,1)
z3n_RB= np.poly1d(z3_RB)


lx_Eq = np.log(tSL)
ly_Eq = np.log(SDAS_Eq*1e06)
z3_Eq= np.polyfit(lx_Eq,ly_Eq,1)
z3n_Eq= np.poly1d(z3_Eq)


lx_1st_Eq = np.log(tSL)
ly_1st_Eq = np.log(SDAS_het_1st_Eq*1e06)
z3_1st_Eq= np.polyfit(lx_1st_Eq,ly_1st_Eq,1)
z3n_1st_Eq= np.poly1d(z3_1st_Eq)


lx_non_Eq = np.log(tSL_non_Eq)
ly_non_Eq = np.log(SDAS_non_Eq*1e06)
z3_non_Eq= np.polyfit(lx_non_Eq,ly_non_Eq,1)
z3n_non_Eq= np.poly1d(z3_non_Eq)

fig1,ax1 = plt.subplots(1)
ax1.set_yscale('log')
ax1.set_xscale('log')
#ax1.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax1.set_title(r'SDAS Models Predictions against Experimental Scatter',{'color': 'black', 'fontsize': 14})
###ax1.plot(tSL,(SDAS_RB*1e06),'r>',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_RB[0]))*tSL**z3n_RB[1],'-k',linewidth = 1.5)
###ax1.plot(tSL,(SDAS_Eq*1e06),'g^',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_Eq[0]))*tSL**z3n_Eq[1],'-.k',linewidth = 1.5)
###ax1.plot(tSL,(SDAS_het_1st_Eq*1e06),'b^',linewidth = 1.5)
ax1.plot(tSL,(np.exp(z3n_1st_Eq[0]))*tSL**z3n_1st_Eq[1],'--k',linewidth = 1.5)


ax1.plot(xsc, ysc, 'sk',linewidth = 1.5)
#ax1.errorbar(xsc, ysc, uplims=yscp,linestyle=ls)



ax1.errorbar(xsc, ysc, yerr=yscp, fmt='.k')
ax1.errorbar(xsc, ysc, yerr=yscm, fmt='.k')

ax1.set_ylabel(u'Secondary Dendrite Arm Spacing, $SDAS$ [$\mathrm{{\mu}m}$]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax1.set_xlabel(r'Local solidification time, $\mathrm{t_{SL}}$ [s]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.text(15,380,'Alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe Horizontal', {'color': 'black', 'fontsize': 16})
###plt.text(1.5,210,r'$\mathit{SDAS_{RB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_RB[0]), z3n_RB[1]), {'color': 'red', 'fontsize': 12})
###plt.text(1.5,120,r'$\mathit{SDAS_{MBR}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_1st_Eq[0]), z3n_1st_Eq[1]), {'color': 'green', 'fontsize': 12})
###plt.text(1.5,70,r'$\mathit{SDAS_{MRB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_Eq[0]), z3n_Eq[1]), {'color': 'blue', 'fontsize': 12})
#plt.text(1.5,2000,r'$\mathit{V_{L}\,=\,5.0\,{P}^{-0.41} \,\, \mathrm{[mm.s^{-1}]}}$', {'color': 'black', 'fontsize': 12})
#plt.text(1.5,1100,r'$\mathit{t_{SL}\,=\,0.0036\,{P}^{\, 2.32} \,\, \mathrm{[s]}}$', {'color': 'black', 'fontsize': 12})
#plt.text(1.5,500,r'$\mathit{\overline{\Gamma } _{het} ^{\,1st-Order} \,=\, %g \,\, \mathrm{[m.K]}}$' % (Gibbs_Thomson_het_1st), {'color': 'green', 'fontsize': 12})
#plt.text(1.5,250,r'$\mathit{\theta  \,=\, %g \,\, \mathrm{[rad]}}$' % (thet), {'color': 'black', 'fontsize': 12})
#plt.text(1.5,130,r'$\mathit{\overline{\Gamma } ^{2nd-Order}_{het} \,=\, %g \,\, \mathrm{[m.K]}}$' % (Gibbs_Thomson_het), {'color': 'blue', 'fontsize': 12, 'weight':'bold' })
#plt.text(1.5,60,r'$\mathit{\theta ^{2nd-Order}  \,=\, %g \,\, \mathrm{[rad]}}$' % (thet_2nd), {'color': 'blue', 'fontsize': 12})
#plt.text(20.0,12000,r'$\mathit{\overline{\nabla T}_{EXP} \,= %g \, \mathrm{[K.m^{-1}]}}$' % (dTdr), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,7000,r'$\mathit{r_{C,hom}\,= %g \, \mathrm{[m]}}$' % (r_hom), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,4000,r'$\mathit{r_{C,het}\,= %g \, [m]}$' % (r_het), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,2000,r'$\mathit{\sigma _{SL}\,= %g \, \mathrm{[N.m^{-1}]}}$' % (sigma), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,1100,r'$\mathit{s \, = %g \, [N.m^{-1}]}$' % (surface_stress), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,500,r'$\mathit{\gamma _{SL}\,= %g \, \mathrm{[J.m^{-2}]}}$' % (gamsl), {'color': 'black', 'fontsize': 12})
#plt.text(20.0,250,r'$\mathit{\frac {\partial \gamma _{SL}} {\partial r} \,= %g \, \mathrm{[J.m^{-3}]}}$' % (dgamsldr), {'color': 'black', 'fontsize': 12})
#plt.text(150,12000,r'$\mathit{\sigma _{0}\,= %g \, \mathrm{[N.m^{-1}]}}$' % (-sigmasl), {'color': 'black', 'fontsize': 12})
#plt.text(150,7000,r'$\mathit{\gamma _{0}\,= %g \, \mathrm{[J.m^{-2}]}}$' % (gamma_0), {'color': 'black', 'fontsize': 12})
#plt.text(1.5,20000,r'$^*$Thermal gradient $\mathrm{\nabla T}$ for horizontal experiments were only measured at the symmetry line along the X axis.' , {'color': 'black','weight':'normal', 'fontsize': 10})



plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'RB SDAS Model (Rappaz and Boettinger, 1999), $\mathit{SDAS_{RB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_RB[0]), z3n_RB[1]),r'MRB Model (Ferreira, 2024) with $\Gamma ^{1nd-Order} _{Het}$ (Present work), $\mathit{SDAS_{MBR}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_1st_Eq[0]), z3n_1st_Eq[1]),r'MRB Model (Ferreira, 2024) with $\Gamma ^{2nd-Order} _{Het}$ (Present work), $\mathit{SDAS_{MRB}\,=\,%.2f\,t_{SL}^{%.3f} \,\, \mathrm{[\mu m]}}$' % (np.exp(z3n_Eq[0]), z3n_Eq[1]), r'Experimental (Marques et al., 2024)'),frameon=True,edgecolor='k',fontsize=10,shadow=True, loc= 'lower right')
###plt.legend((r'RB SDAS Model (Rappaz and Boettinger, 1999)','RB Model Curve Fit',r'MRB Model (Ferreira et al.,2018) with $\Gamma ^{1nd-Order} _{Het}$ (Present work)',r'MRB Model with $\Gamma ^{1nd-Order} _{Het}$ Curve Fit ',r'MRB Model (Ferreira et al.,2018) with $\Gamma ^{2nd-Order} _{Het}$ (Present work)',r'MRB Model with $\Gamma ^{2nd-Order} _{Het}$ Curve Fit ', r'Experimental (Costa et al., 2015)'),frameon=True,edgecolor='k',fontsize=10,shadow=True, loc= 'lower right')
#plt.legend(('SDAS Predictions Equilibrium, Rappaz and Boettinger (1999)','SDAS Fit, Rappaz and Boettinger (1999)', r'MRB, Ferreira et al. (2018) and Ferreira (2022)' ,r'SDAS Fit , Ferreira et al. (2018) and Ferreira (2022)','Experimental'),frameon=True,edgecolor='k',shadow=True, loc= 'lower right')
#plt.axis([0,0.35,-0.2e7,0])
plt.axis([1e0,1e3,1e0,1e3])


fig2,ax2 = plt.subplots(1)
ax2.set_xscale('log')
#ax2.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax2.set_title(r'Horizontal Solidification Back-Diffusion Parameter, $\beta _i$',{'color': 'black', 'fontsize': 14})
ax2.plot(tSL,beta_CU,'k-',linewidth = 1.5)
ax2.plot(tSL,beta_SI,'k-.',linewidth = 1.5)
ax2.plot(tSL,beta_FE,'k--',linewidth = 1.5)
plt.text(5,0.080,'Alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe Horizontal', {'color': 'black', 'fontsize': 16})
#plt.text(5,0.150,r'$\mathit{T_{L}\,= %g \, [K]}$' % (TL), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.137,r'$\mathit{T_{S}\,= %g \, [K]}$' % (Ts), {'color': 'black', 'fontsize': 14})
#plt.text(5,.125,r'$\mathit{D _{L,Si}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Si,Qol_Si), {'color': 'black', 'fontsize': 14})
#plt.text(5,.1125,r'$\mathit{D _{S,Si}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Si,Qos_Si), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.100,r'$\mathit{D _{L,Cu}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Cu,Qol_Cu), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.087 ,r'$\mathit{D _{S,Cu}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Cu,Qos_Cu), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.075,r'$\mathit{D _{L,Fe}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dol_Fe,Qol_Fe), {'color': 'black', 'fontsize': 14})
#plt.text(5,0.063 ,r'$\mathit{D _{S,Fe}\, = %.3e  e^{\frac {- %g} {R T}} \, \mathrm{[m^{2}.s^{-1}]}}$' % (Dos_Fe,Qos_Fe), {'color': 'black', 'fontsize': 14})


ax2.set_xlabel(r'Local solidification time, $t_{SL}$ [s]', {'color': 'black','fontweight':'normal', 'fontsize': 17})
ax2.set_ylabel(r'Back-Diffusion Parameter, $\beta _i$  ', {'color': 'black','fontweight':'normal', 'fontsize': 17})
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'$\beta _{Mg}$',r'$\beta _{Si}$', r'$\beta _{Fe}$'),frameon=True,edgecolor='k',shadow=True, loc= 'upper right',fontsize= 17)
plt.axis([0,1e3,0,0.1])




plt.show()

