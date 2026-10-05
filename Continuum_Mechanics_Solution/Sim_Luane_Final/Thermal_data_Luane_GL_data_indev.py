#!/bin/bash
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, exp, log, pi, cos, sin
from math import pi
import csv
from scipy import integrate
from scipy.optimize import fsolve, brentq
from derivative_shim import derivative
#


"""
Prof. Ivaldo Leão Ferreira
Universidade Federal do Pará - UFPa
Instituto de Tecnologia - ITEC
Faculdade de Engenharia Mecânica - FEM
Av. Augusto Correa, 1
Belém - PA
Brazil
CEP 66075-110
ileao@ufpa.br
"""

n = 100

#P,t,VL,TR,GL
t_v   = zeros(n)
POS   = zeros(n)
POS_m = zeros(n)
VL    = zeros(n)
VL_m  = zeros(n)
TR    = zeros(n)
GL    = zeros(n)
#Experimental 
POS_EXP  = zeros(20)
TIME_EXP = zeros(20)
VL_EXP   = zeros(20)
POS1_EXP = zeros(20)
VEL_EXP = zeros(20)
POS2_EXP = zeros(20)
TR_EXP = zeros(20)
POS3_EXP = zeros(20)
GL_EXP = zeros(20)




# Experimental Phase Change
# Position versus time
TIME_EXP[1] = 3.116
POS_EXP[1] = 4.664 
TIME_EXP[2] = 10.198
POS_EXP[2] = 9.717 
TIME_EXP[3] = 23.796
POS_EXP[3] = 19.823 
TIME_EXP[4] = 36.827
POS_EXP[4] = 27.208  
TIME_EXP[5] = 60.057
POS_EXP[5] = 40.035 
TIME_EXP[6] = 86.119
POS_EXP[6] = 54.806 
TIME_EXP[7] = 120.68
POS_EXP[7] = 75.018 
TIME_EXP[8] = 153.824
POS_EXP[8] = 100.283


POS1_EXP[1] = 4.207
VEL_EXP[1] = 0.734
POS1_EXP[2] = 9.972
VEL_EXP[2] = 0.691
POS1_EXP[3] = 19.632
VEL_EXP[3] = 0.649
POS1_EXP[4] = 26.799
VEL_EXP[4] = 0.635
POS1_EXP[5] = 39.887
VEL_EXP[5] = 0.615
POS1_EXP[6] = 54.844
VEL_EXP[6] = 0.596
POS1_EXP[7] = 75.099
VEL_EXP[7] = 0.581
POS1_EXP[8] = 100.807
VEL_EXP[8] = 0.566


POS2_EXP[1] = 4.207
TR_EXP[1] = 29.506
POS2_EXP[2] = 9.66
TR_EXP[2] = 8.025
POS2_EXP[3] = 19.632
TR_EXP[3] = 4.691
POS2_EXP[4] = 26.799
TR_EXP[4] = 4.074
POS2_EXP[5] = 40.198
TR_EXP[5] = 2.716
POS2_EXP[6] = 55.156
TR_EXP[6] = 1.111
POS2_EXP[7] = 75.411
TR_EXP[7] = 0.864
POS2_EXP[8] = 99.873
TR_EXP[8] = 0.617

resTR_EXP = integrate.simpson(TR_EXP, x=POS2_EXP, dx=1.0, axis=-1)
print(f'integral = {resTR_EXP:g}')
TR_EXP_mean = 1/(POS2_EXP[8]-POS2_EXP[1]) * resTR_EXP
print(f'TR_EXP_mean = {TR_EXP_mean:g}[K.s^-1]')

scale = 1000.

POS3_EXP[1] = 4.207
GL_EXP[1] = (TR_EXP[1]/VEL_EXP[1])*scale
POS3_EXP[2] = 9.66
GL_EXP[2] = (TR_EXP[2]/VEL_EXP[2])*scale
POS3_EXP[3] = 19.632
GL_EXP[3] = (TR_EXP[3]/VEL_EXP[3])*scale
POS3_EXP[4] = 26.799
GL_EXP[4] = (TR_EXP[4]/VEL_EXP[4])*scale
POS3_EXP[5] = 40.198
GL_EXP[5] = (TR_EXP[5]/VEL_EXP[5])*scale
POS3_EXP[6] = 55.156
GL_EXP[6] = (TR_EXP[6]/VEL_EXP[6])*scale
POS3_EXP[7] = 75.411
GL_EXP[7] = (TR_EXP[7]/VEL_EXP[7])*scale
POS3_EXP[8] = 99.873
GL_EXP[8] = (TR_EXP[8]/VEL_EXP[8])*scale



#print(f'{POS3_EXP[1]:g},{GL_EXP[1]:g}')


# Position vs time
# Lscale milimeter to meters scale factor
Lscale = 1000.0    


aP = 0.92
bP = 0.92
def Pft(t,aP,bP):
    Pt = aP*t**(bP)
    return Pt

# Liquidus Velocity vs position
aVL = 0.84
bVL = -0.086
def VLfp(P,aVL, bVL):
    VLp = aVL*exp(bVL*log(P))
    return VLp

# Liquidus cooling rate
aTR = 305.7
bTR = -1.46 
def TRfp(P,aTR, bTR):
    TRp = aTR*exp(bTR*log(P))
    return TRp



t_total = 200

dt_ = t_total / n


for i in range(1,n):
    t_ = dt_ * i 
    t_v[i] = t_
    pos = Pft(t_,aP,bP)
    POS[i] = pos
    posm = pos/Lscale
    POS_m[i] = posm
    VL_  = VLfp(pos,aVL,bVL)
    VL[i] = VL_
    VLm = VL_/Lscale
    VL_m[i] = VLm
    TR_  = TRfp(pos,aTR,bTR)
    TR[i] = TR_
    GL_  = TR_ / VLm 
    GL[i] = GL_
    print(f'{t_:g}s  ----  {(posm):g}[m] --- {(VLm):g}[m/s] --- {TR_:g}[K/s] --- {GL_:g}[K/m]')

res = integrate.simpson(GL, x=POS, dx=1.0, axis=-1)
print(f'integral = {res:g}')
GL_mean = 1/(110-0) * res 
print(f'GL_mean = {GL_mean:g}[K.m^-1]')





fig1,ax1 = plt.subplots(1)
ax1.set_title(r'Liquidus Position $P$ function of $t$  $\mathrm{S_{L} = f(t)}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe   ',{'color': 'black', 'fontsize': 14})
ax1.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax1.plot(t_v[1:n],POS[1:n],'k-',linewidth = 1.5)
ax1.plot(TIME_EXP[1:8],POS_EXP[1:8],'ks',linewidth = 1.5)
ax1.set_ylabel(r'Liquidus Isothermal Position, $\mathit{P} \,\,\, \mathrm{ [mm]}$',{'color': 'black', 'fontsize': 18})
ax1.set_xlabel(r'Time, $\,\mathit{t} \,\, \mathrm{ [s]}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend(('Experimental Curve Fit','Experimental' ),frameon=True,edgecolor='k',shadow=True, loc= 'best', fontsize = 16 )
plt.axis([0,200,0,120])

fig2,ax2 = plt.subplots(1)
ax2.set_title(r'Liquidus Growth Rate $V_{L}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe  ',{'color': 'black', 'fontsize': 14})
ax2.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax2.plot(POS[1:n],VL[1:n],'k-',linewidth = 1.5)
ax2.plot(POS1_EXP[1:8],VEL_EXP[1:8],'ks',linewidth = 1.5)
#ax2.plot(TIME_EXP[0:8],POS_EXP[0:8],'ks',linewidth = 1.5)
ax2.set_ylabel(r'Liquidus Growth Rate, $\mathit{V_{L}} \,\,\, \mathrm{ [mm.s^{-1}]}$',{'color': 'black', 'fontsize': 18})
ax2.set_xlabel(r'Position, $\,\mathit{P} \,\, \mathrm{ [mm]}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend(('Experimental Curve Fit','Experimental' ),frameon=True,edgecolor='k',shadow=True, loc= 'best', fontsize = 16 )
plt.axis([0,120,0.5,0.8])

fig3,ax3 = plt.subplots(1)
ax3.set_title(r'Liquidus Cooling Rate $T_{R}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe  ',{'color': 'black', 'fontsize': 14})
ax3.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax3.plot(POS[1:n],TR[1:n],'k-',linewidth = 1.5)
ax3.plot(POS2_EXP[1:8],TR_EXP[1:8],'ks',linewidth = 1.5)
ax3.axhline(TR_EXP_mean,0,120,color='r', linestyle='--',lw=1.5)
#ax3.plot(TIME_EXP[0:8],POS_EXP[0:8],'ks',linewidth = 1.5)
ax3.set_ylabel(r'Liquidus Cooling Rate, $\mathit{T_{R}} \,\,\, \mathrm{ [K.s^{-1}]}$',{'color': 'black', 'fontsize': 18})
ax3.set_xlabel(r'Position, $\,\mathit{P} \,\, \mathrm{ [mm]}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax3.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax3.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Experimental Curve Fit','Experimental','Experimental Mean Integral $\overline{T}_R$ = %g $\mathrm{[K.s^{-1}]}$' % (TR_EXP_mean) ),frameon=True,edgecolor='k',shadow=True, loc= 'best', fontsize = 16 )
plt.axis([0,120,0,30])

fig4,ax4 = plt.subplots(1)
ax4.set_title(r'Liquidus Thermal Gradient $G_{L}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe  ',{'color': 'black', 'fontsize': 14})
ax4.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax4.plot(POS[1:n],GL[1:n],'k-',linewidth = 1.5)
#ax4.plot(POS[1:n],GL[1:n],'k-',linewidth = 1.5)
ax4.plot(POS3_EXP[1:8],GL_EXP[1:8],'ks',linewidth = 1.5)
ax4.axhline(GL_mean,0,10000,color='r', linestyle='--',lw=1.5)
#ax3.plot(TIME_EXP[0:8],POS_EXP[0:8],'ks',linewidth = 1.5)
ax4.set_ylabel(r'Liquidus Thermal Gradient, $\mathit{G_{L}} \,\,\, \mathrm{ [K.m^{-1}]}$',{'color': 'black', 'fontsize': 18})
ax4.set_xlabel(r'Position, $\,\mathit{P} \,\, \mathrm{ [mm]}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax4.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax4.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Experimental Curve Fit','Experimental','Experimental Mean Integral $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_mean) ),frameon=True,edgecolor='k',shadow=True, loc= 'best', fontsize = 16 )
plt.axis([0,120,0,120000])



plt.show()


