#!/bin/bash
from re import X
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, exp, log, pi, cos, sin, diff
from math import pi
import csv
from scipy import integrate
from scipy.optimize import fsolve, brentq, curve_fit
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

Tv    = zeros(n)
Cv     = zeros(n)


kvCu  = zeros(n)
CvCu  = zeros(n)

TvNb    = zeros(n)
CvNb     = zeros(n)


kvNb  = zeros(n)
CvNb  = zeros(n)



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
TIME_EXP[1] = 5.95596
POS_EXP[1] = 5 
TIME_EXP[2] = 9.70219
POS_EXP[2] = 10 
TIME_EXP[3] = 9.9052
POS_EXP[3] = 15 
TIME_EXP[4] = 18.13736
POS_EXP[4] = 20  
TIME_EXP[5] = 44.13588
POS_EXP[5] = 30
TIME_EXP[6] = 95.1389
POS_EXP[6] = 50
TIME_EXP[7] = 128.26357
POS_EXP[7] = 70 
TIME_EXP[8] = 230.1016
POS_EXP[8] = 90




POS1_EXP[1] = 5 #4.207
VEL_EXP[1] =  1.31773 ##0.734
POS1_EXP[2] = 10 #9.972
VEL_EXP[2] =  0.89396 ##0.691
POS1_EXP[3] = 15 #19.632
VEL_EXP[3] =  0.71244 ##0.649
POS1_EXP[4] = 20 #26.799
VEL_EXP[4] =  0.60647 ##0.635
POS1_EXP[5] = 30 #39.887
VEL_EXP[5] =  0.48333 ##0.615
POS1_EXP[6] = 50 #54.844
VEL_EXP[6] =  0.36312 ##0.596
POS1_EXP[7] = 70 #75.099
VEL_EXP[7] =  0.30079 ##0.581
POS1_EXP[8] = 90 #100.807
VEL_EXP[8] =  0.26131 ##0.566



POS2_EXP[1] = 5.0 #4.207
TR_EXP[1] = 11.81538 #29.506
POS2_EXP[2] = 10.0 #9.66
TR_EXP[2] = 7.58399 #8.025
POS2_EXP[3] = 14.0 #19.632
TR_EXP[3] = 7.67473 #4.691
POS2_EXP[4] = 19.0 #26.799
TR_EXP[4] = 4.65463 #4.074
POS2_EXP[5] = 29.0 # 40.198
TR_EXP[5] = 2.53968 #2.716
POS2_EXP[6] = 50.0 #55.156
TR_EXP[6] = 1.04679 # 1.111
POS2_EXP[7] = 69.0 #75.411
TR_EXP[7] = 0.33553 # 0.864
POS2_EXP[8] = 90.0 # 99.873
TR_EXP[8] = 0.29251 #0.617



resTR_EXP = integrate.simpson(TR_EXP, x=POS2_EXP, dx=1.0, axis=-1)
print(f'integral = {resTR_EXP:g}')
TR_EXP_mean = 1/(POS2_EXP[8]-POS2_EXP[1]) * resTR_EXP
print(f'TR_EXP_mean = {TR_EXP_mean:g}[K.s^-1]')

scale = 1000.

POS3_EXP[1] = 5.0
GL_EXP[1] = (TR_EXP[1]/VEL_EXP[1])*scale
POS3_EXP[2] = 10.0
GL_EXP[2] = (TR_EXP[2]/VEL_EXP[2])*scale
POS3_EXP[3] = 14.0
GL_EXP[3] = (TR_EXP[3]/VEL_EXP[3])*scale
POS3_EXP[4] = 19.0
GL_EXP[4] = (TR_EXP[4]/VEL_EXP[4])*scale
POS3_EXP[5] = 29.0
GL_EXP[5] = (TR_EXP[5]/VEL_EXP[5])*scale
POS3_EXP[6] = 50.0
GL_EXP[6] = (TR_EXP[6]/VEL_EXP[6])*scale
POS3_EXP[7] = 69.0
GL_EXP[7] = (TR_EXP[7]/VEL_EXP[7])*scale
POS3_EXP[8] = 90.0
GL_EXP[8] = (TR_EXP[8]/VEL_EXP[8])*scale



#print(f'{POS3_EXP[1]:g},{GL_EXP[1]:g}')


# Position vs time
# Lscale milimeter to meters scale factor
Lscale = 1000.0    


aP = 2.83 #0.92
bP = 0.64 #0.92
def Pft(t,aP,bP):
    Pt = aP*t**(bP)
    return Pt

# Liquidus Velocity vs position
aVL = 3.4 #2.40 #0.84
bVL = -0.56 #-0.305 #-0.086
def VLfp(P,aVL, bVL):
    VLp = aVL*exp(bVL*log(P))
    return VLp

# Liquidus cooling rate
aTR = 43.82 #263.0 #305.7
bTR = -0.78 #-1.44 #-1.46 

def TRfp(P,aTR, bTR):
    TRp = aTR*exp(bTR*log(P))
    return TRp


def grad_exp_fit(P):
    
    A = 10405.54
    B = -200.22
    C = 1.04692
    GLp = A + B*P + C*P**2

    return GLp


t_total = 230

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
    GL[i] = grad_exp_fit(pos)

    print(f'{t_:g}s  ----  {(posm):g}[m] --- {(VLm):g}[m/s] --- {TR_:g}[K/s] --- {GL_:g}[K/m]')

res = integrate.simpson(GL, x=POS, dx=1.0, axis=-1)
print(f'integral = {res:g}')
GL_mean = 1/(110-0) * res 
print(f'GL_mean = {GL_mean:g}[K.m^-1]')



# Calculation of mean mlCu

Tv[0] = 941.75
Cv[0] = 0
Tv[1] = 940
Cv[1] = 0.9313
Tv[2] = 938
Cv[2] = 2.173
Tv[3] = 936
Cv[3] = 3.425
Tv[4] = 934
Cv[4] = 5.029
Tv[5] = 933
Cv[5] = 6.022
Tv[6] = 932
Cv[6] = 7.326
Tv[7] = 931.6
Cv[7] = 7.928
Tv[8] = 931.3
Cv[8] = 9.888
Tv[9] = 932.35
Cv[9] = 11.933
Tv[10] = 920
Cv[10] = 17
Tv[11] = 915
Cv[11] = 18
Tv[12] = 910
Cv[12] = 19.82
Tv[13] = 905
Cv[13] = 21.2
Tv[14] = 900
Cv[14] = 22.5
Tv[15] = 880
Cv[15] = 26.82
Tv[16] = 860
Cv[16] = 30.3
Tv[17] = 840
Cv[17] = 33.22
Tv[18] = 831.84
Cv[18] = 34.347


dydx = diff(Tv) / diff(Cv)

for i in range(18):
    print(f'{i}, C = {Cv[i]} wt%Cu, mL = {dydx[i]:g}')


resmL = integrate.simpson(dydx[0:17], x=Cv[0:17], dx=1.0, axis=-1)
print(f'integral = {resmL:g}')
mL_Mean = 1/(Cv[17]-Cv[1]) * resmL
print(f'mLMean = {mL_Mean:g} [K.wt%^-1]')


kvCu[0] = 0.1236
CvCu[0] = 0.0
kvCu[1] = 0.1174
CvCu[1] = 5
kvCu[2] = 0.1078
CvCu[2] = 10
kvCu[3] = 0.09709
CvCu[3] = 15
kvCu[4] = 0.0936
CvCu[4] = 19.95
kvCu[5] = 0.0970
CvCu[5] = 25
kvCu[6] = 0.1075
CvCu[6] = 30
kvCu[7] = 0.12619
CvCu[7] = 34.35

print('\n Partition Coefficient of Cu')


for i in range(7):
    print(f'{i}, C = {CvCu[i]:g} wt%Cu, ko_Cu = {kvCu[i]:g}')


print('\n# Mean Partition Coefficient of Cu')
reskvCu= integrate.simpson(kvCu[0:6], x=CvCu[0:6], dx=1.0, axis=-1)
print(f'\nintegral = {reskvCu:g}')
kvCu_Mean = 1/(CvCu[6]-CvCu[1]) * reskvCu
print(f'kvCu_Mean = {kvCu_Mean:g} [-]')


# TL Nb

TvNb[0] = 925.82
CvNb[0] = 0
TvNb[1] = 951.1
CvNb[1] = 11.280
TvNb[2] = 950.1
CvNb[2] = 20
TvNb[3] = 948.5
CvNb[3] = 30
TvNb[4] = 946.16
CvNb[4] = 40
TvNb[5] = 942.5
CvNb[5] = 50
TvNb[6] = 939.41
CvNb[6] = 55
TvNb[7] = 934.95
CvNb[7] = 60
TvNb[8] = 927.4
CvNb[8] = 65
TvNb[9] = 911.9
CvNb[9] = 70
TvNb[10] = 874.3
CvNb[10] = 75
TvNb[11] = 848.0
CvNb[11] = 77
TvNb[12] = 831.85
CvNb[12] = 78.025


dydxNb = diff(TvNb) / diff(CvNb)

for i in range(12):
    print(f'{i}, C = {CvNb[i]} wt%Nb, mL = {dydxNb[i]:g}')


resmLNb = integrate.simpson(dydxNb[0:12], x=CvNb[0:12], dx=1.0, axis=-1)
print(f'integral = {resmLNb:g}')
mL_MeanNb = 1/(Cv[17]-Cv[1]) * resmL
print(f'mLMeanNb = {mL_MeanNb:g} [K.wt%^-1]')




# Nb Partition coefficient

kvNb[0] = 2.8338
CvNb[0] = 0
kvNb[1] = 2.188
CvNb[1] = 5
kvNb[2] = 2.049
CvNb[2] = 10
kvNb[3] = 2.059
CvNb[3] = 15
kvNb[4] = 2.106
CvNb[4] = 20
kvNb[5] = 2.162
CvNb[5] = 25
kvNb[6] = 2.238
CvNb[6] = 30
kvNb[7] = 2.322
CvNb[7] = 35
kvNb[8] = 2.458
CvNb[8] = 40
kvNb[9] = 2.620
CvNb[9] = 45
kvNb[10] = 2.878
CvNb[10] = 50
kvNb[11] = 3.264
CvNb[11] = 55
kvNb[12] = 3.969
CvNb[12] = 60
kvNb[13] = 5.565
CvNb[13] = 65
kvNb[14] = 10.00
CvNb[14] = 70
kvNb[15] = 39
CvNb[15] = 75
kvNb[16] = 86
CvNb[16] = 77
kvNb[17] = 132.08
CvNb[17] = 78.03

print('\n# Mean Partition Coefficient of Nb')
reskvNb= integrate.simpson(kvNb[0:17], x=CvNb[0:17], dx=1.0, axis=-1)
print(f'\nintegral = {reskvNb:g}')
kvNb_Mean = 1/(CvNb[6]-CvNb[1]) * reskvNb
print(f'kvNb_Mean = {kvNb_Mean:g} [-]')


"""
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

"""

fig2,ax2 = plt.subplots(1)
ax2.set_title(r'Liquidus Growth Rate $V_{L}$, alloy Al-3wt%Cu-0.5wt%Nb-0.1wt%Fe  ',{'color': 'black', 'fontsize': 14})
ax2.grid(color='k', which='both', linestyle='-', linewidth=0.5)
ax2.plot(POS[1:n],VL[1:n],'k-',linewidth = 1.5)
ax2.plot(POS1_EXP[1:9],VEL_EXP[1:9],'ks',linewidth = 1.5)
#ax2.plot(TIME_EXP[0:8],POS_EXP[0:8],'ks',linewidth = 1.5)
ax2.set_ylabel(r'Liquidus Growth Rate, $\mathit{V_{L}} \,\,\, \mathrm{ [mm.s^{-1}]}$',{'color': 'black', 'fontsize': 18})
ax2.set_xlabel(r'Position, $\,\mathit{P} \,\, \mathrm{ [mm]}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend(('Experimental Curve Fit','Experimental' ),frameon=True,edgecolor='k',shadow=True, loc= 'best', fontsize = 16 )
plt.axis([0,120,0,4])

fig3,ax3 = plt.subplots(1)
ax3.set_title(r'Liquidus Cooling Rate $T_{R}$, alloy Al-3wt%Cu-5wt%Nb-0.1wt%Fe  ',{'color': 'black', 'fontsize': 14})
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
ax4.set_title(r'Liquidus Thermal Gradient $G_{L}$, alloy Al-3wt%Cu-5wt%Nb-0.1wt%Fe  ',{'color': 'black', 'fontsize': 14})
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
plt.axis([0,120,0,20000])



plt.show()


