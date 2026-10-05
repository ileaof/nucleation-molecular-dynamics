#!/bin/bash
import matplotlib.pyplot as plt
from numpy import zeros, sqrt, exp, log, pi, cos, sin
from math import pi
import csv
from midpoint    import midpoint
# only to calculate enthalpy
from scipy import integrate
from scipy.optimize import fsolve, brentq
from derivative_shim import derivative
#

temp_int = zeros(100)
cp_int = zeros(100)


"""
I.L.Ferreira, A.L.S. Moreira. On the Continuous Mechanics First and Second‑Order Formulations for Nonequilibrium Nucleation: Derivation and Applications. Int. J. Thermophys. (2023)
DOI: 10.1007/s10765-023-03178-2

I.L. Ferreira. A Non-Equilibrium Nucleation Model to Calculate the Density of State and Its
   Application to the HeatCapacity of Stoichiometric UO2. Int. J. Thermophys. 42:148 (2021).
DOI: 10.1007/s10765-021-02903-z

I.L. Ferreira. Non-equilibrium Nucleation: Application to Solidification and Molar- Specific Heat Capacity of Pure Metals and Phases. Int. J. Thermophys. 43:33:1-25 (2022)
DOI: 10.1007/s10765-021-02956-0

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
rc_v = zeros(n)
dfgamdr_ana_v = zeros(n)
dfgamdr_num_v = zeros(n)
dTdr_v = zeros(n)
DT_v = zeros(n)
GT_hom_v = zeros(n)
GT_hom_v2 = zeros(n)
r_hom_v = zeros(n)
r_het_v = zeros(n)
theta_v = zeros(n)
termoA_v = zeros(n)
termoA_het_v = zeros(n)
GT_het_v = zeros(n)
V_hom_v = zeros(n)
V_het_v = zeros(n)
DoS_hom_v = zeros(n)
DoS_het_v = zeros(n)
gam_hom_v = zeros(n)
gam_het_v = zeros(n)

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



# Avogadro's number
Nav =  6.022E23
# Planck's Constant - [J.s]
h   =  6.626E-34
h_  =  h / (2.0*pi)
# Gas Universal Constant - [J/(mol.K)]
R   = 8.3144 
# Constant alfa  
alpha = 0.155
# Boltzmann's Constant - [J/K]
kB = 1.380658e-23 
# Avogadro Na - [atomos/mol]
Na = 6.022e+23
# speed of light - [m/s]
c = 299792458.
# Electron rest mass me
me = 9.1093837015e-31 # kg
# Wiedemann-Franz Law Constant
L = 2.44e-8 # W.Ohm/K^2
#small number
small = 1.e-31


# Debye Temperature - M.P. Mader, Condensed Matter Physics, 2000.
ThetaD_Al = 433.0
ThetaD_Cu = 347.0
ThetaD_Si = 645.0
ThetaD_Mg = 403.0
ThetaD_Au = 162.0
ThetaD_Ag = 227.0
ThetaD_C_gra = 1700.0
ThetaD_Fe = 477.0
ThetaD_Cr = 606.0
ThetaD_Ni = 477.0
# atomic mass
M_Al = 26.982e-3
M_Cu = 63.54e-3
M_Si = 28.086e-3
M_Mg = 24.312e-3
M_Au = 196.97e-3
M_Ag = 107.87e-3
M_C  = 12.011e-3
M_Fe = 55.847e-3
M_Cr = 51.906e-3
M_Ni = 58.71e-3
radius_Al = 143.1e-12 #118.0e-12
radius_Cu = 128.0e-12
radius_Si = 111.0e-12
radius_Mg = 160.0e-12
radius_Au = 144.0e-12
radius_Ag = 144.0e-12
radius_C  =  67.0e-12
radius_Fe = 156.0e-12
radius_Cr = 166.0e-12
radius_Ni = 149.0e-12
#mu0 = 1.25663706144e-6
Z_Al =+3.0 # real +3
Z_Cu =+2.0 # real +2 # can be 2
Z_Si =+4.0 # real +4
Z_Mg =+2.0 # real +2
Z_Au =+1.0 # real +3
Z_Ag =+1.0 # real +3
Z_C_gra = 0.0 #2.0, 4.0
Z_Fe = 2.0 #2, 3
Z_Cr = 6.0 #2, 3, 6
Z_Ni = 4.0 #+2 # # may be 0, +2 and +4

T_Al = 933.15
T_Si = 1687.15
T_Cu = 1357.77
T_Mg = 923.15
T_Au = 1337.33
T_Ag = 1234.93
T_C_gra = 3300.0
T_Fe = 1811.15
T_Cr = 2180.15
T_Ni = 1728.15

vDia = 12000.0 # [m/s]
vPhase  = 6400.0  #vAl = 5240
TL = 901.18 #
DeltaH = 260300.0
rhol = 2380.0 # a calcular # kg/m^3
rhos = 2549.14
Delta_Sv_hom = DeltaH * rhos / TL
#Equilibrium surface tension and surface stress
sigma_0 = 1.09 #0.914 #0.821 #1.1 #0.821
sigma = sigma_0
gamma_0 = 0.183
mu0 = 1.6
lambda0 = -5.0
rho0 = 2.549e-6
mu = 2.594e10
lambda_ = 5.034e10

rhet = 0.0


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




#    print(f'Pft({t_:.2g}) = {value:.2g}[mm]')

    


def molarfrac(T):

    if T < 300:
       xCu = 3.06615e-5
       xSi = 5.74931e-8
       xAl = 1 - xCu - xSi
    elif (T >= 300.0) and (T < 790.0):
       xCu = 1.256e-5   * exp(T/106.95095) - 5.22441e-4
       xSi = 1.79789e-6 * exp(T/91.018020) - 2.77861e-4
       xAl = 1 - xCu - xSi

    elif(T >= 790.0 ) and (T <= 890.0):
       xCu =    10458.04364 * 0.98348 ** T
       xSi =   -0.84061 + 0.00209 * T - 1.2869e-6 * T**2
       xAl = 1 - xCu - xSi
       
    elif T > 890.0:
       xCu = 0.00351
       xSi = 0.00431
       xAl = 1 - xCu - xSi
 
    return xAl, xCu, xSi



def func_r(r, req, lambda_, mu, sigma, lambda0, mu0, delta_H, rho, Tf, gamma_0, dTdx):
    delta_Sv = delta_H * rho/Tf
#    y  = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*r   + 2.0 * (lambda0 + mu0)
#    y0 = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*req + 2.0 * (lambda0 + mu0)
    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    surf_stress  =  -2.0 * sigma  * (3.0*lambda_ + 2.0*mu)/ ( (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)  ) * log ( y0 / y )
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress
    
    alpha = 2 * (sigma + 2*(lambda0 + mu0)) / ( r*(3*lambda_ + 2*mu) )
    dgammadr = -2*gamma_0 * req ** 2 / (r ** 3) - 2*sigma / (r*(1+alpha))

    DT = dgammadr / (delta_Sv - gamma_f / (4.0 * pi * r * r * dTdx))
 
#    func = gamma_f / (4.0 * pi * r * r * (delta_Sv - 1/DT * dgammadr) ) - dTdx


    func = gamma_f / (4.0 * pi * r * r * (delta_Sv ) ) - dTdx
    
    print(f'__________________________________')
    print(f'lambda_ = {lambda_} ')
    print(f'lambda0 = {lambda0} ')
    print(f'mu = {mu} ')
    print(f'mu0 = {mu0} ')

    print(f'surf_stress = {surf_stress} [N/m]')
    print(f'gamma_0 = {gamma_0} [J/m^2]')
    print(f'gamma_f = {gamma_f} [J/m^2]')
    print(f'dgammadr = {dgammadr} [J/m^3]')
    print(f'DT = {DT} [K]')
    print(f'r = {r} [m]')
    
    return func


def func_r_new(r, req, lambda_, mu, sigma, lambda0, mu0, delta_H, rho, Tf, gamma_0, dTdx):
    delta_Sv = delta_H * rho/Tf
#    y  = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*r   + 2.0 * (lambda0 + mu0)
#    y0 = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*req + 2.0 * (lambda0 + mu0)
    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    surf_stress  =  -2.0 * sigma  * (3.0*lambda_ + 2.0*mu)/ ( (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)  ) * log ( y0 / y )
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress
    
    alpha = 2 * (sigma + 2*(lambda0 + mu0)) / ( r*(3*lambda_ + 2*mu) )
    dgammadr = -2*gamma_0 * req ** 2 / (r ** 3) - 2*sigma / (r*(1+alpha))

    DT = dgammadr / (delta_Sv - gamma_f / (4.0 * pi * r * r * dTdx))
 
    func = gamma_f / (delta_Sv - 1/DT * dgammadr) - dTdx * 4.0 * pi * r * r 
    
    return func




def func_r_theta(r, req, lambda_, mu, sigma, lambda0, mu0, delta_H, rho, Tf, gamma_0, theta, dTdx):
    delta_Sv = delta_H * rho/Tf
#    y  = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*r   + 2.0 * (lambda0 + mu0)
#    y0 = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*req + 2.0 * (lambda0 + mu0)
    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    surf_stress  =  -2.0 * sigma  * (3.0*lambda_ + 2.0*mu)/ ( (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)  ) * log ( y0 / y )
    surf_stress_theta = surf_stress * ( 2 ) / (1 - cos(theta))
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress_theta
    
    #func_theta = gamma_f / (4.0 * pi * r * r * delta_Sv) - dTdx
    
    return func_theta



def gibbs_thomson(gamma, delta, req, delta_H, rho, Tf, theta):
    """Returns de Gibbs-Thomson coefficient for Equilibrium and non-Equilibrium nucleation. For the equilibrium Gibbs-Thomson, simply make delta = 0, req = 1e6.
    Input data:
    gamma  -  equilibrium surface energy;
    req     - Equilibrium radius;
    delta_H - Latent heat;
    rho     - Density of the nuclating phase in the nucleation temperature;
    Tf      - Transformation temperature, normally temperature of fusion;
    """
    delta_Sv = delta_H * rho/Tf * (2-3*cos(theta)+cos(theta)**3)/4
    delta_Sv_f_delta = (((1.0 - delta/req) ** 2.0) * delta_Sv)
    GT = (gamma / delta_Sv_f_delta)
    return GT, delta_Sv
    

def sigma_func(req, r, lambda_, mu, sigma, lambda0, mu0):

    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    tau =  sigma  +  2 * sigma  * log ( y0 / y )
    surf_stress  =  -2 * sigma  * log ( y0 / y )

    return tau, surf_stress


def gamma_func(gamma_0, delta, req, surf_stress):
    """Return the surface energy as a function of surface energy in equilibrium, delta, equilibrium radius and the surface stress."""
    gamma_f  = gamma_0   / ((1.0 - delta/req) ** 2.0) - surf_stress
    return gamma_f


def gamma_func_r(gamma_0, r, req, surf_stress):
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress
    return gamma_f

# function test for derivative
# means better order the function variables
#
def fgam(r, req, lambda_, mu, sigma, lambda0, mu0):
    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    surf_stress  =  -2 * sigma  * log ( y0 / y )
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress
    return gamma_f
#



def dgammadr_func_r(req, r, gamma_0, lambda_, mu, sigma, lambda0, mu0):
    alpha = 2 * (sigma + 2*(lambda0 + mu0)) / ( r*(3*lambda_ + 2*mu) )
    dgammadr = -2*gamma_0 * req ** 2 / (r ** 3) - 2*sigma / (r*(1+alpha))
    return dgammadr
 

def ftheta(theta):
    f_theta = (2 - 3 * cos(theta) + cos(theta) ** 3)
    return f_theta

def dfthetadr(theta):
    dfdr = (-3*(2 - 3*cos(theta) + cos(theta)**3) - (1-cos(theta))*(2-cos(theta)-cos(theta)**2) )
    return dfdr


def rc_het(rc_hom, theta):
    from scipy.integrate import quad
    def f(theta):
        f_ = ( sin(theta) * (1+cos(theta)) ) / ( 2 - cos(theta) - cos(theta)**2 )
        return f_
    sol, err = quad(f,pi,theta)
#    rc_star = rc_hom / exp(-sol)
    rc_het = rc_hom / exp(-sol)
    if rc_het > 1e-15:
       rhet = rc_het
    print(f'integral = {sol}, error = {err}, rc_het = {rc_het}m')
    
    
    return rhet
    

def gibbs_thomson_r(gammasl, delta_Sv_hom, DT_, dgammadr, ftheta, dfthetadr, DTdx):
    """
    Returns de Gibbs-Thomson coefficient for Equilibrium and non-Equilibrium nucleation. For the equilibrium Gibbs-Thomson, simply make delta = 0, req = 1e6.
    Input data:
    gammasl        - Surface energy [J/m^2];
    delta_Sv_hom   - Volume entropy [J/(kg.K)];
    DT_        - Supercooling [K];
    dgammadr   - derivative of gamma in relation to r [J/m^3];
    ftheta     - Function f(theta) = 2 - 3*cos(theta) + cos(theta)**3 ;
    dfthetadr  - Derivative of f(theta) in respect to r
                 dfthetadr = (-3*(2 - 3*cos(theta) + cos(theta)**3) - (1-cos(theta))*(2-cos(theta)-cos(theta)**2) ) / 4;
    """
    delta_Ss_hom = 1/DT_ * dgammadr
    GT = gammasl*ftheta / (( delta_Sv_hom - delta_Ss_hom  )*ftheta - gammasl/DT_ * dfthetadr) - DTdx
    delta_Ss_het = 1/DT_ * dgammadr * ftheta
    delta_Sv_het =  delta_Sv_hom * ftheta
      
    return GT, delta_Sv_hom, delta_Ss_hom, delta_Sv_het, delta_Ss_het


def gibbs_thomson_hom_r(gammasl, delta_Sv_hom, DT, dgammadr):
    """
    Returns de Gibbs-Thomson coefficient for Equilibrium and non-Equilibrium nucleation. For the equilibrium Gibbs-Thomson, simply make delta = 0, req = 1e6.
    Input data:
    gammasl        - Surface energy [J/m^2];
    delta_Sv_hom   - Volume entropy [J/(kg.K)];
    DT        - Supercooling [K];
    dgammadr   - derivative of gamma in relation to r [J/m^3];
    ftheta     - Function f(theta) = 2 - 3*cos(theta) + cos(theta)**3 ;
    dfthetadr  - Derivative of f(theta) in respect to r
                 dfthetadr = (-3*(2 - 3*cos(theta) + cos(theta)**3) - (1-cos(theta))*(2-cos(theta)-cos(theta)**2) ) / 4;
    """
    
    delta_Ss_hom = 1/DT * dgammadr
    GT  = gammasl / ( delta_Sv_hom - delta_Ss_hom )
    return GT

def gibbs_thomson_het_r(gammasl, delta_Sv_hom, DT, dgammadr, theta):
    """
    Returns de Gibbs-Thomson coefficient for Equilibrium and non-Equilibrium nucleation. For the equilibrium Gibbs-Thomson, simply make delta = 0, req = 1e6.
    Input data:
    gammasl        - Surface energy [J/m^2];
    delta_Sv_hom   - Volume entropy [J/(kg.K)];
    DT        - Supercooling [K];
    dgammadr   - derivative of gamma in relation to r [J/m^3];
    ftheta     - Function f(theta) = 2 - 3*cos(theta) + cos(theta)**3 ;
    dfthetadr  - Derivative of f(theta) in respect to r
                 dfthetadr = (-3*(2 - 3*cos(theta) + cos(theta)**3) - (1-cos(theta))*(2-cos(theta)-cos(theta)**2) ) / 4;
    """
    
    delta_Ss_hom = 1/DT * dgammadr
    GT  = gammasl*ftheta(theta) / (( delta_Sv_hom - delta_Ss_hom )*ftheta(theta)/4 - gammasl/DT * dfthetadr(theta))
    return GT



def DA(T):
    EA =  25800    # [J/mol]
    DoA = 1.78e-7  # [m^2/s]
    D = DoA * exp(-EA/(R * T))
    return D
    
def CL(rhol, M):
    Cl = rhol * Na / M
    return Cl
    
def a(Cl):
    aL = (4/3 * pi / Cl)**(1/3)
    return aL

def IHom(r_, gamsl, DSv, dgamdr, M, rhol , dT_, Temp_):

    DA_ = DA(Temp_)
    CL_ = CL(rhol, M)
    a_  = a(CL_)
    termN1 = -16/3 * pi
    termN2 = gamsl**3 * (DSv - 3/dT_ * dgamdr)
    termD =  dT_**2 * kB * Temp_ * (DSv - 3/dT_ * dgamdr)**3
    #Ihom = (DA_/(a_**2))*(4 * pi * (r_**2)/(a_**2)) * CL_ * exp(termN1 * termN2 / termD)
    Ihom = (DA_/(a_**2))*(4 * pi * (r_**2)/(a_**2)) * CL_ * exp(-16/3*pi * gamsl**3 * (DSv - 3/dT_ * dgamdr)/(dT_**2*((DSv - 3/dT_ * dgamdr)**3 * kB * Temp_)))
    multiplier = (DA_/(a_**2))*(4 * pi * (r_**2)/(a_**2)) * CL_
    Ihom = 10e40 * exp(-16/3*pi * gamsl**3 * (DSv - 3/dT_ * dgamdr)/(dT_**2*((DSv - 3/dT_ * dgamdr)**3 * kB * Temp_)))
    print('_______________________')
    print(f'pi = {pi}')
    print(f'kB = {kB}')
    print(f'multiplier = {multiplier}')
    print(f'r_ = {r}')
    print(f'gamsl = {gamsl}')
    print(f'DSv = {DSv}')
    print(f'dgamdr = {dgamdr}')
    print(f'M = {M}')
    print(f'dT_ = {dT_}')
    print(f'Temp_ = {Temp_}')
    print(f'DA_ = {DA_}')
    print(f'CL_ = {CL_}')
    print(f'a_ = {a_}')
    
    return Ihom

def f(G, gam, dgamdr, DSv, r ):
    fu = dgamdr / (DSv - gam/(4 * pi * r * r * G)) - 8 * pi * r * G
    return fu


def f_DT(DT, gam, dgamdr, DSv, r ):
    f_ = dgamdr / (DSv - 2*gam/(r * DT)) - DT
    return f_

def f_het(theta, rc_hom, termoa, termob, termoc ):
    fA= termoc / (termoa * ftheta(theta)/4 - termob*dfthetadr(theta)) - rc_het(rc_hom, theta)

    return fA


#Phase mean molar composition or use a curve as a function of temperature provided
x_Cu = 6.33256e-3
x_Si = 3.93038e-3
x_Mg = 0.0
x_C  = 0.0
x_Fe = 0.0
x_Cr = 0.0
x_Ni = 0.0

x_Al    = 1.0 - x_Cu - x_Si - x_Mg - x_C - x_Fe - x_Cr - x_Ni
Z_Phase = x_Al * Z_Al + x_Cu * Z_Cu + x_Si * Z_Si + x_Mg * Z_Mg + x_C * Z_C_gra +  x_Fe * Z_Fe + x_Cr * Z_Cr + x_Ni * Z_Ni
M       = x_Al * M_Al + x_Cu * M_Cu + x_Si * M_Si + x_Mg * M_Mg + x_C * M_C + x_Fe * M_Fe + x_Cr * M_Cr + x_Ni * M_Ni
# Warning! Debye temperature of the phase! must be calculated at low temperatures, as it is so simple!
# Some situation could be approximated!
ThetaD_Phase = x_Al * ThetaD_Al + x_Cu * ThetaD_Cu + x_Si * ThetaD_Si + x_Mg * ThetaD_Mg + x_C * ThetaD_C_gra + x_Cr * ThetaD_Cr + x_Ni * ThetaD_Ni

wD_Phase =  kB * ThetaD_Phase / h_
wD_Al    =  kB * ThetaD_Al / h_
wD_Cu    =  kB * ThetaD_Cu / h_
wD_Si    =  kB * ThetaD_Si / h_
wD_Mg    =  kB * ThetaD_Mg / h_
wD_Au    =  kB * ThetaD_Au / h_
wD_Ag    =  kB * ThetaD_Ag / h_
wD_C_gra =  kB * ThetaD_C_gra / h_
wD_Fe    =  kB * ThetaD_Fe / h_
wD_Cr    =  kB * ThetaD_Cr / h_
wD_Ni    =  kB * ThetaD_Ni / h_
kD = wD_Phase/vPhase

GT_ref = 2.08e-7 # Phase
req = 6.3e-6 #5.96e-6 #3.79687e-06 #5.56e-6
dreq = 0.01 * req
dfdx = 0.0
#DT = 0.09
### Solid state calculation
N_V  = (1.0/(6.0*pi*pi))*(ThetaD_Phase*kB/(h_*vPhase))**3.0
N_Vi = Z_Phase * N_V
N_V_Al = x_Al * N_V
N_V_Cu = x_Cu * N_V
N_V_Si = x_Si * N_V
P = 0.5
kF_Phase = (3 * pi ** 2 * Z_Phase * N_V)**(1/3)
kF_Al = (3 * pi ** 2 * Z_Al * N_V_Al)**(1/3)
kF_Cu = (3 * pi ** 2 * Z_Cu * N_V_Cu)**(1/3)
kF_Si = (3 * pi ** 2 * Z_Si * N_V_Si)**(1/3)
vF_Phase = h_ * kF_Phase / me
vF_Al = h_ * kF_Al / me
vF_Cu = h_ * kF_Cu / me
vF_Si = h_ * kF_Si / me
EF_Phase = 1/P * h_ ** 2 * kF_Phase**2 / (2*me)
EF_Al = 1/P * h_ ** 2 * kF_Al**2 / (2*me)
EF_Cu = 1/P * h_ ** 2 * kF_Cu**2 / (2*me)
Eg_Si = 1.17 /  6.241509e18
EF_Si = (Eg_Si / 2) * (N_V_Si)**(2/3) / (N_V)**(2/3)
T_F_Phase = EF_Phase / kB
T_F_Al = EF_Al / kB
T_F_Cu = EF_Cu / kB
T_F_Si = EF_Si / kB
term_cve_Al = Z_Al*(ThetaD_Al**3.0)/T_F_Al
term_cve_Cu = Z_Cu*(ThetaD_Cu**3.0)/T_F_Cu
term_cve_Si = Z_Si*(ThetaD_Si**3.0)/T_F_Si




#N_Vi = N_V * Z_Al
#wD = kB * ThetaD_Al / h_
#kD = wD/vAl
### Fermi calculations
#kF = (3 * pi ** 2 * Z_Al * N_V)**(1/3)
#vF = h_ * kF / me
#P = 0.5
#EF = 1/P * h_ ** 2 * kF**2 / (2*me)
#T_F = EF / kB
#term_cve = Z_Al*(ThetaD_Al**3.0)/T_F

### Solid state

#args = (req, lambda_, mu, sigma_0, lambda0, mu0)
#args =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
#args2 =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
j = 6
dTdx = GL[j]
r    = req - dreq
args = (req, lambda_, mu, sigma_0, lambda0, mu0, DeltaH, rhos, TL, gamma_0, dTdx)
print(f'req = {req:g}m')
print(f'lambda_ = {lambda_:g}, mu = {mu:g}, simga_0 = {sigma_0:g}, lambda0 = {lambda0:g}')
print(f'mu0 = {mu0:g},  DeltaH = {DeltaH:g}, rhos = {rhos:g}, TL = {TL:g}, gamma_0 = {gamma_0:g},dTdx = {dTdx:g}')
print(args)
#func_r_new(r, req, lambda_, mu, sigma, lambda0, mu0, delta_H, rho, Tf, gamma_0, dTdx)
r_hom = brentq(func_r_new,r,req, args=args, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)
sigm, surfstress = sigma_func(req, r_hom, lambda_, mu, sigma_0, lambda0, mu0)
gam = gamma_func_r(gamma_0, r_hom, req, surfstress)
args = (req, lambda_, mu, sigma_0, lambda0, mu0)
dfgamdr_num = derivative(fgam, r_hom, dx=1e-10,args=args)
dfgamdr_ana = dgammadr_func_r(req, r_hom, gamma_0, lambda_, mu, sigma_0, lambda0, mu0)
args2 =   (gam, dfgamdr_ana, Delta_Sv_hom, r_hom)
DT = brentq(f_DT,0.0001,1, args=args2, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)
GT  = gibbs_thomson_hom_r(gam, Delta_Sv_hom, DT, dfgamdr_ana)
print(f'GL[{j}] = {GL[j]:g}[K/m], r_hom = {r_hom:g}[m],sigma_hom = {sigm:g}[N.m^-1], gamma = {gam:g}[J.m^-2]')
print(f'dfgamdr_num = {dfgamdr_num:g}[J.m^-3], dfgamdr_ana = {dfgamdr_ana:g}[J.m^-3]')
print(f'DT = {DT:g}[K], GT_hom = {GT:g}[m.K]')
#Heterogeneous nucleation
TERMOA = Delta_Sv_hom - 1/DT * dfgamdr_ana
TERMOB = gam/DT
TERMOC = 2 * gam / DT
args3 = (r, TERMOA, TERMOB, TERMOC)
thet = brentq(f_het,0.01,pi, args=args3, xtol=2e-14, rtol=8.881784197001252e-14, maxiter=16, full_output=False, disp=True)
rhet = rc_het(r_hom, thet)
gam_het = gam * ftheta(thet)
print(f'r = {r_hom:.4e} [m], σ = {sigm:.5f} [N/m], s = {surfstress:.5f} [N/m], γ = {gam:.5f} [J/m^2], γ_het = {gam_het:.5f} [J/m^2], dγdr(num) = {dfgamdr_num:.5f} [J/m^3], dγdr(ana) = {dfgamdr_ana:.5f} [J/m^3], dTdx = {dTdx:.2f} [K/m], DT = {DT:.5f} [K], GT_Hom = {GT:.5e} [K.m], TERMOA = {TERMOA}, TERMOB = {TERMOB},thet = {thet:.4f} [radian], {thet*180/pi:.4f} [degrees], rhet = {rhet:.4e} [m]')
GT_het = gibbs_thomson_het_r(gam, Delta_Sv_hom, DT,dfgamdr_ana,thet)
print(f'GT_hom = {GT:g}[m.K],GT_het = {GT_het:g}[m.K]')


fig1,ax1 = plt.subplots(1)
ax1.set_title(r'Liquidus Position $P$ function of $t$  $\mathrm{S_{L} = f(t)}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.65wt%Fe   ',{'color': 'black', 'fontsize': 14})
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
ax2.set_title(r'Liquidus Growth Rate $V_{L}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.65wt%Fe  ',{'color': 'black', 'fontsize': 14})
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
ax3.set_title(r'Liquidus Cooling Rate $T_{R}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.65wt%Fe  ',{'color': 'black', 'fontsize': 14})
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
ax4.set_title(r'Liquidus Thermal Gradient $G_{L}$, alloy Al-0.8wt%Si-0.6wt%Mg-0.65wt%Fe  ',{'color': 'black', 'fontsize': 14})
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

print("### Script has been interrupted here ###")
quit()

for i in range(0,36):
    r = req - i * dreq
    rc_v[i] = r
    sigm, surfstress = sigma_func(req, r, lambda_, mu, sigma_0, lambda0, mu0)
    gam = gamma_func_r(gamma_0, r, req, surfstress)
    gam_hom_v[i] = gam
    args = (req, lambda_, mu, sigma_0, lambda0, mu0)
    dfgamdr_num = derivative(fgam, r, dx=1e-10,args=args)
    dfgamdr_ana = dgammadr_func_r(req, r, gamma_0, lambda_, mu, sigma_0, lambda0, mu0)
    dfgamdr_num_v[i] = dfgamdr_num
    dfgamdr_ana_v[i] = dfgamdr_ana
    args =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
#    print(f'gam = {gam}, dfgamdr_ana = {dfgamdr_ana}, Delta_Sv_hom = {Delta_Sv_hom}, r = {r}')
#    G = f(1000, gam, dfgamdr_ana, Delta_Sv_hom, r )
    dTdr = brentq(f,1,30000, args=args, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)
    
    args2 =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
    DT = brentq(f_DT,0.0001,1, args=args2, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)
#    DT = dfgamdr_ana / (Delta_Sv_hom - gam/(4 * pi * r * r * dTdr)) # DT = 8 * pi * r  * dTdr
    DT_v[i] = DT
    GT  = gibbs_thomson_hom_r(gam, Delta_Sv_hom, DT, dfgamdr_ana)
    GT2 = r * DT / 2 #dTdr * 4 * pi * r ** 2
    GT_hom_v[i]  = GT
    GT_hom_v2[i] = GT2

##    dTdr = 3/(8 * pi * r) * DT
    dTdr_v[i] = dTdr
#    DT_  = dfgamdr_ana / (Delta_Sv_hom - gam / (4 * pi * r * r * dTdr))
    
#    print(f'counter = {i}, r = {r:.4e} [m], σ = {sigm:.5f} [N/m], s = {surfstress:.5f} [N/m], γ = {gam:.5f} [J/m^2], dγdr(num) = {dfgamdr_num:.5f} [J/m^3], dγdr(ana) = {dfgamdr_ana:.5f} [J/m^3], dTdr = {dTdr:.2f} [K/m], DT_ = {DT_:.2f} [K]')
    TERMOA = Delta_Sv_hom - 1/DT * dfgamdr_ana
    TERMOB = gam/DT
    TERMOC = 2 * gam / DT
    args3 = (r, TERMOA, TERMOB, TERMOC)
    thet = brentq(f_het,0.01,pi, args=args3, xtol=2e-14, rtol=8.881784197001252e-14, maxiter=16, full_output=False, disp=True)
    rhet = rc_het(r, thet)
    gam_het_v[i] = gam * ftheta(thet)
    print(f'counter = {i}, r = {r:.4e} [m], σ = {sigm:.5f} [N/m], s = {surfstress:.5f} [N/m], γ = {gam:.5f} [J/m^2], dγdr(num) = {dfgamdr_num:.5f} [J/m^3], dγdr(ana) = {dfgamdr_ana:.5f} [J/m^3], dTdr = {dTdr:.2f} [K/m], DT = {DT:.5f} [K], GT_Hom = {GT:.5e} [K.m],GT_Hom2 = {GT2:.5e} [K.m], TERMOA = {TERMOA}, TERMOB = {TERMOB},thet = {thet:.4f} [radian], {thet*180/pi:.4f} [degrees], rhet = {rhet:.4e} [m]')
    r_hom_v[i] = r
    r_het_v[i] = rhet
    theta_v[i] = thet
    termoA_v[i] = TERMOA
    termoA_het_v[i] = ( TERMOA * ftheta(thet) - gam/DT * dfgamdr_ana )
    GT_het_v[i] = gibbs_thomson_het_r(gam, Delta_Sv_hom, DT,dfgamdr_ana,thet)
    V_hom_v[i] = (4.0 / 3.0 * pi * r ** 3.0)
    V_het_v[i] = (1.0 / 3.0 * pi * rhet ** 3.0) * (2 - 3 * cos(thet) + cos(thet) ** 3 )
    DoS_hom_v[i] = V_hom_v[i] * wD_Phase ** 2 / (2 * pi * pi * vPhase ** 3)
    DoS_het_v[i] = V_het_v[i] * wD_Phase ** 2 / (2 * pi * pi * vPhase ** 3)

print("\nN_V = %g [modes/m^3] \n" % (N_V))
print("\n(N_V)i = N_V*Z_Phase = %g [conducting electrons/m^3] \n" % (N_V*Z_Phase))
print("\nwD_Phase = %g [Hz] \n" % (wD_Phase))
print("\nkD = %g [cycles/m] \n" % (kD))
print("\nkF_Al = %g \n" % (kF_Al))
print("\nkF_Cu = %g \n" % (kF_Cu))
print("\nkF_Si = %g \n" % (kF_Si))
print("\nvF_Al = h_ * kF_Al / me == > vF_Al = %g m/s \n" % (vF_Al))
print("\nvF_Cu = h_ * kF_Cu / me == > vF_Cu = %g m/s \n" % (vF_Cu))
print("\nvF_Si = h_ * kF_Si / me == > vF_Si = %g m/s \n" % (vF_Si))
print("\nEF_Al = %g J = %g eV \n" % (EF_Al,EF_Al*6.241509e18))
print("\nEF_Cu = %g J = %g eV \n" % (EF_Cu,EF_Cu*6.241509e18))
print("\nEF_Si = %g J = %g eV \n" % (EF_Si,EF_Si*6.241509e18))
print("\nT_F_Al = %g K\n" % (T_F_Al))
print("\nT_F_Cu = %g K\n" % (T_F_Cu))
print("\nT_F_Si = %g K\n" % (T_F_Si))

print("\nterm_cve_Al = %g " % (term_cve_Al))
print("\nterm_cve_Cu = %g " % (term_cve_Cu))
print("\nterm_cve_Si = %g " % (term_cve_Si))

fig1,ax1 = plt.subplots(1)
#ax1.set_yscale('log')
ax1.set_title(r'Surface energy derivative, phase FCC_A1 $\frac{d\gamma}{dr}$ ',{'color': 'black', 'fontsize': 14})
ax1.plot(1-rc_v[1:91]/req,dfgamdr_ana_v[1:91],'k-',linewidth = 1.5)
ax1.plot(1-rc_v[1:91]/req,dfgamdr_num_v[1:91],'b.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume21[1:94]/1e-6,'g-.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume22[1:94]/1e-6,'b-.',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Al[1:89]/EF_Al,'k-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Cu[1:89]/EF_Cu,'g-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Si[1:89]/EF_Si,'r-',linewidth = 1.5)
ax1.set_ylabel(r'$\frac {d\gamma} {dr}\,\, \mathrm{ [J.m^{-3}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax1.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
plt.axis([0,0.35,-0.2e7,0])
#plt.axis([0,1,-7e7,0])
#plt.axis([0,900,8.5,13.0])



fig2,ax2 = plt.subplots(1)
#ax1.set_yscale('log')
ax2.set_title(r'Thermal gradient $\mathrm{\nabla T}$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax2.plot(1-rc_v[1:91]/req,dTdr_v[1:91],'k-',linewidth = 1.5)
#ax2.plot(1-rc_v[1:91]/req,dfgamdr_num_v[1:91],'b.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume21[1:94]/1e-6,'g-.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume22[1:94]/1e-6,'b-.',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Al[1:89]/EF_Al,'k-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Cu[1:89]/EF_Cu,'g-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Si[1:89]/EF_Si,'r-',linewidth = 1.5)
ax2.set_ylabel(r'$\nabla T \,\, \mathrm{ [K.m^{-1}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax2.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
#plt.axis([0,1,0,12000])
plt.axis([0,0.35,0,10000])
#plt.axis([0,900,8.5,13.0])


fig3,ax3 = plt.subplots(1)
#ax1.set_yscale('log')
ax3.set_title(r'Undercooling $\mathrm{\Delta T}$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax3.plot(1-rc_v[1:91]/req,DT_v[1:91],'b-',linewidth = 1.5)
#ax2.plot(1-rc_v[1:91]/req,dfgamdr_num_v[1:91],'b.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume21[1:94]/1e-6,'g-.',linewidth = 1.5)
#ax1.plot(T_v[1:94],molarVolume22[1:94]/1e-6,'b-.',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Al[1:89]/EF_Al,'k-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Cu[1:89]/EF_Cu,'g-',linewidth = 1.5)
#ax1.plot(T_v[1:89],1-chem_pot_Si[1:89]/EF_Si,'r-',linewidth = 1.5)
ax3.set_ylabel(r'$\mathrm{\Delta T}\,\, \mathrm{ [K]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax3.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax3.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax3.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
#plt.axis([0,1,0,12000])
plt.axis([0,0.35,0,0.9])
#plt.axis([0,900,8.5,13.0])


fig4,ax4 = plt.subplots(1)
#ax4.set_yscale('log')
ax4.set_title(r'Homogeneous Gibbs-Thomson, phase FCC_A1 $\mathrm{\Gamma ^{Hom}}$ ',{'color': 'black', 'fontsize': 14})
ax4.plot(1-rc_v[1:35]/req,GT_hom_v[1:35],'k-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax4.set_ylabel(r'$\mathrm{\Gamma ^{Hom}}\,\, \mathrm{ [m.K]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax4.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax4.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
#plt.axis([0,1,0,12000])
plt.axis([0,0.35,1e-8,1e-6])
#plt.axis([0,900,8.5,13.0])


fig5,ax5 = plt.subplots(1)
#ax4.set_yscale('log')
ax5.set_title(r'Nucleation radius r, phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax5.plot(dTdr_v[1:35],r_hom_v[1:35]*1e6,'b-',linewidth = 1.5)
ax5.plot(dTdr_v[1:35],r_het_v[1:35]*1e6,'r-',linewidth = 1.5)

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax5.set_ylabel(r'Nucleation radius, r $\mathrm{[\mu m]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax5.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax5.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax5.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.legend((r'Homogeneous',r'Heterogeneous' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 16 )
#plt.axis([0,1,0,12000])
plt.axis([0,7000,3,8])
#plt.axis([0,900,8.5,13.0])


fig6,ax6 = plt.subplots(1)
#ax4.set_yscale('log')
ax6.set_title(r'Nucleation angle $\mathrm{\theta}$, phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax6.plot(dTdr_v[1:35],theta_v[1:35],'b-',linewidth = 1.5)
##ax5.plot(dTdr_v[1:36],r_het_v[1:36]*1e6,'r-',linewidth = 1.5)

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax6.set_ylabel(r'Nucleation angle, $\mathrm{\theta \,\, [Radian]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax6.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax6.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax6.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
#plt.legend((r'Homogeneous',r'Heterogeneous' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 16 )
#plt.axis([0,1,0,12000])
plt.axis([0,7000,1.1,1.7])
#plt.axis([0,900,8.5,13.0])


fig7,ax7 = plt.subplots(1)
#ax4.set_yscale('log')
ax7.set_title(r'Nucleation Volume and Surface Entropies, $\mathrm{\Delta S_{V} \, , \Delta S_{S} }$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax7.plot(dTdr_v[1:35],termoA_v[1:35],'b-',linewidth = 1.5)
ax7.plot(dTdr_v[1:35],termoA_het_v[1:35],'r-',linewidth = 1.5)

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax7.set_ylabel(r' $\mathrm{\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right ) \,\, and \,\, \left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )\,f \left(\theta \right )  - \frac {\gamma _{SL} } {\Delta T} \frac {\partial f\left ( \theta \right ) } {\partial r} \,\, [J.m^{3}.K^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax7.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax7.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax7.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.legend((r'Homogeneous',r'Heterogeneous' ),frameon=True,edgecolor='k',shadow=True, loc= 'center right', fontsize = 16 )
#plt.axis([0,1,0,12000])
#plt.axis([0,7000,0.9e6,6e6])
#plt.axis([0,900,8.5,13.0])


fig8,ax8 = plt.subplots(1)
#ax4.set_yscale('log')
ax8.set_title(r'Homogeneous and Heterogeneous Gibbs-Thomson, $\mathrm{\Gamma ^{Hom} \, , \Gamma ^{Het} }$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax8.plot(dTdr_v[1:35],GT_hom_v[1:35],'b-',linewidth = 1.5)
ax8.plot(dTdr_v[1:35],GT_het_v[1:35],'r-',linewidth = 1.5)
ax8.axhline(GT_ref,0,6475,color='r', linestyle='--',lw=1.5)

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax8.set_ylabel(r'$\mathrm{Gibbs-Thomson, \,\, \Gamma \,\, [m.K]}$',{'color': 'black', 'fontsize': 16})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax8.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax8.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax8.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.legend((r'Homogeneous, $\mathrm{\Gamma ^{Hom} = \frac {\gamma _{SL}} {\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )}}$',r'Heterogeneous, $\mathrm{\Gamma ^{Het} = \frac {\gamma _{SL} \, f \left(\theta \right ) } {\left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )\,f \left(\theta \right ) - \frac {\gamma _{SL} } {\Delta T} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 15 )
#plt.axis([0,1,0,12000])
plt.axis([0,7000,1e-9,3.0e-6])
#plt.axis([0,900,8.5,13.0])


fig9,ax9 = plt.subplots(1)
#ax4.set_yscale('log')
ax9.set_title(r'Homogeneous and Heterogeneous Nuclei Volume, $\mathrm{V ^{Hom} \, , V ^{Het} }$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax9.plot(dTdr_v[1:35],V_hom_v[1:35],'b-',linewidth = 1.5)
ax9.plot(dTdr_v[1:35],V_het_v[1:35],'r-',linewidth = 1.5)
#ax8.axhline(0.9e-7,0,6475,color='r', linestyle='--',lw=1.5)

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax9.set_ylabel(r'$\mathrm{Nuclei\,\,Volume, \,\, V \,\, [m^{3}]}$',{'color': 'black', 'fontsize': 16})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax9.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax9.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax9.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.legend((r'Homogeneous',r'Heterogeneous'),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 16 )
ax9.yaxis.get_offset_text().set_fontsize(14)
plt.ticklabel_format(axis='y', style='sci')
#plt.axis([0,1,0,12000])
plt.axis([0,7000,1e-17,12e-16])
#plt.axis([0,900,8.5,13.0])


fig10,ax10 = plt.subplots(1)
#ax4.set_yscale('log')
ax10.set_title(r'Density of State, $D \left ( \omega \right ) $, phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax10.plot(dTdr_v[1:35],DoS_hom_v[1:35],'b-',linewidth = 1.5)
ax10.plot(dTdr_v[1:35],DoS_het_v[1:35],'r-',linewidth = 1.5)

ax10.set_ylabel(r'$\mathrm{Density\,of\,State,\, D \left ( \omega \right )}$',{'color': 'black', 'fontsize': 16})

ax10.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax10.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax10.get_yticklabels(), fontsize=16, fontweight="normal")
plt.text(4060,0.7, r' $\mathrm{\frac{N}{V}\,\, = %g \,[modes.m^{-3}]} $' % (N_V), {'color': 'black', 'fontsize': 14})
plt.text(4000,0.65, r' $\mathrm{\left ( \frac{N}{V} \right )_{i} = %g \,[conducting\, e^{-}.m^{-3}]} $' % (N_Vi), {'color': 'black', 'fontsize': 14})
plt.text(4060,0.6, r' $\mathrm{E^{Fermi}_{Phase}}$ = %g [J] = %g [eV]' % (EF_Phase,EF_Phase*6.241509e18), {'color': 'black', 'fontsize': 14})
plt.text(4060,0.55, r' $\mathrm{T^{Fermi}_{Phase}}$ = %g [K] ' % (T_F_Phase), {'color': 'black', 'fontsize': 14})
plt.text(4060,0.5, r' $\mathrm{V^{Fermi}_{Phase}}$ = %g [$\mathrm{m.s^{-1}}$] ' % (vF_Phase), {'color': 'black', 'fontsize': 14})

plt.text(4060,0.45, r' $\mathrm{\omega _{D}}$ = %g [Hz]' % (wD_Phase), {'color': 'black', 'fontsize': 14})

plt.legend((r'Homogeneous, $D \left ( \omega \right ) = \frac {V^{Hom} \omega ^{2}} {2 \pi ^{2} \nu ^{3} }$',r'Heterogeneous, $D \left ( \omega \right ) = \frac {V^{Het} \omega ^{2}} {2 \pi ^{2} \nu ^{3} }$'),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 16 )

plt.axis([0,7000,0,1])


fig11,ax11 = plt.subplots(1)
#ax4.set_yscale('log')
ax11.set_title(r'Homogeneous and Heterogeneous Surface Energies, $  \gamma^{Hom}_{SL}\,and\,\gamma^{Het}_{SL}   $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax11.plot(dTdr_v[1:36],gam_hom_v[1:36],'b-',linewidth = 1.5)
ax11.plot(dTdr_v[1:36],gam_het_v[1:36],'r-',linewidth = 1.5)
ax11.axhline(gamma_0,0,6475,color='r', linestyle='--',lw=1.5)


ax11.set_ylabel(r'$\mathrm{Surface\,Energy} \,\, ,\gamma_{SL}^{Hom,Het}\,\, \mathrm{[J.m^{-1}]} $',{'color': 'black', 'fontsize': 16})
ax11.set_xlabel(r'Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax11.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax11.get_yticklabels(), fontsize=16, fontweight="normal")

#plt.text(4060,0.7, r' $\mathrm{\frac{N}{V}\,\, = %g \,[modes.m^{-3}]} $' % (N_V), {'color': 'black', 'fontsize': 16})
#plt.text(4000,0.6, r' $\mathrm{\left ( \frac{N}{V} \right )_{i} = %g \,[conducting\, e^{-}.m^{-3}]} $' % (N_Vi), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.5, r' $\mathrm{E^{Fermi}}$ = %g [J] = %g [eV]' % (EF,EF*6.241509e18), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.4, r' $\mathrm{\omega _{D}}$ = %g [Hz]' % (wD), {'color': 'black', 'fontsize': 16})
plt.legend((r'Homogeneous',r'Heterogeneous',r'Literature reference value, $\gamma _{0}$ = %g $\mathrm{[J.m^{-1}]}$ ' % (gamma_0)),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )

plt.axis([0,7000,0,3])





plt.show()

 
 
 
    
    
#    gam = sigma_func(req, r, lambda_, mu, sigma, lambda0, mu0)




"""

# vectors
rchom_v = zeros(100)
grad_v  = zeros(100)

DTn = 195.000 # [K]
DTf = 0.05 # [K]

# homogenous nucleation
print(f'### Homogeneous Nucleation ###')
thet = pi
print(f'theta = {thet} radians')
fthet = ftheta(thet)
print(f'f(theta) = {fthet}')
# initial approach
dfthetdr = 0.0 #dfthetadr(thet)
print(f'df(theta)dr = {dfthetdr}')
#initial approach dgammadr = 0.0
dgamdr = 0.0
DTdx = 0.0

#GT_Eq_approach, delta_Sv_hom, delta_Ss_hom, delta_Sv_het, delta_Ss_het = gibbs_thomson_r(gamma_0, Delta_Sv_hom, DTf, dgamdr, fthet, dfthetdr, DTdx)
print("\nDTf = %g [K]\n" % (DTf))
GT_Eq_approach = gibbs_thomson_hom_r(gamma_0, Delta_Sv_hom, DTf, dgamdr, DTdx)
print("\nGibbs_Thomson initial approach = %g [m.K]\n" % (GT_Eq_approach))
#print("\nDelta_Sv_hom = %g [J/(K.m^3)]\n" % (delta_Sv_hom))
#print("\nDelta_Ss_hom = %g [J/(K.m^3)]\n" % (delta_Ss_hom))
#print("\nDelta_Sv_het = %g [J/(K.m^3)]\n" % (delta_Sv_het))
#print("\nDelta_Ss_het = %g [J/(K.m^3)]\n" % (delta_Ss_het))
rc_hom_Eq_approach = 2 * GT_Eq_approach / DTf
print("\nrc_Eq_approach = %g [m]\n" % (rc_hom_Eq_approach))
rc_hom_Eq_star = 0.99999 * rc_hom_Eq_approach
print("\nrc_Eq_star = %g [m]\n" % (rc_hom_Eq_star))

req = rc_hom_Eq_approach
r = rc_hom_Eq_approach * 0.8

sigma, surfstress =  sigma_func(req, r, lambda_, mu, sigma, lambda0, mu0)

print(f'\nsurface tension Eq. = {sigma_0} N/m\n')
print(f'\nsurface tension = {sigma} N/m\n')
print(f'\nsurface stress =  {surfstress} N/m\n')

r = rc_hom_Eq_approach

args = (rc_hom_Eq_approach, lambda_, mu, sigma_0, lambda0, mu0, DeltaH, rhos, T_Al, gamma_0, dTdx)

#rc_hom_non_Eq = fsolve(func_r, r, args)

rc_hom_non_Eq = brentq(func_r, 1e-12,rc_hom_Eq_approach, args)


print("\nrc_hom_non_Eq = %.6e [m]" % (rc_hom_non_Eq))



"""


"""
### Calcula o rEq considerando gradient de composição



for i in range(81):
    r = rc_hom_Eq_approach
    dTdx = i * 100
    
    args = (rc_hom_Eq_approach, lambda_, mu, sigma_0, lambda0, mu0, DeltaH, rhos, T_Al, gamma_0, dTdx)
    rc_hom_non_Eq = fsolve(func_r, r, args)[0]
    grad_v[i]  = dTdx
    rchom_v[i] = rc_hom_non_Eq
#    print("\nrc_hom_Eq_approach = %.6e [m]" % (rc_hom_Eq_approach))
#    print("\nrc_hom_non_Eq = %.6e [m] at dTdx = %.2f [K/m]" % (rc_hom_non_Eq,dTdx))
    print("\ndTdx = %.2f [K/m], rc_hom_non_Eq = %.6e [m], rc_hom_Eq_approach %.6e [m] " % (dTdx,rc_hom_non_Eq,rc_hom_Eq_approach))

"""


"""

sigm, surf_stress = sigma_func(rc_Eq_approach, rc_Eq_star, lambda_, mu, sigma_0, lambda0, mu0)
print("\n|sigma_0| = %g [N/m]\n " % (sigma_0))
print("\n|sigma| = %g [N/m]\n " % (sigm))
print("\ns(r = %g m) = %g [N/m] \n" % (rc_Eq_approach, surf_stress))
gamma = gamma_func_r(gamma_0, rc_Eq_star, rc_Eq_approach, surf_stress)
print("\ngamma_0 = %g [N/m]\n " % (gamma_0))
print("\ngamma = %g [J/m^2]\n" % (gamma))
dgamdr = dgammadr_func_r(rc_Eq_approach, rc_Eq_star, gamma_0, lambda_, mu, sigma_0, lambda0, mu0)
print("\ndgammadr = %g [J/m^3]\n" % (dgamdr))
#GT_Eq, delta_Sv_hom, delta_Ss_hom, delta_Sv_het, delta_Ss_het = gibbs_thomson_r(gamma, Delta_Sv_hom, DTf, dgamdr, fthet, dfthetdr, DTdx)
GT_Eq = gibbs_thomson_hom_r(gamma, Delta_Sv_hom, DTf, dgamdr, DTdx)
rc_Eq = 2 * GT_Eq / DTf


rc_hom = rc_Eq
print("\nrc_Eq = %g [m]\n" % (rc_Eq))

"""

"""

# initial guess for thet_het
thet = 0.01 * pi

def func(thet):
    
    fthet = ftheta(thet)
    
    dfthetdr = dfthetadr(thet)
    
    GT_Eq_het, delta_Sv_hom, delta_Ss_hom, delta_Sv_het, delta_Ss_het = gibbs_thomson_r(gamma, Delta_Sv_hom, DTf, dgamdr, fthet, dfthetdr, DTdx)
    
    rchet = rc_het(rc_Eq, thet)
    
    f  = 2 * GT_Eq_het / DTf - rchet
    
    return f
    
    
print(f'\nfunc(0.9*pi) = {func(0.9*pi)}\n')

thet_ = fsolve(func, thet)


print('\nthet_het = %2.14f radians\n' % (thet_))

print('\nthet_het = %2.14f degrees \n' % (thet_ * 180 / pi))



rc_heter = rc_het(rc_hom, thet_)

print(f'\nrc_hom = {rc_hom}m\n' )

print(f'\nrc_het = {rc_heter}m\n' )

fthet_final = ftheta(thet_)

print(f'\nf(thet_) = {fthet_final}\n' )

"""

"""
#Equilibrium Gibbs-Thomson
theta = 2*86 * pi / 180
print("\nNucleation angle = %g [radian]\n" % (theta))
Gibbs_Thomson, Delta_Sv = gibbs_thomson(gamma_0, 0.0, 1e-6, DeltaH, rhos, T_Al, theta)
dgamdr = dgammadr_func_r(1e-6, 1e-6, gamma_0, lambda_, mu, sigma, lambda0, mu0)
fthet  = ftheta(theta)
dfthetdr = dfthetadr(theta)
Gibbs_Thomson_new, Delta_Sv_hom, Delta_Ss_hom, Delta_Sv_het, Delta_Ss_het  = gibbs_thomson_r(gamma_0, Delta_Sv_hom, 1 , 0.0, fthet, 0, 0)
print("\nGibbs_Thomson = %g [m^3]\n" % (Gibbs_Thomson))
print("\nGibbs_Thomson_new = %g [m^3]\n" % (Gibbs_Thomson_new))
print("\nDelta_Sv_hom = %g [J/(K.m^3)]\n" % (Delta_Sv_hom))
print("\nDelta_Ss_hom = %g [J/(K.m^3)]\n" % (Delta_Ss_hom))
print("\nDelta_Sv_het = %g [J/(K.m^3)]\n" % (Delta_Sv_het))
print("\nDelta_Ss_het = %g [J/(K.m^3)]\n" % (Delta_Ss_het))

#DT_ = 0.09 #0.09 #0.105 #0.47 #0.52 #0.5
DT_ = 0.09 #0.209 #0.09 #0.105 #0.47 #0.52 #0.5
req = 2.0 * Gibbs_Thomson / DT_
print("\nrEq = %g [m^3]\n" % (req))
delta = 0.15 * req
r = req - delta
print("\nr_non_Eq = %g [m]\n" % (r))
print("\ndelta = %g [m]\n" % (delta))

sigm, surf_stress = sigma_func(req, r, lambda_, mu, sigma, lambda0, mu0)
gamma_v = gamma_func_r(gamma_0, r, req, surf_stress)

print("\n|sigma0| = %g [N/m]\n " % (sigma))
print("\n|sigma| = %g [N/m]\n " % (sigm))
print("\nsurface stress at r = %g m = %g [N/m] \n" % (r, surf_stress))
print("\ngamma_v = %g [J/m^2]\n" % (gamma_v))

Gibbs_Thomson_non_Eq, Delta_Sv_non_Eq = gibbs_thomson(gamma_v, delta, req, DeltaH, rhos, T_Al,theta)


#dgamdr = dgammadr_func_r(req, r, gamma_0, lambda_, mu, sigma, lambda0, mu0)
#fthet  = ftheta(theta)
#dfthetdr = dfthetadr(theta)
#delta_Sv = DeltaH * rhos / T_Al
#Gibbs_Thomson_non_Eq, Delta_Sv_non_Eq = gibbs_thomson_r(gamma_v, delta_Sv, DT_, dgamdr, fthet, dfthetdr, 0.0)

print("\nNon-Equilibrium Gibbs_Thomson = %g [K.m]\n" % (Gibbs_Thomson_non_Eq))
#DT_non_Eq = 2.0 * Gibbs_Thomson / r
DT_non_Eq = 2.0 * Gibbs_Thomson_non_Eq / r
print("\nDT_non_Eq = %g [K]\n" % (DT_non_Eq))

"""
