#!/bin/bash
import matplotlib.pyplot as plt
import mpl_textguard  # mpl>=3.10: skip off-chart text whose ft2font delta overflows int32
from numpy import zeros, sqrt, exp, log, pi, cos, sin, diff
from math import pi, acos
import csv
from midpoint    import midpoint
# only to calculate enthalpy
from scipy import integrate
from scipy.optimize import fsolve, brentq
from derivative_shim import derivative

# For Horizontal mean gradient calculation

temp_int = zeros(100)
cp_int = zeros(100)



"""
This script was updated on July 7th 2025 to write 1st-order variables on tablefull.csv file
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
dfgamdr_ana_het_v = zeros(n)
dTdr_v = zeros(n)
dTdr_het_v = zeros(n)
DT_v = zeros(n)
dDT_dr_v = zeros(n)
DGc_hom_expr_v = zeros(n)
DGc_het_expr_v = zeros(n)
DGc_hom_2nd_expr_v = zeros(n)
DGc_het_2nd_expr_v = zeros(n)
Delta_Sv_het_v = zeros(n)
Delta_Sv_hom_2nd_v = zeros(n)
Delta_Sv_het_2nd_v = zeros(n)
dDeltaSvhomdr_v = zeros(n)
dDeltaSvhomDTdr_v = zeros(n)
dDeltaSvhetdr_v = zeros(n)
dDeltaSvDTdr_v = zeros(n)
dDeltaSvhomDTfthetadr_v = zeros(n)
GT_hom_v = zeros(n)
GT_hom_v2 = zeros(n)
GT_hom_2nd_v = zeros(n)
r_hom_v = zeros(n)
r_hom_2nd_v = zeros(n)
r_het_v = zeros(n)
r_het_2nd_v = zeros(n)
rc_2nd_v = zeros(n)
r_spherical_cap = zeros(n)
theta_v = zeros(n)
theta_2nd_v = zeros(n)
dfthetadr_2nd_v = zeros(n)
termoA_v = zeros(n)
termoA_het_v = zeros(n)
GT_het_v = zeros(n)
GT_het_v2 = zeros(n)
GT_het_2nd_v = zeros(n)
V_hom_v = zeros(n)
V_het_v = zeros(n)
DoS_hom_v = zeros(n)
DoS_het_v = zeros(n)
gam_hom_v = zeros(n)
gam_het_v = zeros(n)
gam_hom_2nd_v = zeros(n)
gam_het_2nd_v = zeros(n)
gb_v = zeros(n)
gs_v = zeros(n)
gc_v = zeros(n)
sigma_v = zeros(n)
sigma_2nd_v = zeros(n)
sigma_het_v = zeros(n)
surface_stress_v = zeros(n)
surface_stress_2nd_v = zeros(n)
surface_stress_het_v = zeros(n)
surface_stress_abs_v = zeros(n)
surface_stress_het_abs_v = zeros(n)
crystal_order_hom_v = zeros(n)
crystal_order_het_v = zeros(n)
stat_hom_2nd_v = zeros(n)
stat_het_2nd_v = zeros(n)

# Total values
DeltaSV_Total_hom_v = zeros(n)
DeltaSV_Total_het_v = zeros(n)
DSv_hom_v = zeros(n)
dDSv_homdr_v = zeros(n)
DSv_het_v = zeros(n)
DSs_hom_v = zeros(n)
DSs_het_v = zeros(n)
DSc_v = zeros(n)
#

# By termoA Hom and Het
GT_hom_A_v = zeros(n)
GT_het_A_v = zeros(n)
#

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
TL = 922.608 #925.32 #890.65 
DeltaH = 293900.0 #335300 #234700 
rhol = 2433.79 #2378.33 #2493.98 
rhos = 2620.24 #2555.72 #2644.74
Delta_Sv_hom = -DeltaH * rhos / TL
#Equilibrium surface tension and surface stress
sigma_0 = 1.209 #1.09 #0.914 #0.821 #1.1 #0.821
sigma = sigma_0
gamma_0 = 0.169 #0.183
mu0 = 1.6
lambda0 = -5.0
rho0 = 2.549e-6
mu = 2.594e10
lambda_ = 5.034e10

rhet = 0.0

ad = 0.25e-9 # m Al atom jump from liquid to solid interface
CL_hom = 1.17594e28 # modes/m^3

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



def func_r_theta(r, req, lambda_, mu, sigma, lambda0, mu0, delta_H, rho, Tf, gamma_0, theta, dTdx):
    delta_Sv = delta_H * rho/Tf
#    y  = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*r   + 2.0 * (lambda0 + mu0)
#    y0 = (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)*req + 2.0 * (lambda0 + mu0)
    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    surf_stress  =  -2.0 * sigma  * (3.0*lambda_ + 2.0*mu)/ ( (3.0*lambda_ + 2.0*mu)*(1.0 + 2.0 * sigma)  ) * log ( y0 / y )
    surf_stress_theta = surf_stress * ( 2 ) / (1 - cos(theta))
    gamma_f  = gamma_0   / ((r/req) ** 2.0) - surf_stress_theta
    
    func_theta = gamma_f / (4.0 * pi * r * r * delta_Sv) - dTdx
    
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
    tau =  -sigma  -  2 * sigma  * log ( y0 / y )
    surf_stress  =  -2 * sigma  * log ( y0 / y )

    return tau, surf_stress


def sigma_func_theta(req, r, lambda_, mu, sigma, lambda0, mu0, theta):

    y  = (3*lambda_ + 2*mu)*r   + 2 * (lambda0 + mu0)
    y0 = (3*lambda_ + 2*mu)*req + 2 * (lambda0 + mu0)
    tau =  -sigma   -  2 * sigma * 2 / (1-cos(theta))  * log ( y0 / y ) 
    surf_stress  =  -2 * sigma * 2 / (1-cos(theta))  * log ( y0 / y )

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


def dgammadr_func_r(req, r, gamma_0, lambda_, mu, sigma, lambda0, mu0):
    alpha = 2 * (sigma + 2*(lambda0 + mu0)) / ( r*(3*lambda_ + 2*mu) )
    dgammadr = -2*gamma_0 * req ** 2 / (r ** 3) - 2*sigma / (r*(1+alpha))
    return dgammadr


def dgammadr_func_r_het(req, r, gamma_0, lambda_, mu, sigma, lambda0, mu0, theta):
    alpha = 2 * (sigma * 2 / (1-cos(theta)) + 2*(lambda0 + mu0)) / ( r*(3*lambda_ + 2*mu) )
    dgammadr = -2*gamma_0 * req ** 2 / (r ** 3) - 2*sigma * 2/(1-cos(theta)) / (r*(1+alpha))
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
#    print(f'integral = {sol}, error = {err}, rc_het = {rc_het}m')
    
    
    return rhet
    
###
###

def GB(Delta_Sv_hom, DT_):
    
      GB_ = Delta_Sv_hom * DT_ 

      return GB_


def GS(gammasl):

    GS_ = gammasl 

    return GS_


def GC(Delta_Sv_hom, DT_, rc_het_2nd, gammasl, theta_het_2nd):

    GC_ = gammasl / ftheta(theta_het_2nd) * dfthetadr(theta_het_2nd) #( 1/3 * pi * rc_het_2nd ** 3 * GB(Delta_Sv_hom, DT_)  + pi * rc_het_2nd ** 2 ) * ftheta(rc_het_2nd)

    return GC_





###
###


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
    GT = gammasl*ftheta / (( -delta_Sv_hom + delta_Ss_hom  )*ftheta - gammasl/DT_ * dfthetadr) - DTdx
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
    GT  = - gammasl / ( delta_Sv_hom + delta_Ss_hom )
    return GT


def gibbs_thomson_hom_r_2nd(delta_Sv_hom, DT, dgammadr, dDeltaSvhomdr, dDT_dr):
    """
    Returns 2nd-order homogeneous Gibbs-Thomson . For the equilibrium Gibbs-Thomson, simply make delta = 0, req = 1e6.
    Input data:
    delta_Sv_hom   - Volume entropy [J/(kg.K)];
    DT        - Supercooling [K];
    dgammadr   - derivative of gamma in relation to r [J/m^3];
    dDeltaSvhomdr - derivative of DeltaSv to r [J/K.m^4]
    dDT_dr - derivative of DT to r [K/m]
    """
    GT_hom_2nd = -3/4 * (delta_Sv_hom * DT + dgammadr)/(dDeltaSvhomdr + delta_Sv_hom/DT * dDT_dr)
    
    return GT_hom_2nd
 

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
    GT  = - gammasl*ftheta(theta) / (( delta_Sv_hom + delta_Ss_hom )*ftheta(theta) + gammasl/DT * dfthetadr(theta))
#    GT  = - gammasl / (( delta_Sv_hom + delta_Ss_hom ) + gammasl/DT * 1/ftheta(theta) * dfthetadr(theta))

    return GT


def gibbs_thomson_het_r_2nd(gammasl, delta_Sv_hom, DT, dgammadr, dDSvdr,dDTdr, theta):
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
    dDSvdr     - Derivative of bulk entropy in respect to r
    dDTdr      - Derivative of Undercooling in respect to r 
                 """
    
    #delta_Ss_hom = 1/DT * dgammadr
    
    
    #GThet2nd  = gammasl*ftheta(theta) / (( delta_Sv_hom - delta_Ss_hom )*ftheta(theta)/4 - gammasl/DT * dfthetadr(theta))    

    GThet2nd = -3 * ((delta_Sv_hom * DT + dgammadr )*ftheta(theta) + gammasl * dfthetadr(theta) ) / (4 * ((dDSvdr + delta_Sv_hom/DT * dDTdr) + delta_Sv_hom / ftheta(theta) * dfthetadr(theta)))
   
    #GThet2nd = -3/4 * (delta_Sv_hom * DT + dgammadr)/(dDeltaSvhomdr + delta_Sv_hom/DT * dDTdr)
    

    return GThet2nd



def DGc_hom_expr(gammasl,dgammadr, delta_Sv_hom, DT):
   
    DGchomexpr = 16/3 * pi * (gammasl ** 3) / (DT**2) * (delta_Sv_hom + 3/DT*dgammadr)/((delta_Sv_hom + 1/DT*dgammadr)**3)

    return DGchomexpr 


def DGc_hom_2nd_expr(gammasl,delta_Sv_hom, DT, dgammadr, dDSvdr, dDTdr):
   
    termo_1_Num = ( -(delta_Sv_hom * DT + dgammadr) ) ** 2
    termo_2a_Num = (delta_Sv_hom * (delta_Sv_hom * DT + dgammadr))
    termo_2b_Num = (2 * gammasl * (dDSvdr + delta_Sv_hom/DT * dDTdr))
    termo_Den    = ( dDSvdr + delta_Sv_hom/DT * dDTdr ) ** 3
    
    DGchom2ndexpr = -9 * pi / (2*DT**2) * (termo_1_Num * (termo_2a_Num - termo_2b_Num)) / termo_Den

    return DGchom2ndexpr 


def DGc_het_expr(gammasl,dgammadr, delta_Sv, DT, theta):
   
    DGchetexpr = 16/3 * pi * (gammasl ** 3) / (DT**2) * (delta_Sv - 3/DT*dgammadr - 3*gammasl/(DT*ftheta(theta))*dfthetadr(theta))/(((delta_Sv - 1/DT*dgammadr)-gammasl/(DT*ftheta(theta))*dfthetadr(theta)))**3 * ftheta(theta)/4
 
    return DGchetexpr 

def DGc_het_2nd_expr(gammasl,delta_Sv_hom, DT, dgammadr, dDSvdr, dDTdr, theta):
    
    termo_1_Num = ( -(delta_Sv_hom * DT - dgammadr) - gammasl/ftheta(theta)*dfthetadr(theta) ) ** 2
    termo_2a_Num = (delta_Sv_hom * ((delta_Sv_hom * DT - dgammadr) + gammasl/ftheta(theta)*dfthetadr(theta)))
    termo_2b_Num = (2 * gammasl * (dDSvdr + delta_Sv_hom/DT * dDTdr))
    termo_Den    = ( (dDSvdr + delta_Sv_hom/DT * dDTdr) + delta_Sv_hom/ftheta(theta) * dfthetadr(theta) ) ** 3
 
 
    DGchet2ndexpr = 9 * pi * ftheta(theta)/ (8*DT**2) * (termo_1_Num * (termo_2a_Num + termo_2b_Num)) / termo_Den

    return DGchet2ndexpr 

def DGc_2nd(r, DSv, DT, gammasl, theta):
    DGc =  (pi * r ** 2 * gammasl  + 1/3 * pi * r ** 3 * DSv * DT) * ftheta(theta)
    return DGc 

def Dself(T,Tm):

    """ 
    Self-diffusion coefficient of Al
    doi:10.1016/j.jnoncrysol.2006.03.049
    """
    DTm = 7.44e-9 # m^2/s
    p   = 0.10
    q   = -0.11
    D_self = DTm * (T/Tm)**(3/2) / ((1+p-p*T/Tm)*(1-q+q*T/Tm)**(2/3))
    return D_self


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
 
#    fu = dgamdr / (DSv - gam/(4 * pi * r * r * G)) - 8 * pi * r * G   
    fu = dgamdr / (gam/(4 * pi * r * r * G) + DSv) + 8 * pi * r * G  
    
    return fu


def f_DT(DT, gam, dgamdr, DSv, r ):
#    f_ = dgamdr / (DSv - 2*gam/(r * DT)) - DT
    f_ = dgamdr / (2*gam/(r * DT) + DSv) + DT

    return f_

def f_het(theta, rc_hom, termoa, termob, termoc ):
   
    fA= -termoc*ftheta(theta)/4 / (termoa*ftheta(theta)/4  + termob * dfthetadr(theta)) - rc_het(rc_hom, theta) 
    #fA= termoc*ftheta(theta) / (termoa*ftheta(theta)  - termob * dfthetadr(theta)) - rc_het(rc_hom, theta) 
    #fA= termoc / (termoa  - (termob/ftheta(theta)) * dfthetadr(theta)) - rc_het(rc_hom, theta)
    #fA= termoc / (termoa * ftheta(theta)/4 - termob*dfthetadr(theta)) - rc_het(rc_hom, theta)
    #fA= termoc / ( termoa  - termob/ftheta(theta)*dfthetadr(theta) ) - rc_het(rc_hom, theta)
    return fA

### TESTING ON MARCH 19TH THIS NEW FUNCTION FOR HETEROGENEOUS NUCLEATION ###
def f_het_var(theta, rc_hom, req, lambda_, mu, sigma_0, lambda0, mu0, gamma_0, delta_Sv_hom, DT):
    
    sigm_het, surf_stress_het = sigma_func_theta(req, rc_hom, lambda_, mu, sigma_0, lambda0, mu0, theta)
    gam_het = gamma_func_r(gamma_0, rc_hom, req, surf_stress_het)
    dgammadr_het = dgammadr_func_r_het(req, rc_hom, gamma_0, lambda_, mu, sigma, lambda0, mu0, theta)
    termoa = delta_Sv_hom - 1/DT * dgammadr_het
    termob = gam_het/DT
    termoc = 2 * gam_het / DT
    fA = termoc / (termoa  - (termob/ftheta(theta)) * dfthetadr(theta)) - rc_het(rc_hom, theta)

    return fA


    
#def f_het_2nd(theta, rc_hom_2nd, DSv, gam, dgamdr, dDSvdr, dDTdr, DT):
#
#    fB = -3*( (DSv*DT - dgamdr) - gam / ftheta(theta) * dfthetadr(theta) )  - rc_het(rc_hom_2nd, theta) * (2*DT*((dDSvdr + DSv/DT * dDTdr) + DSv/ftheta(theta) * dfthetadr(theta) ))
#
#    return fB

def f_het_2nd(theta, r_hom, DSv, DT, gam, dgamdr, dDSvDTdr):

    fc = - 3/2 * ((DSv*DT + dgamdr)*ftheta(theta)/4 + gam  * dfthetadr(theta)) / ( dDSvDTdr*ftheta(theta)/4 + DSv*DT*dfthetadr(theta) ) - rc_het(r_hom, theta)

    return fc



def f_het2(theta, gam, DSv, DT, dgamdr, GThom):

    fh2 = gam /((DSv - 1/DT*dgamdr) - gam/DT * 1/ftheta(theta) * dfthetadr(theta)) - GThom / (1 - GThom/DT * 1/ftheta(theta) * dfthetadr(theta))
 
    print(f'DSv = {DSv}')
    print(f'DT = {DT}')
    print(f'gam = {gam}')
    print(f'dgamdr = {dgamdr}')
    print(f'GThom = {GThom}')
    print(f'theta = {theta}')
    print(f'tA = {tA}')
    print(f'tB = {tB}')
    print(f'________________')
    return fh2 


def interpolation(Gxm1,Gx,Gxp1,Pxm1,Pxp1):
    
    Px = Pxm1 + (Pxp1-Pxm1) * (Gx - Gxm1) / (Gxp1 - Gxm1)
    
    return Px

#Phase mean molar composition or use a curve as a function of temperature provided
x_Cu = 0.00123132 #6.33256e-3
x_Si = 0.0 #0.00101715 #3.93038e-3
x_Mg = 0.00185518#0.00268834
x_C  = 0.0
x_Fe = 2.19594e-5 #2.8316e-5
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

GT_ref = 1.43e-7 #2.08e-7 # Phase
GL_EXP = 3708.2 #7902.38
GT_EXP = 3.3e-6
req = 6.3e-6 #5.96e-6 #3.79687e-06 #5.56e-6
dreq = 0.005 * req
dfdx = 0.0
#DT = 0.09
### Solid state calculation
N_V  = (1.0/(6.0*pi*pi))*(ThetaD_Phase*kB/(h_*vPhase))**3.0
N_Vi = Z_Phase * N_V
N_V_Al = x_Al * N_V
#N_V_Si = x_Si * N_V
N_V_Cu = x_Cu * N_V
N_V_Mg = x_Mg * N_V
N_V_Fe = x_Fe * N_V

P = 0.5
kF_Phase = (3 * pi ** 2 * Z_Phase * N_V)**(1/3)
kF_Al = (3 * pi ** 2 * Z_Al * N_V_Al)**(1/3)
#kF_Si = (3 * pi ** 2 * Z_Si * N_V_Si)**(1/3)
kF_Cu = (3 * pi ** 2 * Z_Cu * N_V_Cu)**(1/3)
kF_Mg = (3 * pi ** 2 * Z_Mg * N_V_Mg)**(1/3)
kF_Fe = (3 * pi ** 2 * Z_Fe * N_V_Fe)**(1/3)


vF_Phase = h_ * kF_Phase / me
vF_Al = h_ * kF_Al / me
#vF_Si = h_ * kF_Si / me
vF_Cu = h_ * kF_Cu / me
vF_Mg = h_ * kF_Mg / me
vF_Fe = h_ * kF_Fe / me


EF_Phase = 1/P * h_ ** 2 * kF_Phase**2 / (2*me)
EF_Al = 1/P * h_ ** 2 * kF_Al**2 / (2*me)
#Eg_Si = 1.17 /  6.241509e18
EF_Cu = 1/P * h_ ** 2 * kF_Cu**2 / (2*me)
EF_Mg = 1/P * h_ ** 2 * kF_Mg**2 / (2*me)
EF_Fe = 1/P * h_ ** 2 * kF_Fe**2 / (2*me)


#EF_Si = (Eg_Si / 2) * (N_V_Si)**(2/3) / (N_V)**(2/3)
T_F_Phase = EF_Phase / kB
T_F_Al = EF_Al / kB
#T_F_Si = EF_Si / kB
T_F_Cu = EF_Cu / kB
T_F_Mg = EF_Mg / kB
T_F_Fe = EF_Fe / kB


term_cve_Al = Z_Al*(ThetaD_Al**3.0)/T_F_Al
#term_cve_Si = Z_Si*(ThetaD_Si**3.0)/T_F_Si
term_cve_Cu = Z_Si*(ThetaD_Cu**3.0)/T_F_Cu
term_cve_Mg = Z_Mg*(ThetaD_Mg**3.0)/T_F_Mg
term_cve_Fe = Z_Fe*(ThetaD_Fe**3.0)/T_F_Fe




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

#mean unit vec
### Solid state


file = open('tabela_paper_MR_2025.csv','w')
filefull = open('tablefull_paper_MR_2025.csv','w')
# Calculo da TxL
#for i in range(countermax):
#    txsv[i] = gsv_pred[i] * vsv[i]
#    print(f'{tsv[i]},{xsv[i]},{vsv[i]},{gsv_pred[i]},{txsv[i]}')
#    file.write(f'{tsv[i]},{xsv[i]},{vsv[i]},{gsv_pred[i]},{txsv[i]}\n')
    
#file.close()

#file.write(f'i;r_hom;r_het;r_hom_2nd;r_het_2nd;sigma;surface_stress;gam_hom;gam_het;dfgamdr_ana;dfgamdr_ana_het;DT;dTdr;theta;theta_2nd;dDeltaSvhomDTdr;dDeltaSvhomDTfthetadr;GT_hom;GT_het;GT_hom_2nd;GT_het_2nd;DGc_hom_2nd_expr_v;DGc_het_2nd_expr \n')
#file.write(f'i;r_hom;r_het;r_hom_2nd;r_het_2nd;sigma;surface_stress;gam_hom;gam_het;dfgamdr_ana;dfgamdr_ana_het;DT;dTdr;theta;theta_2nd;dDeltaSvhomDTdr;dDeltaSvhomDTfthetadr;GT_hom;GT_het;GT_hom_2nd;GT_het_2nd;DGc_hom_2nd_expr_v;DGc_het_2nd_expr \n')
file.write(    f'i;dTdr;DT;r_hom_2nd;r_het_2nd;theta_2nd;sigma;surface_stress;gam_hom;GT_het;GT_het_2nd; \n')
filefull.write(f'i;dTdr;DT;r_hom_1st;r_het_1st;theta_1st;dfthetadr_1st;r_hom_2nd;r_het_2nd;theta_2nd;dfthetadr_2nd;sigma;surface_stress;gam_hom;dfgamdr_ana;GT_het;GT_het_2nd;DSv_hom;DSv_het;DSs_hom;DSs_het;DSc;GB;GS;GC;dDSv_homdr \n')


npoints = 75
npointsm1 = npoints - 1
for i in range(0,npoints):
    r = req - i * dreq   
    rc_v[i] = r
    r_hom_v[i] = r

    sigm, surfstress = sigma_func(req, r, lambda_, mu, sigma_0, lambda0, mu0)
    sigma_v[i] = sigm
    surface_stress_v[i] = surfstress
    surface_stress_abs_v[i] = abs(surfstress)
    gam = gamma_func_r(gamma_0, r, req, surfstress)
    gam_hom_v[i] = gam
    args = (req, lambda_, mu, sigma_0, lambda0, mu0)
    dfgamdr_num = derivative(fgam, r, dx=1e-10,args=args)
    dfgamdr_ana = dgammadr_func_r(req, r, gamma_0, lambda_, mu, sigma_0, lambda0, mu0)
    dfgamdr_num_v[i] = dfgamdr_num
    dfgamdr_ana_v[i] = dfgamdr_ana
    args =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
    dTdr = brentq(f,1,30000, args=args, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)   
    dTdr_v[i] = dTdr   
    args2 =   (gam, dfgamdr_ana, Delta_Sv_hom, r)
    DT = brentq(f_DT,0.0001,1, args=args2, xtol=2e-16, rtol=8.881784197001252e-16, maxiter=500, full_output=False, disp=True)
    DT_v[i] = DT      
    GT  = gibbs_thomson_hom_r(gam, Delta_Sv_hom, DT, dfgamdr_ana) 
    GT_hom_v[i]  = GT
 #   DGc_hom_expr_v[i] = DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)
    DGc_hom_expr_v[i] = ( pi * gam_hom_v[i] * gam_hom_v[i] * r_hom_v[i] ** 2 + pi/3 * Delta_Sv_hom * DT_v[i] * r_hom_v[i] ** 3 )  * 4   #DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)
    

## Second order-calculations
    if i > 0:
       DSv_hom_v[0] = Delta_Sv_hom
       dDT_dr_v[i] = (DT_v[i]-DT_v[i-1])/(r_hom_v[i]-r_hom_v[i-1]) #(DT_v[i]-DT_v[i-1])/dreq
       dDeltaSvhomdr = -3* (Delta_Sv_hom * DT_v[i] + dfgamdr_ana)/(2*r_hom_v[i]*DT) - Delta_Sv_hom  /DT_v[i] * dDT_dr_v[i]
       dDeltaSvhomDTdr_v[i] = dDeltaSvhomdr
 #      dDeltaSvhomDTdr_v[i] = -3/(2*r_hom_v[i]) * (Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i] ) #dDeltaSvhomdr
       r_hom_2nd_v[i] = r_hom_v[i] #(3/4) * r_hom_v[i] * ((Delta_Sv_hom*DT_v[i] + dfgamdr_ana_v[i]) ** 2) / (gam_hom_v[i]*(dDeltaSvhomdr * DT_v[i] + Delta_Sv_hom * dDT_dr_v[i]) )
#       DGc_hom_2nd_expr_v[i] = DGc_hom_2nd_expr(gam_hom_v[i], Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr_v[i], dDT_dr_v[i])
###       DGc_hom_2nd_expr_v[i] = ( pi * gam_hom_v[i] * gam_hom_v[i] * r_hom_2nd_v[i] ** 2 + pi/3 * Delta_Sv_hom * DT_v[i] * r_hom_2nd_v[i] ** 3 )  * 4   #DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)

#       GT_hom_2nd_v[i] = DT_v[i] * r_hom_2nd_v[i] / 2 #gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana, dDeltaSvhomdr, dDT_dr_v[i])  
#       GT_hom_2nd_v[i] =  -3/4 * (Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i]) / dDeltaSvhomDTdr_v[i] ##gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr, dDT_dr_v[i])  
       GT_hom_2nd_v[i] =  gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr, dDT_dr_v[i]) 

       dDeltaSvDTdr_v[i] = -3/2 * (Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i]) / r_hom_v[i]
       DGc_hom_2nd_expr_v[i] = DGc_2nd(r_hom_v[i], Delta_Sv_hom, DT_v[i], gam_hom_v[i], pi)
    
       DGc_CNT = -16/3 * pi * gam_hom_v[i] ** 3 * ftheta(theta_v[i])/ (Delta_Sv_hom**2 * (DT_v[i] * 1e2)**2)
       Tref = TL - DT_v[i]
       stat_hom_2nd_v[i] = (Dself(Tref,TL)/ad**2) * (4*pi*(r_hom_2nd_v[i]**2) / ad**2) * CL_hom * exp(DGc_hom_2nd_expr_v[i]/DGc_hom_2nd_expr_v[1])
       args4 = (r_hom_v[i], Delta_Sv_hom, DT_v[i], gam_hom_v[i],dfgamdr_ana_v[i],dDeltaSvDTdr_v[i])
       thet_2nd = brentq(f_het_2nd,0.01,pi, args=args4, xtol=2e-14, rtol=8.881784197001252e-14, maxiter=100, full_output=False, disp=True)
       theta_2nd_v[i] = thet_2nd
       r_het_2nd_v[i] = rc_het(r_hom_v[i],theta_2nd_v[i])
       dTdr_v[i] = DT_v[i] / (8*pi*r_hom_v[i])
       dTdr_het_v[i] = DT_v[i] /(4*pi*r_het_2nd_v[i]*(1-cos(theta_2nd_v[i])))
       DGc_het_2nd_expr_v[i] = DGc_2nd(r_hom_v[i], Delta_Sv_hom, DT_v[i], gam_hom_v[i], theta_2nd_v[i])
       

       stat_het_2nd_v[i] = (Dself(Tref,TL)/ad**2) * (2*pi*(r_het_2nd_v[i]**2) * (1-cos(theta_2nd_v[i])) / ad**2) * CL_hom * exp(DGc_het_2nd_expr_v[i]/DGc_het_2nd_expr_v[1])

       GT_hom_2nd_v[i] = 4 * pi * r_hom_v[i] ** 2 * dTdr_v[i]
       GT_het_2nd_v[i] = 2 * pi * r_hom_v[i] ** 2 * (1-cos(theta_2nd_v[i])) * dTdr_het_v[i]
       gam_hom_2nd_v[i] = gam_hom_v[i]
       gam_het_2nd_v[i] = gam_hom_v[i] * ftheta(theta_2nd_v[i])/4
       dDeltaSvhomDTfthetadr_v[i] = dDeltaSvhomDTdr_v[i] * ftheta(theta_2nd_v[i]) + Delta_Sv_hom * DT_v[i] * dfthetadr(theta_2nd_v[i])
       # Begin code change, added DeltaSV_Total_hom_v[i-1] +
       DeltaSV_Total_hom_v[i] =  dDeltaSvhomDTdr_v[i] * dreq / DT_v[i]
       DSv_hom_v[i] = DSv_hom_v[i-1] + dDeltaSvhomDTdr_v[i] * dreq / DT_v[i]
       # 
       dDSv_homdr_v[i] = (DSv_hom_v[i] - DSv_hom_v[i-1]) / (rc_v[i]-rc_v[i-1])
       # End code change
       dfthetadr_2nd_v[i] = dfthetadr(theta_2nd_v[i])


    TERMOA = Delta_Sv_hom + 1/DT_v[i] * dfgamdr_ana_v[i]  #Delta_Sv_hom + 1/DT * dfgamdr_ana
    TERMOB = gam_hom_v[i]  #gam/DT
    TERMOC = 2 * gam_hom_v[i] / DT_v[i] #2 * gam / DT
    args3 = (r, TERMOA, TERMOB, TERMOC)
    thet = brentq(f_het,0.01,pi, args=args3, xtol=2e-14, rtol=8.881784197001252e-14, maxiter=16, full_output=False, disp=True)
    theta_v[i] = thet 
    rhet = rc_het(r, thet)
    r_het_v[i] = rhet
    sigm_het, surfstress_het = sigma_func_theta(req, r, lambda_, mu, sigma_0, lambda0, mu0, thet)
    sigma_het_v[i] = sigm_het 
    surface_stress_het_v[i] = surfstress_het
    surface_stress_het_abs_v[i] = abs(surfstress_het)
 #   gam_het = gam_hom_v[i] * ftheta(theta_v[i]) #gamma_func_r(gamma_0, r, req, surfstress_het)
    gam_het = gam_hom_v[i] * ftheta(theta_v[i]) #gamma_func_r(gamma_0, r, req, surfstress_het)
 
    gam_het_v[i] = gam_het       
    # Begin change exchanging Delta_Sv_hom by eltaSV_Total_hom_v[i]
    Delta_Sv_het = Delta_Sv_hom * ftheta(theta_v[i])
    DSv_het_v[i] = DSv_hom_v[i] * ftheta(theta_2nd_v[i])/4.0
    # end change
    Delta_Sv_het_v[i] = Delta_Sv_het
    dfgamdr_ana_het = dfgamdr_ana_v[i] * ftheta(thet) +  gam_hom_v[i] * dfthetadr(thet) # To be check after
    dfgamdr_ana_het_v[i] = dfgamdr_ana_het # check after
##    dfgamdr_ana_het_v[i] = dgammadr_func_r_het(req, r, gamma_0, lambda_, mu, sigma_0, lambda0, mu0, thet)
    GT_het_v[i] =  gibbs_thomson_het_r(gam_het_v[i], Delta_Sv_het_v[i], DT_v[i],dfgamdr_ana_het_v[i],thet)  
    #GT_het_v[i] =  gibbs_thomson_het_r(gam_hom_v[i], Delta_Sv_hom, DT_v[i],dfgamdr_ana_v[i],thet)  
    DGc_het_expr_v[i] = ( pi * gam_het_v[i] * gam_het_v[i] * r_het_v[i] ** 2 + pi/3 * Delta_Sv_het_v[i] * DT_v[i] * r_het_v[i] ** 3 )  * ftheta(theta_v[i])
    DeltaSV_Total_het_v[i] = dDeltaSvhomDTfthetadr_v[i] * dreq/ ( ftheta(theta_v[i])*DT_v[i])

    DSs_hom_v[i] = dfgamdr_ana_v[i] / DT_v[i]
    DSs_het_v[i] = DSs_hom_v[i] * ftheta(theta_2nd_v[i])/4.0
    DSc_v[i] = gam_hom_v[i] / DT_v[i] * 1.0/ftheta(theta_2nd_v[i]) * dfthetadr(theta_2nd_v[i])

    gb_v[i] = GB(DSv_hom_v[i], DT_v[i])
    gs_v[i] = GS( gam_hom_v[i])
    gc_v[i] = GC(DSv_hom_v[i], DT_v[i],r_het_2nd_v[i], gam_hom_v[i],  theta_2nd_v[i])


    if i > 0:
       
 #      dDeltaSvhetdr_v[i] = dDeltaSvhomdr_v[i] * ftheta(theta_v[i]) + Delta_Sv_hom * dfthetadr(theta_v[i])    
 #      GT_het_2nd = gibbs_thomson_het_r_2nd(gam_hom_v[i],Delta_Sv_hom, DT_v[i], dfgamdr_ana, dDeltaSvhomdr, dDT_dr_v[i], theta_v[i])      
 #      r_het_2nd_v[i] = 3/4 * r_het_v[i] * ((Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i])*ftheta(theta_v[i]) + gam_hom_v[i] * dfthetadr(theta_v[i]))**2 / (gam_hom_v[i]*(dDeltaSvhomdr_v[i]*DT_v[i]+Delta_Sv_hom*dDT_dr_v[i])*ftheta(theta_v[i])+Delta_Sv_hom*DT_v[i]*dfthetadr(theta_v[i]))    
 #      #      GT_het_2nd = gibbs_thomson_het_r_2nd(gam_het_v[i],Delta_Sv_het_v[i], DT_v[i], dfgamdr_ana_het_v[i], dDeltaSvhetdr_v[i], dDT_dr_v[i], theta_v[i])      
 #      GT_het_2nd_v[i] = GT_het_2nd
 #      DGc_het_2nd_expr_v[i] = ( pi * gam_het_v[i] * gam_het_v[i] * r_het_2nd_v[i] ** 2 + pi/3 * Delta_Sv_het_v[i] * DT_v[i] * r_het_2nd_v[i] ** 3 )  * ftheta(theta_v[i])


 
    #print(f'r_hom_v[{i}] = {r_hom_v[i]:g},r_het_v[{i}] = {r_het_v[i]:g},r_hom_2nd_v[{i}] = {r_hom_2nd_v[i]:g}, sigma_v[{i}] = {sigma_v[i]:g},sigma_het_v[{i}] = {sigma_het_v[i]:g}, surface_stress_v[{i}]={surface_stress_v[i]:g}, surface_stress_het_v[{i}]={surface_stress_het_v[i]:g}, gam_hom_v[{i}]= {gam_hom_v[i]:g},gam_het_v[{i}]= {gam_het_v[i]:g}, dfgamdr_num_v[{i}] = {dfgamdr_num_v[i]:g}, dfgamdr_ana_v[{i}] = {dfgamdr_ana_v[i]:g},dfgamdr_ana_het_v[{i}]={dfgamdr_ana_het_v[i]:g},Delta_Sv[{i}] = {Delta_Sv_hom:g},  Delta_Sv_het[{i}]={Delta_Sv_het_v[i]:g},dTdr_v[{i}] = {dTdr_v[i]:g}, DT_v[{i}] = {DT_v[i]:g}, GT_hom_v[{i}] = {GT_hom_v[i]:g},GT_hom_2nd_v[{i}] = {GT_hom_2nd_v[i]:g},GT_het_v[{i}] = {GT_het_v[i]:g},GT_het_2nd_v[{i}] = {GT_het_2nd_v[i]:g}, DGc_hom_expr_v[{i}] = {DGc_hom_expr_v[i]:g},dDT_dr_v[{i}] = {dDT_dr_v[i]:g}, dDeltaSvhomdr_v[{i}] = {dDeltaSvhomdr_v[i]:g},  DGc_hom_2nd_expr_v[{i}]= {DGc_hom_2nd_expr_v[i]:g},theta_v[{i}]={theta_v[i]:g} \n')
     print(f'r_hom_v[{i}] = {r_hom_v[i]:g},theta_2nd_v[{i}]= {theta_2nd_v[i]:g},r_het_2nd_v[{i}] = {r_het_2nd_v[i]:g},dTdr_v[{i}] = {dTdr_v[i]:g},dTdr_het_v[{i}] = {dTdr_het_v[i]:g},GT_het_v[{i}]={GT_het_v[i]:g},stat_hom_2nd_v[{i}]={stat_hom_2nd_v[i]:g};{Delta_Sv_hom:g}; {DeltaSV_Total_hom_v[i]}')
#     file.write(f'{i};{r_hom_v[i]:g};{r_het_v[i]:g};{r_hom_2nd_v[i]:g};{r_het_2nd_v[i]:g};{sigma_v[i]:g};{surface_stress_v[i]:g};{gam_hom_v[i]:g};{gam_het_v[i]:g};{dfgamdr_ana_v[i]:g};{dfgamdr_ana_het_v[i]:g};{DT_v[i]:g};{dTdr_v[i]:g};{theta_v[i]:g};{theta_2nd_v[i]:g};{dDeltaSvhomDTdr_v[i]:g};{dDeltaSvhomDTfthetadr_v[i]:g};{GT_hom_v[i]:g};{GT_het_v[i]:g};{GT_hom_2nd_v[i]:g};{GT_het_2nd_v[i]:g};{DGc_hom_2nd_expr_v[i]:g};{DGc_het_2nd_expr_v[i]:.5e}\n')
     file.write(f'{i};{dTdr_v[i]:g};{DT_v[i]:g};{r_hom_2nd_v[i]:g};{r_het_2nd_v[i]:g};{theta_2nd_v[i]:g};{sigma_v[i]:g};{surface_stress_v[i]:g};{gam_hom_v[i]:g};{GT_het_v[i]:g};{GT_het_2nd_v[i]:g}\n')
     filefull.write(f'{i};{dTdr_v[i]:g};{DT_v[i]:g};{r_hom_v[i]:.5e};{r_het_v[i]:.5e};{theta_v[i]:g};{dfthetadr(theta_v[i]):g};{r_hom_2nd_v[i]:.5e};{r_het_2nd_v[i]:.5e};{theta_2nd_v[i]:g};{dfthetadr(theta_2nd_v[i]):g};{sigma_v[i]:g};{surface_stress_v[i]:g};{gam_hom_v[i]:g};{dfgamdr_ana_v[i]:g};{GT_het_v[i]:.5e};{GT_het_2nd_v[i]:.5e};{DSv_hom_v[i]:g};{DSv_het_v[i]:g};{DSs_hom_v[i]:g};{DSs_het_v[i]:g};{DSc_v[i]:g};{gb_v[i]:g};{gs_v[i]:g};{gc_v[i]:g};{dDSv_homdr_v[i]:.5e}\n')


     print(f'\n')

file.close()
#quit()
j = 45 #73 #45
jp1 = j + 1
dTdrj  = 3708.2 #7902.38 #4267.76 #4205 #4267.76
GL_EXP = dTdrj



r_hom_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],r_hom_v[j],r_hom_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m],r_hom_v[{j}] = {r_hom_v[j]:g} [m]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], r_hom_int[{(j+0.5):g}] = {r_hom_int:g} [m]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], r_hom_v[{jp1}] = {r_hom_v[jp1]:g} [m]')

r_hom_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],r_hom_2nd_v[j],r_hom_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m],r_hom_2nd_v[{j}] = {r_hom_2nd_v[j]:g} [m]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], r_hom_2nd_int[{(j+0.5):g}] = {r_hom_2nd_int:g} [m]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], r_hom_2nd_v[{jp1}] = {r_hom_2nd_v[jp1]:g} [m]')



r_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],r_het_v[j],r_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m],r_het_v[{j}] = {r_het_v[j]:g} [m]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], r_het_int[{(j+0.5):g}] = {r_het_int:g} [m]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], r_het_v[{jp1}] = {r_het_v[jp1]:g} [m]')


r_het_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],r_het_2nd_v[j],r_het_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m],r_het_2nd_v[{j}] = {r_het_2nd_v[j]:g} [m]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], r_het_2nd_int[{(j+0.5):g}] = {r_het_2nd_int:g} [m]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], r_het_2nd_v[{jp1}] = {r_het_2nd_v[jp1]:g} [m]')


theta_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],theta_v[j],theta_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], theta_v[{j}] = {theta_v[j]:g} [rad]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], theta_int[{(j+0.5):g}] = {theta_int:g} [rad]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], theta_v[{jp1}] = {theta_v[jp1]:g} [rad]')


theta_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],theta_2nd_v[j],theta_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], theta_2nd_v[{j}] = {theta_2nd_v[j]:g} [rad]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], theta_2nd_int[{(j+0.5):g}] = {theta_2nd_int:g} [rad]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], theta_2nd_v[{jp1}] = {theta_2nd_v[jp1]:g} [rad]')

#dfthetadr_2nd_v[i] 

dfthetadr_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dfthetadr_2nd_v[j],dfthetadr_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dfthetadr_2nd_v[{j}] = {dfthetadr_2nd_v[j]:g} [1/m]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dfthetadr_2nd_int[{(j+0.5):g}] = {dfthetadr_2nd_int:g} [1/m]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dfthetadr_2nd_v[{jp1}] = {dfthetadr_2nd_v[jp1]:g} [1/m]')


DT_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DT_v[j],DT_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DT_v[{j}] = {DT_v[j]:g} [K]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DT_int[{(j+0.5):g}] = {DT_int:g} [K]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DT_v[{jp1}] = {DT_v[jp1]:g} [K]')

GT_hom_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],GT_hom_2nd_v[j],GT_hom_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], GT_hom_2nd_v[{j}] = {GT_hom_2nd_v[j]:g} [m.K]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], GT_hom_2nd_int[{(j+0.5):g}] = {GT_hom_2nd_int:g} [m.K]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], GT_hom_2nd_v[{jp1}] = {GT_hom_2nd_v[jp1]:g} [m.K]')

GT_hom_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],GT_hom_v[j],GT_hom_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], GT_hom_v[{j}] = {GT_hom_v[j]:g} [m.K]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], GT_hom_int[{(j+0.5):g}] = {GT_hom_int:g} [m.K]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], GT_hom_v[{jp1}] = {GT_hom_v[jp1]:g} [m.K]')


GT_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],GT_het_v[j],GT_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], GT_het_v[{j}] = {GT_het_v[j]:g} [m.K]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], GT_het_int[{(j+0.5):g}] = {GT_het_int:g} [m.K]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], GT_het_v[{jp1}] = {GT_het_v[jp1]:g} [m.K]')



GT_het_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],GT_het_2nd_v[j],GT_het_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], GT_het_2nd_v[{j}] = {GT_het_2nd_v[j]:g} [m.K]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], GT_het_2nd_int[{(j+0.5):g}] = {GT_het_2nd_int:g} [m.K]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], GT_het_2nd_v[{jp1}] = {GT_het_2nd_v[jp1]:g} [m.K]')

sigma_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],sigma_v[j],sigma_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], sigma_v[{j}] = {-sigma_v[j]:g} [N.m^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], sigma_int[{(j+0.5):g}] = {sigma_int:g} [N.m^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], sigma_v[{jp1}] = {-sigma_v[jp1]:g} [N.m^-1]')



surface_stress_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],surface_stress_v[j],surface_stress_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], surface_stress_v[{j}] = {surface_stress_v[j]:g} [N.m^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], surface_stress_int[{(j+0.5):g}] = {surface_stress_int:g} [N.m^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], surface_stress_v[{jp1}] = {surface_stress_v[jp1]:g} [N.m^-1]')


gam_hom_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gam_hom_v[j],gam_hom_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gam_hom_v[{j}] = {gam_hom_v[j]:g} [J.m^-2]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gam_hom_v[{(j+0.5):g}] = {gam_hom_int:g} [J.m^-2]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gam_hom_v[{jp1}] = {gam_hom_v[jp1]:g} [J.m^-2]')


gam_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gam_het_v[j],gam_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gam_het_v[{j}] = {gam_het_v[j]:g} [J.m^-2]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gam_het_v[{(j+0.5):g}] = {gam_het_int:g} [J.m^-2]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gam_het_v[{jp1}] = {gam_het_v[jp1]:g} [J.m^-2]')

gam_het_2nd_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gam_het_2nd_v[j],gam_het_2nd_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gam_het_2nd_v[{j}] = {gam_het_2nd_v[j]:g} [J.m^-2]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gam_het_2nd_v[{(j+0.5):g}] = {gam_het_2nd_int:g} [J.m^-2]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gam_het_2nd_v[{jp1}] = {gam_het_2nd_v[jp1]:g} [J.m^-2]')


dfgamdr_ana_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dfgamdr_ana_v[j],dfgamdr_ana_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dfgamdr_hom_v[{j}] = {dfgamdr_ana_v[j]:g} [J.m^-3]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dfgamdr_hom_v[{(j+0.5):g}] = {dfgamdr_ana_int:g} [J.m^-3]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dfgamdr_hom_v[{jp1}] = {dfgamdr_ana_v[jp1]:g} [J.m^-3]')


dfgamdr_ana_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dfgamdr_ana_het_v[j],dfgamdr_ana_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dfgamdr_het_v[{j}] = {dfgamdr_ana_het_v[j]:g} [J.m^-3]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dfgamdr_het_v[{(j+0.5):g}] = {dfgamdr_ana_het_int:g} [J.m^-3]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dfgamdr_het_v[{jp1}] = {dfgamdr_ana_het_v[jp1]:g} [J.m^-3]')


dDeltaSvhomDTdr_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dDeltaSvhomDTdr_v[j],dDeltaSvhomDTdr_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dDeltaSvhomDTdr_v[{j}] = {dDeltaSvhomDTdr_v[j]:g} [J.m^-4]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dDeltaSvhomDTdr_int[{(j+0.5):g}] = {dDeltaSvhomDTdr_int:g} [J.m^-4]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dDeltaSvhomDTdr_v[{jp1}] = {dDeltaSvhomDTdr_v[jp1]:g} [J.m^-4]')


dDeltaSvhomDTfthetadr_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dDeltaSvhomDTfthetadr_v[j],dDeltaSvhomDTfthetadr_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dDeltaSvhomDTdr_v[{j}] = {dDeltaSvhomDTfthetadr_v[j]:g} [J.m^-4]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dDeltaSvhomDTfthetadr_int[{(j+0.5):g}] = {dDeltaSvhomDTfthetadr_int:g} [J.m^-4]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dDeltaSvhomDTfthetadr_v[{jp1}] = {dDeltaSvhomDTfthetadr_v[jp1]:g} [J.m^-4]')



#DGc_hom_2nd_expr_v


DGc_hom_2nd_expr_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DGc_hom_2nd_expr_v[j],DGc_hom_2nd_expr_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DGc_hom_2nd_expr_v[{j}] = {DGc_hom_2nd_expr_v[j]:g} [J]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dDGc_hom_2nd_expr_int[{(j+0.5):g}] = {DGc_hom_2nd_expr_int:g} [J]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DGc_hom_2nd_expr_v[{jp1}] = {DGc_hom_2nd_expr_v[jp1]:g} [J]')


DGc_het_2nd_expr_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DGc_het_2nd_expr_v[j],DGc_het_2nd_expr_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DGc_het_2nd_expr_v[{j}] = {DGc_het_2nd_expr_v[j]:g} [J]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dDGc_het_2nd_expr_int[{(j+0.5):g}] = {DGc_het_2nd_expr_int:g} [J]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DGc_het_2nd_expr_v[{jp1}] = {DGc_het_2nd_expr_v[jp1]:g} [J]')
print(f'\n')


DSv_hom_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DSv_hom_v[j],DSv_hom_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DSv_hom_v[{j}] = {DSv_hom_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DSv_hom_int[{(j+0.5):g}] = {DSv_hom_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DSv_hom_v[{jp1}] = {DSv_hom_v[jp1]:g} [J.m^-3.K^-1]')


dDSv_homdr_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],dDSv_homdr_v[j],dDSv_homdr_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], dDSv_homdr_v[{j}] = {dDSv_homdr_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], dDSv_homdr_int[{(j+0.5):g}] = {dDSv_homdr_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], dDSv_homdr_v[{jp1}] = {dDSv_homdr_v[jp1]:g} [J.m^-3.K^-1]')



DSv_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DSv_het_v[j],DSv_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DSv_het_v[{j}] = {DSv_het_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DSv_het_int[{(j+0.5):g}] = {DSv_het_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DSv_het_v[{jp1}] = {DSv_het_v[jp1]:g} [J.m^-3.K^-1]')


DSs_hom_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DSs_hom_v[j],DSs_hom_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DSs_hom_v[{j}] = {DSs_hom_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DSs_hom_int[{(j+0.5):g}] = {DSs_hom_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DSs_hom_v[{jp1}] = {DSs_hom_v[jp1]:g} [J.m^-3.K^-1]')

DSs_het_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DSs_het_v[j],DSs_het_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DSs_het_v[{j}] = {DSs_het_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DSs_het_int[{(j+0.5):g}] = {DSs_het_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DSs_het_v[{jp1}] = {DSs_het_v[jp1]:g} [J.m^-3.K^-1]')


DSc_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],DSc_v[j],DSc_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], DSc_v[{j}] = {DSc_v[j]:g} [J.m^-3.K^-1]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], DSc_int[{(j+0.5):g}] = {DSc_int:g} [J.m^-3.K^-1]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], DSc_v[{jp1}] = {DSc_v[jp1]:g} [J.m^-3.K^-1]')


gb_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gb_v[j],gb_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gb_v[{j}] = {gb_v[j]:g} [J/m^3]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gb_int[{(j+0.5):g}] = {gb_int:g} [J/m^3]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gb_v[{jp1}] = {gb_v[jp1]:g} [J/m^3]')
print(f'\n')


gs_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gs_v[j],gs_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gs_v[{j}] = {gs_v[j]:g} [J/m^2]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gs_int[{(j+0.5):g}] = {gs_int:g} [J/m^2]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gs_v[{jp1}] = {gs_v[jp1]:g} [J/m^2]')
print(f'\n')


gc_int = interpolation(dTdr_v[j],dTdrj,dTdr_v[jp1],gc_v[j],gc_v[jp1])
print(f'\n')
print(f'dTdr[{j}] = {dTdr_v[j]:g}[K/m], gc_v[{j}] = {gc_v[j]:g} [J]' )
print(f'dTdr[{(j+0.5):g}] = {dTdrj:g}[K/m], gc_int[{(j+0.5):g}] = {gc_int:g} [J]')
print(f'dTdr[{jp1}] = {dTdr_v[jp1]:g}[K/m], gs_v[{jp1}] = {gc_v[jp1]:g} [J]')
print(f'\n')


print(f'r_hom_int[{(j+0.5):g}] = {r_hom_int:g} [m], r_het_int[{(j+0.5):g}] = {r_het_int:g} [m],  theta_int[{(j+0.5):g}] = {theta_int:g} [rad], dTdr[{(j+0.5):g}] = {dTdrj:g} [K/m],  GT_hom_int[{(j+0.5):g}] = {GT_hom_2nd_int:g} [m.K], GT_het_int[{(j+0.5):g}] = {GT_het_2nd_int:g} [m.K], sigma_int[{(j+0.5):g}] = {sigma_int:g} [N.m^-1], surface_stress_int[{(j+0.5):g}] = {surface_stress_int:g} [N.m^-1], gam_hom_v[{(j+0.5):g}] = {gam_hom_int:g} [J.m^-2],  dfgamdr_hom_v[{(j+0.5):g}] = {dfgamdr_ana_int:g} [J.m^-3]')

filefull.write(f'{(j+0.5):g};{dTdrj:g};{DT_int:g};{r_hom_int:.5e};{r_het_int:.5e};{theta_int:g};{dfthetadr(theta_int):g};{r_hom_2nd_int:.5e};{r_het_2nd_int:.5e};{theta_2nd_int:g};{dfthetadr(theta_2nd_int):g};{sigma_int:g};{surface_stress_int:g};{gam_hom_int:g};{dfgamdr_ana_int:g};{GT_het_int:.5e};{GT_het_2nd_int:.5e};{DSv_hom_int:g};{DSv_het_int:g};{DSs_hom_int:g};{DSs_het_int:g};{DSc_int:g};{gb_int:g};{gs_int:g};{gc_int:g};{dDSv_homdr_int:.5e}\n')
filefull.close()


#quit()



fig0,ax0 = plt.subplots(1)
#ax0.set_yscale('log')
ax0.set_title(r'Surface Tension, $  \sigma_{SL}   $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax0.plot(1-rc_v[0:npointsm1]/req,sigma_v[0:npointsm1],'b-',linewidth = 1.5)
ax0.plot(1-rc_v[0:npointsm1]/req,sigma_het_v[0:npointsm1],'r-',linewidth = 1.5)
ax0.plot(1-rc_v[0:npointsm1]/req,surface_stress_v[0:npointsm1],'b--',linewidth = 1.5)
ax0.plot(1-rc_v[0:npointsm1]/req,surface_stress_het_v[0:npointsm1],'r--',linewidth = 1.5)
ax0.axhline(sigma_0,0,6475,color='c', linestyle='--',lw=1.5)
ax0.axvline((1-r_hom_int/req),0,1,color='g', linestyle='--',lw=1.5)
ax0.axhline(sigma_int,0,6475,color='k', linestyle='--',lw=1.5)
plt.text(0.025,-5.5, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax0.set_ylabel(r'Surface Tension, $\sigma_{SL}\,\, \mathrm{[N.m^{-1}]} $',{'color': 'black', 'fontsize': 16})
ax0.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax0.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax0.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax0.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Surface Tension, $\sigma _{SL}$', r'Surface Tension, $\sigma _{SL}^{Het}$',r'Surface Stress, $s^{Hom}$ ',r'Surfsace Stress, $s^{Het}$ ',r'Literature reference value, $\sigma _{0}$ = %g $\mathrm{[N.m^{-1}]}$ ' % (sigma_0), r'Heterogeneous radius for $\overline {\nabla T}_{EXP} \, =\,%g \, \mathrm{[K.m^{-1}]}$, $\mathit{\overline{r} _{C,Het}^{\,\,2nd-Order} }$ = %g $\mathrm{[m]}$ ' % (GL_EXP,r_het_2nd_int),r'Homogeneous Surface Tension, $\overline{\sigma} _{SL}$ = %g $\mathrm{[N.m^{-1}]}$ ' % (sigma_int)),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )


plt.axis([0,0.35,0,-10])




gam_EXP = 2.767
fig1,ax1 = plt.subplots(1)
#ax4.set_yscale('log')
ax1.set_title(r'Homogeneous and Heterogeneous Surface Energies, $  \gamma _{SL,Hom,Het}   $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
#ax1.plot(1-rc_v[0:npointsm1]/req,gam_hom_v[0:npointsm1],'b-',linewidth = 1.5)
#ax1.plot(1-rc_v[0:npointsm1]/req,gam_het_v[0:npointsm1],'r-',linewidth = 1.5)
#ax1.plot(1-rc_v[1:npointsm1]/req,gam_hom_2nd_v[1:npointsm1],'b+',linewidth = 1.5)
#ax1.plot(1-rc_v[1:npointsm1]/req,gam_het_2nd_v[1:npointsm1],'r*',linewidth = 1.5)


ax1.plot(dTdr_v[0:npointsm1],gam_hom_v[0:npointsm1],'b-',linewidth = 1.5)
ax1.plot(dTdr_v[0:npointsm1],gam_het_v[0:npointsm1],'r-',linewidth = 1.5)
ax1.plot(dTdr_v[1:npointsm1],gam_hom_2nd_v[1:npointsm1],'b+',linewidth = 1.5)
ax1.plot(dTdr_v[1:npointsm1],gam_het_2nd_v[1:npointsm1],'r*',linewidth = 1.5)


ax1.axhline(gamma_0,0,6475,color='r', linestyle='--',lw=1.5)
#ax1.axvline((1-r_hom_int/req),0,1,color='g', linestyle='--',lw=1.5)
ax1.axhline(gam_hom_int,0,6475,color='k', linestyle='--',lw=1.5)
ax1.axhline(gam_het_int,0,6475,color='c', linestyle='--',lw=1.5)
ax1.axhline(gam_het_2nd_int,0,6475,color='y', linestyle='--',lw=1.5)


ax1.axvline(GL_EXP,0,1,color='g', linestyle='--',lw=1.0)
#ax1.axhline(gam_EXP,0,6475,color='k', linestyle='--',lw=1.0)
plt.text(200,0.7, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})

ax1.set_ylabel(r'Surface Energy, $\gamma_{SL}\,\, \mathrm{[J.m^{-2}]} $',{'color': 'black', 'fontsize': 16})
#ax1.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})

ax1.set_xlabel(r'Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.text(4060,0.7, r' $\mathrm{\frac{N}{V}\,\, = %g \,[modes.m^{-3}]} $' % (N_V), {'color': 'black', 'fontsize': 16})
#plt.text(4000,0.6, r' $\mathrm{\left ( \frac{N}{V} \right )_{i} = %g \,[conducting\, e^{-}.m^{-3}]} $' % (N_Vi), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.5, r' $\mathrm{E^{Fermi}}$ = %g [J] = %g [eV]' % (EF,EF*6.241509e18), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.4, r' $\mathrm{\omega _{D}}$ = %g [Hz]' % (wD), {'color': 'black', 'fontsize': 16})
##plt.legend((r'Homogeneous',r'Heterogeneous',r'Literature reference value, $\gamma _{0}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gamma_0),r'Experimental $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP),r'Simulated $\overline{\gamma }_{SL}$ = %g $\mathrm{[J.m^{-2}]}$' %(gam_EXP)),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 11 )
plt.legend((r'Homogeneous, $\gamma _{SL,Hom}$',r'Heterogeneous, $\gamma _{SL,Het}$', r'Homogeneous 2nd-Order, $\gamma _{SL,Het} ^{\,\, 2nd-Order}$', r'Heterogeneous 2nd-Order, $\gamma _{SL,Het} ^{\,\, 2nd-Order}$',r'Literature reference value, $\gamma _{0}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gamma_0),r'Heterogeneous radius for $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$, $\overline{r} _{C,Hom} ^{\,\, 2nd-Order}$ = %g $\mathrm{[m]}$ ' % (GL_EXP,r_het_2nd_int),r'Homogeneous Surface Energy, $\overline{\gamma} _{SL}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_hom_int), r'Hoterogeneous Surface Energy, $\overline{\gamma} _{SL}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_het_int), r'Hoterogeneous Surface Energy 2nd-Order, $\overline{\gamma} _{SL}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_het_2nd_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 10 )


plt.axis([0,10000,0,2])

#plt.axis([0,0.35,0,2])



    

fig2,ax2 = plt.subplots(1)
ax2.set_title(r'Thermal gradient $\mathrm{\nabla T}$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax2.plot(1-rc_v[1:npointsm1]/req,dTdr_v[1:npointsm1],'b-',linewidth = 1.5)
ax2.plot(1-rc_v[1:npointsm1]/req,dTdr_het_v[1:npointsm1],'r+',linewidth = 1.5)
plt.text(0.05,7000, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax2.set_ylabel(r'Thermal Gradient, $\mathrm{\nabla T} \,\, \mathrm{ [K.m^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax2.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 16})
plt.setp(ax2.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax2.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Homogeneous 2nd-Order $\mathrm{  \hat n \, \cdot \, \nabla T = \frac {\Delta T} {8 \, \pi \, r_{C,Hom}}}$',r'Heterogeneous 2nd-Order, $\mathrm{  \hat n \, \cdot \, \nabla T = \frac {\Delta T} {4 \, \pi \, r_{C,Het} \left ( 1 - cos \left ( \theta \right ) \right ) }}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 ) #plt.axis([0,1,0,12000])
plt.axis([0,0.35,0,10000])


RB_reference = 0.2
fig3,ax3 = plt.subplots(1)
#ax1.set_yscale('log')
ax3.set_title(r'Undercooling $\mathit{\Delta T}$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax3.plot(1-rc_v[1:npointsm1]/req,DT_v[1:npointsm1],'b-',linewidth = 1.5)
plt.text(0.05,0.8, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax3.axhline(DT_int,0,1,color='b', linestyle='-.',lw=1.5)
ax3.axhline(RB_reference,0,1,color='r', linestyle='--',lw=1.5)
ax3.set_ylabel(r'$\mathit{\Delta T}\,\, \mathrm{ [K]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax3.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax3.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax3.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Undercooling, $\Delta T$',r'Undercooling for Mean Thermal Gradient $\overline{\nabla T}_{EXP}$ = %.0f $\mathrm{[K.m^{-1}]}$, $\overline{\Delta T}$ = %.4f [K]' % (dTdrj,DT_int),r'Literature Reference $\Delta T$ = %g $\mathrm{[K]}, $Rappaz and Boettinger (1999)' % (RB_reference) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 11 )
plt.axis([0,0.35,0,1.0])



fig4,ax4 = plt.subplots(1)
#ax1.set_yscale('log')
ax4.set_title(r'Homogeneous $\mathrm{\frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$ and Heterogeneous $\mathrm{\frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$ Derivatives of Bulk Free Energy $\Delta G _{V}$ of phase FCC_A1',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax4.plot(dTdr_v[1:npointsm1],dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
ax4.plot(dTdr_v[1:npointsm1],dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)


plt.text(4000,-0.4e12, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax4.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)

ax4.set_ylabel(r'$\mathrm{\frac {\partial G _{V}  } {\partial r}}$,  $\mathrm{ [J.m^{-4}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax4.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax4.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax4.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Homogeneous, $\frac {\partial G_{V, Hom}} {\partial r} \, = \, \mathit{\frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$', r'Heterogeneous, $\frac {\partial G_{V, Het}} {\partial r}  \, = \, \mathit{\frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.axis([0,10000,-1.2e12,0.25e12])


fig5,ax5 = plt.subplots(1)
#ax1.set_yscale('log')
ax5.set_title(r'Critical Free Energy $\mathrm{\Delta G_{C} }$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax5.plot(1-rc_v[1:npointsm1]/req,DGc_hom_2nd_expr_v[1:npointsm1],'b-',linewidth = 1.5)
ax5.plot(1-rc_v[1:npointsm1]/req,DGc_het_2nd_expr_v[1:npointsm1],'r+',linewidth = 1.5)

plt.text(0.025,9e-11, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax5.set_ylabel(r'Critical Free Energy, $\mathrm{\Delta G_{C} }$  $\mathrm{[J]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax5.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax5.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax5.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Homogeneneous 2nd-Order, $\Delta G_{C,Hom}$',r'Heterogeneous 2nd-Order, $\Delta G_{C,Het}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 16 )
plt.axis([0,0.35,3e-11,1e-10])


fig6,ax6 = plt.subplots(1)
#ax1.set_yscale('log')
ax6.set_title(r'2nd-Order Solution of Nucleation Angle, $\mathit{\theta }$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax6.plot(1-rc_v[1:npointsm1]/req,theta_2nd_v[1:npointsm1],'b-',linewidth = 1.5)
plt.text(0.1,3.06, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax6.set_ylabel(r'Nucleation Angle, $\mathit{\theta }$  $\mathrm{[rad]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax6.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax6.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax6.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
plt.axis([0,0.35,3.045,3.0650])


fig7,ax7 = plt.subplots(1)
#ax1.set_yscale('log')
ax7.set_title(r'1st-Order Solution of Nucleation Angle Approximation, $\mathrm{\theta }$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax7.plot(1-rc_v[1:npointsm1]/req,theta_v[1:npointsm1],'b-',linewidth = 1.5)
plt.text(0.05,0.9, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax7.set_ylabel(r'Nucleation Angle, $\mathit{\theta }$  $\mathrm{[rad]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax7.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax7.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax7.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
plt.axis([0,0.35,0.3,1.0])


"""
fig7,ax7 = plt.subplots(1)
#ax1.set_yscale('log')
ax7.set_title(r'Gibbs-Thomson, $\mathrm{\Gamma }$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
#ax7.plot(1-rc_v[1:npointsm1]/req,GT_hom_2nd_v[1:npointsm1],'b-',linewidth = 1.5)
#ax7.plot(1-rc_v[1:npointsm1]/req,GT_het_2nd_v[1:npointsm1],'r+',linewidth = 1.5)
ax7.plot(dTdr_v[1:npointsm1],GT_hom_2nd_v[1:npointsm1],'b-',linewidth = 1.5)
ax7.plot(dTdr_v[1:npointsm1],GT_het_2nd_v[1:npointsm1],'r+',linewidth = 1.5)
ax7.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T} \,\, \mathrm{ [K.m^{-1}]}$',{'color': 'black', 'fontsize': 16})

plt.text(0.05,0.8, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax7.set_ylabel(r'Gibbs-Thomson, $\mathrm{\Gamma }$,  $\mathrm{[K.m]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax7.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
plt.setp(ax7.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax7.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend(('Analytical','Numerical' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower left', fontsize = 16 )
#plt.axis([0,0.35,0,1.0])

"""

fig8,ax8 = plt.subplots(1)
#ax4.set_yscale('log')
ax8.set_title(r'Homogeneous and Heterogeneous Gibbs-Thomson, $\mathrm{\Gamma }$, phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax8.plot(dTdr_v[2:npointsm1],GT_het_v[2:npointsm1],'k*',linewidth = 2.5)
ax8.plot(dTdr_v[2:npointsm1],GT_het_2nd_v[2:npointsm1],'k-',linewidth = 2.5)
ax8.axhline(GT_ref,0,6475,color='r', linestyle='--',lw=1.5)
ax8.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)
#ax8.axhline(GT_hom_2nd_int,0,6475,color='b', linestyle='-.',lw=1.5)
ax8.axhline(GT_het_int,0,6475,color='r', linestyle='--',lw=1.0)
ax8.axhline(GT_het_2nd_int,0,6475,color='r', linestyle='-.',lw=1.0)

plt.text(200,1.35e-6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
plt.text(200,1.25e-6, r'$\mathrm{\Delta H_{FCC\_A1} \,\, = %g \,[J.kg]} $' % (DeltaH), {'color': 'black', 'fontsize': 12})
plt.text(200,1.15e-6, r'$\mathrm{T_{L} \,\, = \,\, = %g \,[K]} $' % (TL), {'color': 'black', 'fontsize': 12})
plt.text(200,1.05e-6, r'$\mathrm{\rho _{L} \,\, = \,\, = %g \,[kg.m^{-3}]} $' % (rhol), {'color': 'black', 'fontsize': 12})
plt.text(200,0.95e-6, r'$\mathrm{\rho _{FCC\_A1} \,\, = \,\, = %g \,[kg.m^{-3}]} $' % (rhos), {'color': 'black', 'fontsize': 12})

#ax4.plot(1-rc_v[1:91]/req,GT_hom_v2[1:91],'b-',linewidth = 1.5)
ax8.set_ylabel(r'$\mathrm{Gibbs-Thomson, \,\, \Gamma \,\, [m.K]}$',{'color': 'black', 'fontsize': 16})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax8.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax8.set_xlabel(r'Thermal Gradient, $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax1.legend((r'Al',r'Cu',r'Si'  ),loc="upper left",frameon=True,edgecolor='k',fontsize=14,shadow=True)
plt.setp(ax8.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax8.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'$\mathrm{\Gamma ^{hom} = \frac{\gamma _{SL}}{\Delta S_{V} - \frac{1}{\Delta T}\, \frac {\partial \gamma} {\partial r}} }$',r'$\mathrm{\Gamma ^{Hom} = \frac {r_{C}\,\Delta T}{2}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
#plt.legend((r'Homogeneous 1st-Order, $\mathrm{\Gamma ^{Hom} = \frac {\gamma _{SL}} {\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )}}$',r'Homogeneous 2nd-Order',r'Heterogeneous 1st-Order, $\mathrm{\Gamma ^{Het} = \frac {\gamma _{SL} } {\left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )  - \frac {\gamma _{SL} } {\Delta T} \frac{1} {f \left(\theta \right )} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref), r'Experimental $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP), r'Simulated $\overline{\Gamma}$ = %g [m.K]' %(GT_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
#plt.legend((r'Homogeneous 1st-Order, $\mathrm{\Gamma ^{Hom} = \frac {\gamma _{SL}} {\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )}}$',r'Homogeneous 2nd-Order, $\mathrm{\Gamma ^{Hom}\,=\,- \frac {3} {4} \, \frac {\left ( \Delta S_{V} \, \Delta T - \frac {\partial \gamma _{SL}} {\partial r}  \right ) } {\left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right)}  }$',r'Heterogeneous 1st-Order, $\mathrm{\Gamma ^{Het} = \frac {\gamma _{SL} } {\left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )  - \frac {\gamma _{SL} } {\Delta T} \frac{1} {f \left(\theta \right )} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Heterogeneous 2nd-Order, $\mathrm{\Gamma ^{Het}\,=\,- \frac {3} {4} \,  \frac { \left [ \left ( \Delta S_{V} \, \Delta T - \frac {\partial \gamma _{SL}} {\partial r}  \right ) - \frac {\gamma _{SL}} {f \left ( \theta \right ) } \frac {\partial f \left ( \theta \right )} {\partial r} \right ]} { \left [ \left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right) + \frac {\Delta S_{V}} {f \left ( \theta \right )  } \frac {f \left ( \theta \right ) } {\partial r} \right ] } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
#plt.legend((r'Homogeneous 1st-Order, $\mathrm{\Gamma ^{Hom} = \frac {\gamma _{SL}} {\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )}}$',r'Homogeneous 2nd-Order, $\mathrm{\Gamma ^{Hom}\,=\,- \frac {3} {4} \, \frac {\left ( \Delta S_{V} \, \Delta T - \frac {\partial \gamma _{SL}} {\partial r}  \right ) } {\left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right)}  }$',r'Heterogeneous 1st-Order, $\mathrm{\Gamma ^{Het} = \frac {\gamma _{SL} } {\left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )  - \frac {\gamma _{SL} } {\Delta T} \frac{1} {f \left(\theta \right )} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref), r'Expression for $\mathrm{}\Gamma _{Het}^{1st-order}}$' ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
#plt.legend((r'Homogeneous 1st-Order, $\mathrm{\Gamma ^{Hom} = \frac {\gamma _{SL}} {\left ( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )}}$',r'Homogeneous 2nd-Order, $\mathrm{\Gamma ^{Hom}\,=\,- \frac {3} {4} \, \frac {\left ( \Delta S_{V} \, \Delta T - \frac {\partial \gamma _{SL}} {\partial r}  \right ) } {\left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right)}  }$',r'Heterogeneous 1st-Order, $\mathrm{\Gamma ^{Het} = \frac {\gamma _{SL} } {\left( \Delta S_{V} - \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )  - \frac {\gamma _{SL} } {\Delta T} \frac{1} {f \left(\theta \right )} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref),r'Homogeneous from mean Gradient, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_hom_int) , r'Heterogeneous from mean Gradient, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
#plt.legend((r'Homogeneous 2nd-Order, $\mathrm{\Gamma ^{Hom}\,=\,- \frac {3} {4} \, \frac {\left ( \Delta S_{V} \, \Delta T - \frac {\partial \gamma _{SL}} {\partial r}  \right ) } {\left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right)}  }$','Heterogeneous 2nd-Order',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref),r'Homogeneous from mean Gradient, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_hom_2nd_int) , r'Heterogeneous from mean Gradient, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_het_2nd_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
plt.legend((r'Heterogeneous 1st-Order, $\mathrm{\Gamma _{Het} \, = \, - \, \frac {\gamma _{SL} } {\left( \Delta S_{V} + \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r} \right )  + \frac {\gamma _{SL} } {\Delta T} \frac{1} {f \left(\theta \right )} \frac {\partial f\left ( \theta \right ) } {\partial r} } } $',r'Heterogeneous 2nd-Order, $\mathrm{\Gamma _{Het}\,=\,- \frac {3} {4} \,  \frac { \left [ \left ( \Delta S_{V} \, \Delta T + \frac {\partial \gamma _{SL}} {\partial r}  \right ) + \frac {\gamma _{SL}} {f \left ( \theta \right ) } \frac {\partial f \left ( \theta \right )} {\partial r} \right ]} { \left [ \left ( \frac {\partial {\Delta S_{V}} } {\partial r} \, + \, \frac {\Delta S_{V}} {\Delta T} \, \frac {\partial {\Delta T} } {\partial r} \right) + \frac {\Delta S_{V}} {f \left ( \theta \right )  } \frac {f \left ( \theta \right ) } {\partial r} \right ] } } $',r'Literature reference value, $\Gamma$ = %.1e [m.K]' % (GT_ref), r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %.0f $\mathrm{[K.m^{-1}]}$' % (dTdrj) ,r'1nd-Order Heterogeneous Gibbs-Thomson, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_het_int), r'2nd-Order Heterogeneous Gibbs-Thomson, $\overline{\Gamma}$ = %.2e [m.K]' % (GT_het_2nd_int)) ,frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )


#plt.axis([0,1,0,12000])
#plt.axis([0,0.35,0,5e-6])
plt.axis([0,10000,0.5e-7,2.5e-6])
#plt.axis([0,900,8.5,13.0])




rhet_EXP = r_het_int
fig15,ax15 = plt.subplots(1)
#ax15.set_yscale('log')
ax15.set_title(r'Homogeneous and Heterogeneous Nucleation Radius, $  \mathrm{r_{C,Hom}}$ and $\mathrm{r_{C,Het}}$ of phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax15.plot(dTdr_v[1:npointsm1],rc_v[1:npointsm1],'b-',linewidth = 1.5)
ax15.plot(dTdr_v[1:npointsm1],r_hom_2nd_v[1:npointsm1],'b*',linewidth = 2.5)
ax15.plot(dTdr_v[1:npointsm1],r_het_v[1:npointsm1],'r-',linewidth = 1.5)
ax15.plot(dTdr_v[1:npointsm1],r_het_2nd_v[1:npointsm1],'r+',linewidth = 2.5)
ax15.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)
ax15.axhline(req,0,6475,color='g', linestyle='--',lw=1.5)
ax15.axhline(r_hom_int,0,6475,color='b', linestyle='--',lw=1.5)
ax15.axhline(r_het_int,0,6475,color='r', linestyle='--',lw=1.5)
ax15.axhline(r_hom_2nd_int,0,6475,color='b', linestyle='-.',lw=1.5)
ax15.axhline(r_het_2nd_int,0,6475,color='r', linestyle='-.',lw=1.5)

plt.text(5000,5.5e-6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax15.set_ylabel(r'Nucleation Radius, $\, \mathrm {r_{C}}\,\, \mathrm{[m]} $',{'color': 'black', 'fontsize': 16})
ax15.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax15.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax15.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.text(4060,0.7, r' $\mathrm{\frac{N}{V}\,\, = %g \,[modes.m^{-3}]} $' % (N_V), {'color': 'black', 'fontsize': 16})
#plt.text(4000,0.6, r' $\mathrm{\left ( \frac{N}{V} \right )_{i} = %g \,[conducting\, e^{-}.m^{-3}]} $' % (N_Vi), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.5, r' $\mathrm{E^{Fermi}}$ = %g [J] = %g [eV]' % (EF,EF*6.241509e18), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.4, r' $\mathrm{\omega _{D}}$ = %g [Hz]' % (wD), {'color': 'black', 'fontsize': 16})
#plt.legend((r'Homogeneous',r'Heterogeneous',r'Spherical cap',r'Literature reference value , $r^{Hom}_{eq}$ = %g $\mathrm{[m]}$ ' % (req),r'Experimental $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP),r'Simulated $\overline{r^{Het}_{C}}$ = %g $\mathrm{[m]}$' %(rhet_EXP)),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 11 )
#plt.legend((r'Homogeneous',r'Heterogeneous',r'Literature reference value , $r^{Hom}_{eq}$ = %g $\mathrm{[m]}$ ' % (req),r'Experimental $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP),r'Simulated $\overline{r^{Het}_{C}}$ = %g $\mathrm{[m]}$' %(rhet_EXP)),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 11 )
plt.legend((r'Homogeneous 1st-Order',r'Hemogeneous 2nd-Order',r'Heterogeneous 1st-Order',r'Heterogenous 2nd-Order',r'Experimental $\overline{\nabla T} _{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP),r'Literature reference value , $ r _{C,Hom,Eq}$ = %g $\mathrm{[m]}$ ' % (req) ,r'Homogeneous radius , $ r_{C,Hom}$ = %g $\mathrm{[m]}$ ' % (r_hom_int), r'Heterogeneous radius , $ r_{C,Het}$ = %g $\mathrm{[m]}$ ' % (r_het_int),r'Homogeneous 2nd-Order radius , $ r_{C,Hom}$ = %g $\mathrm{[m]}$ ' % (r_hom_2nd_int), r'Heterogeneous 2nd-Order radius , $r_{C,Het}$ = %g $\mathrm{[m]}$ ' % (r_het_2nd_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'center left', fontsize = 9 )

plt.axis([0,10000,0.25e-6,6.5e-6])


RB_reference = 0.2
fig16,ax16 = plt.subplots(1)
#ax1.set_yscale('log')
ax16.set_title(r'Undercooling $\mathit{\Delta T}$ phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax16.plot(dTdr_v[1:npointsm1],DT_v[1:npointsm1],'k-',linewidth = 1.5)
plt.text(800,0.90, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax16.axvline(dTdrj,0,1,color='k', linestyle='-.',lw=1.5)
ax16.axhline(DT_int,0,9000,color='b', linestyle='-.',lw=1.5)
ax16.axhline(RB_reference,0,9000,color='r', linestyle='--',lw=1.5)
ax16.set_ylabel(r'$\mathit{\Delta T}\,\, \mathrm{ [K]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
ax16.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax16.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax16.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Undercooling, $\Delta T$',r'Mean Thermal Gradient $\overline{\nabla T}_{EXP}$ = %.0f $\mathrm{[K.m^{-1}]}$' % (dTdrj), r'Mean Undercooling, $\overline{\Delta T}$ = %.4f [K] ' % (DT_int),r'Literature Reference $\Delta T$ = %g $\mathrm{[K]}, $Rappaz and Boettinger (1999)' % (RB_reference) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 11 )
plt.axis([0,10000,0,1.0])





#gam_EXP = 2.767
fig17,ax17 = plt.subplots(1)
#ax4.set_yscale('log')
ax17.set_title(r'Homogeneous and Heterogeneous Surface Energies Derivatives , $ \frac {\partial \gamma _{SL}} {\partial r}    $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
ax17.plot(dTdr_v[0:npointsm1],dfgamdr_ana_v[0:npointsm1],'b-',linewidth = 1.5)
ax17.plot(dTdr_v[0:npointsm1],dfgamdr_ana_het_v[0:npointsm1],'r-',linewidth = 1.5)
#ax17.plot(1-rc_v[1:npointsm1]/req,gam_hom_2nd_v[1:npointsm1],'b+',linewidth = 1.5)
#ax17.plot(1-rc_v[1:npointsm1]/req,gam_het_2nd_v[1:npointsm1],'r*',linewidth = 1.5)
#ax17.axhline(gamma_0,0,6475,color='r', linestyle='--',lw=1.5)
#ax17.axvline((1-r_hom_int/req),0,1,color='g', linestyle='--',lw=1.5)
ax17.axvline(dTdrj,0,1,color='k', linestyle='-.',lw=1.5)
ax17.axhline(dfgamdr_ana_int,0,6475,color='b', linestyle='--',lw=1.5)
ax17.axhline(dfgamdr_ana_het_int,0,6475,color='r', linestyle='--',lw=1.5)
#ax17.axhline(gam_het_2nd_int,0,6475,color='y', linestyle='--',lw=1.5)
#ax1.axhline(gam_EXP,0,6475,color='k', linestyle='--',lw=1.0)
plt.text(5000,-800000, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})

ax17.set_ylabel(r'Partial Dervitive of Surface Energy, $ \frac {\partial \gamma _{SL}} {\partial r}    \,\, \mathrm{[J.m^{-3}]} $',{'color': 'black', 'fontsize': 16})
ax17.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax17.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})

#ax1.set_xlabel(r'Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax17.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax17.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.text(4060,0.7, r' $\mathrm{\frac{N}{V}\,\, = %g \,[modes.m^{-3}]} $' % (N_V), {'color': 'black', 'fontsize': 16})
#plt.text(4000,0.6, r' $\mathrm{\left ( \frac{N}{V} \right )_{i} = %g \,[conducting\, e^{-}.m^{-3}]} $' % (N_Vi), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.5, r' $\mathrm{E^{Fermi}}$ = %g [J] = %g [eV]' % (EF,EF*6.241509e18), {'color': 'black', 'fontsize': 16})
#plt.text(4060,0.4, r' $\mathrm{\omega _{D}}$ = %g [Hz]' % (wD), {'color': 'black', 'fontsize': 16})
##plt.legend((r'Homogeneous',r'Heterogeneous',r'Literature reference value, $\gamma _{0}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gamma_0),r'Experimental $\overline{G}_L$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP),r'Simulated $\overline{\gamma }_{SL}$ = %g $\mathrm{[J.m^{-2}]}$' %(gam_EXP)),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 11 )
#plt.legend((r'Homogeneous, $\gamma _{SL,Hom}$',r'Heterogeneous, $\gamma _{SL,Het}$', r'Homogeneous 2nd-Order, $\gamma _{SL,Het} ^{\,\, 2nd-Order}$', r'Heterogeneous 2nd-Order, $\gamma _{SL,Het} ^{\,\, 2nd-Order}$',r'Literature reference value, $\gamma _{0}$ = %g $\mathrm{[J.m^{-2}]}$ ' % (gamma_0),r'Heterogeneous radius for $\Delta T\, =\,%g \, \mathrm{[K.m^{-1}]}$, |$\overline{r} _{C,Hom} ^{\,\, 2nd-Order}$| = %g $\mathrm{[m]}$ ' % (GL_EXP,r_het_2nd_int),r'Homogeneous Surface Energy, |$\overline{\gamma} _{SL}$| = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_hom_int), r'Hoterogeneous Surface Energy, |$\overline{\gamma} _{SL}$| = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_het_int), r'Hoterogeneous Surface Energy 2nd-Order, |$\overline{\gamma} _{SL}$| = %g $\mathrm{[J.m^{-2}]}$ ' % (gam_het_2nd_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 10 )
plt.legend((r'Homogeneous, $ \frac {\partial \gamma _{SL,Hom} } {\partial r}$', r'Heterogeneous, $ \frac {\partial \gamma _{SL,Het} } {\partial r}$', r'Mean Thermal Grandient $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$' % (GL_EXP), r'$\overline{\frac {\partial \gamma _{SL,Hom}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_int), r'$\overline{\frac {\partial \gamma _{SL,Het}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )

plt.axis([0,10000,0,-9e5])


func_v = zeros(10001)
x_v    = zeros(10001)
dx = 0.1
nv = 100
for j in range(nv):
#    x_v[j] = dx * j
#    func_v[j] = exp(x_v[j])    #exp(1e-11/(kB*TL))
#    print(f'x[{j}] = {x_v[j]}, func_v[{x_v[j]}]={func_v[j]}\n')
     Ttest = TL - 0.01 * j
     D_ = Dself(Ttest,TL)
     print(f'{Ttest:g}K, {D_:g} m^2/s')
     

    


fig18,ax18 = plt.subplots(1)
#ax1.set_yscale('log')
ax18.set_title(r'Nucleation Statistics phase FCC_A1',{'color': 'black', 'fontsize': 14})
ax18.plot(DT_v[1:npointsm1],stat_hom_2nd_v[1:npointsm1],'k-',linewidth = 1.5)
ax18.plot(DT_v[1:npointsm1],stat_het_2nd_v[1:npointsm1],'k-.',linewidth = 2.5)
#ax18.plot(1-rc_v[1:npointsm1]/req,stat_hom_2nd_v[1:npointsm1],'b-',linewidth = 1.5)
#ax18.plot(1-rc_v[1:npointsm1]/req,stat_het_2nd_v[1:npointsm1],'r+',linewidth = 1.5)


#plt.text(0.025,9e-11, r'Alloy Al-0.8wt%Si-0.6wt%Mg-0.2wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax18.set_ylabel(r'Nucleation Sites',{'color': 'black', 'fontsize': 18})
#ax18.set_ylabel(r'Critical Free Energy, $\mathrm{\Delta G_{C} }$  $\mathrm{[J]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax18.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax18.set_xlabel(r'Undercooling $\mathit{\Delta T}\,\, \mathrm{ [K]}$',{'color': 'black', 'fontsize': 18})
plt.text(0.5,1.05e50, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
plt.setp(ax18.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax18.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Homogeneneous 2nd-Order Statistics',r'Heterogeneous 2nd-Order Statistics' ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 16 )
plt.axis([0,1,0.2e50,1.75e50])



fig19,ax19 = plt.subplots(1)
#ax1.set_yscale('log')
ax19.set_title(r'Homogeneous $\mathrm{\frac {dr} {\Delta T} \frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$ and Heterogeneous $\mathrm{\frac {dr} {\Delta T \, f\left ( \theta \right )} \frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$ Total Entropy  $\Delta S_{T}$ of phase FCC_A1',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_hom_v[1:npointsm1],'b-',linewidth = 1.5)
ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_het_v[1:npointsm1],'r-',linewidth = 1.5)

plt.text(4000,-1.6e7, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax19.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)

ax19.set_ylabel(r'Total Entropy $\mathrm{ \Delta S_{T} }$,  $\mathrm{ [J.m^{-3}.K^{-1}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax19.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})

plt.setp(ax19.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax19.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Homogeneous, $\mathrm{\frac {dr} {\Delta T} \frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$', r'Heterogeneous, $\mathrm{\frac {dr} {\Delta T \, f\left ( \theta \right )} \frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
#plt.axis([0,10000,-3.0e7,1e4])


fig20,ax20 = plt.subplots(1)
#ax1.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax20.set_title(r'2nd-Order Heterogeneous Gibbs Bulk Free-Energy, $ \Delta G_B $',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax20.plot(dTdr_v[1:npointsm1],gb_v[1:npointsm1],'b-',linewidth = 1.5)
#ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_het_v[1:npointsm1],'r-',linewidth = 1.5)

plt.text(4000,-0.75e6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax20.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)

ax20.set_ylabel(r'Gibbs Bulk Free-Energy, $\mathrm{ \Delta G_{B} }$  $\mathrm{[J.m^{-3}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax20.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax20.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax20.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Gibbs Bulk Free-Energy, $\mathrm{\Delta G_{B}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 14 )
plt.axis([0,10000,-2.5e6,0])



fig21,ax21 = plt.subplots(1)
#ax1.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax21.set_title(r'2nd-Order Heterogeneous Gibbs Surface Free-Energy, $ \Delta G_S $',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax21.plot(dTdr_v[1:npointsm1],gs_v[1:npointsm1],'b-',linewidth = 1.5)
#ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_het_v[1:npointsm1],'r-',linewidth = 1.5)
plt.text(4000,1.7, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax21.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)
ax21.set_ylabel(r'Gibbs Surface Free-Energy, $\mathrm{ \Delta G_{S} }$  $\mathrm{[J.m^{-2}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax21.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax21.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax21.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'Homogeneous, $\mathrm{\frac {dr} {\Delta T} \frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$', r'Heterogeneous, $\mathrm{\frac {dr} {\Delta T \, f\left ( \theta \right )} \frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.legend((r'Gibbs Surface Free-Energy, $\mathrm{\Delta G_{S}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.axis([0,10000,0,2.0])


fig22,ax22 = plt.subplots(1)
#ax1.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax22.set_title(r'2nd-Order Heterogeneous Bulk Entropy, $ \Delta S_{V} $',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax22.plot(dTdr_v[1:npointsm1],DSv_het_v[1:npointsm1],'b-',linewidth = 1.5)
#ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_het_v[1:npointsm1],'r-',linewidth = 1.5)
plt.text(4000,-1.5e6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax22.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)
ax22.set_ylabel(r'Bulk Entropy, $\mathrm{ \Delta S_{V} }$  $\mathrm{[J.m^{-3}.K^{-1}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax22.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax22.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax22.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'Homogeneous, $\mathrm{\frac {dr} {\Delta T} \frac {\partial \left ( \Delta S_{V} \, \Delta T \right ) } {\partial r}}$', r'Heterogeneous, $\mathrm{\frac {dr} {\Delta T \, f\left ( \theta \right )} \frac {\partial \left ( \Delta S_{V} \, \Delta T \, f \left ( \theta \right ) \right ) } {\partial r}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.legend((r'Bulk Entropy, $\mathrm{\Delta S_{V}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 14 )
plt.axis([0,10000,-2.6e6,-1e6])


fig23,ax23 = plt.subplots(1)
#ax4.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax23.set_title(r'Surface Entropy , $ \Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r}    $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
#ax23.plot(dTdr_v[1:npointsm1],DSs_hom_v[1:npointsm1],'b.-',linewidth = 1.5)
ax23.plot(dTdr_v[1:npointsm1],DSs_het_v[1:npointsm1],'b-',linewidth = 1.5)
ax23.axvline(dTdrj,0,1,color='k', linestyle='-.',lw=1.5)
plt.text(5000,-2e6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax23.set_ylabel(r'Surface Entropy, $ \Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL}} {\partial r}    \,\, \mathrm{[J.m^{-3}.K^{-1}]} $',{'color': 'black', 'fontsize': 16})
ax23.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax17.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.set_xlabel(r'Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax23.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax23.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'Homogeneous, $ \Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL,Hom} } {\partial r}$', r'Heterogeneous, $\Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL,Het} } {\partial r}$', r'Mean Thermal Grandient $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$' % (GL_EXP), r'$\overline{\frac {\partial \gamma _{SL,Hom}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_int), r'$\overline{\frac {\partial \gamma _{SL,Het}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
plt.legend((r'Heterogeneous, $\Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL,Het} } {\partial r}$', r'Mean Thermal Grandient $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$' % (GL_EXP), r'$\overline{\frac {\partial \gamma _{SL,Hom}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_int), r'$\overline{\frac {\partial \gamma _{SL,Het}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.axis([0,10000,-6.5e6,-0.5e6])



fig24,ax24 = plt.subplots(1)
#ax4.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax24.set_title(r'Configurational Entropy, $ \Delta S_{C} = \frac {\gamma _{SL}} {\Delta T} \, \frac {1} {f(\theta)} \, \frac {\partial f(\theta) } {\partial r}    $ phase FCC_A1 ',{'color': 'black', 'fontsize': 14})
#ax23.plot(dTdr_v[1:npointsm1],DSs_hom_v[1:npointsm1],'b.-',linewidth = 1.5)
ax24.plot(dTdr_v[1:npointsm1],DSc_v[1:npointsm1],'b-',linewidth = 1.5)
ax24.axvline(dTdrj,0,1,color='k', linestyle='-.',lw=1.5)
plt.text(5000,-2e6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax24.set_ylabel(r'Configurational Entropy, $ \Delta S_{C} = \frac {\gamma _{SL}} {\Delta T} \, \frac {1} {f(\theta)} \, \frac {\partial f(\theta)} {\partial r}    \,\, \mathrm{[J.m^{-3}.K^{-1}]} $',{'color': 'black', 'fontsize': 16})
ax24.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
#ax17.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
#ax1.set_xlabel(r'Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax24.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax24.get_yticklabels(), fontsize=16, fontweight="normal")
#plt.legend((r'Homogeneous, $ \Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL,Hom} } {\partial r}$', r'Heterogeneous, $\Delta S_{S} = \frac {1} {\Delta T} \frac {\partial \gamma _{SL,Het} } {\partial r}$', r'Mean Thermal Grandient $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$' % (GL_EXP), r'$\overline{\frac {\partial \gamma _{SL,Hom}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_int), r'$\overline{\frac {\partial \gamma _{SL,Het}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 11 )
plt.legend((r'$\Delta S_{C} = \frac {\gamma _{SL}} {\Delta T} \, \frac {1} {f(\theta)} \, \frac {\partial f(\theta) } {\partial r}$', r'Mean Thermal Grandient $\overline{\nabla T}_{EXP}\, =\,%g \, \mathrm{[K.m^{-1}]}$' % (GL_EXP), r'$\overline{\frac {\partial \gamma _{SL,Hom}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_int), r'$\overline{\frac {\partial \gamma _{SL,Het}} {\partial r}}$ = %.3e $\mathrm{[J.m^{-3}]}$' % (dfgamdr_ana_het_int) ),frameon=True,edgecolor='k',shadow=True, loc= 'lower right', fontsize = 14 )
plt.axis([0,10000,-11,-6.0])


fig25,ax25 = plt.subplots(1)
#ax1.set_yscale('log')
plt.grid(True, linestyle='--', color='darkgray', linewidth=0.5)
ax25.set_title(r'2nd-Order Heterogeneous Gibbs Configurational Free-Energy, $ \Delta G_B $',{'color': 'black', 'fontsize': 12})
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTdr_v[1:npointsm1],'b-',linewidth = 1.5)
#ax4.plot(1-rc_v[1:npointsm1]/req,dDeltaSvhomDTfthetadr_v[1:npointsm1],'r-',linewidth = 1.5)
ax25.plot(dTdr_v[1:npointsm1],gc_v[1:npointsm1],'b-',linewidth = 1.5)
#ax19.plot(dTdr_v[1:npointsm1],DeltaSV_Total_het_v[1:npointsm1],'r-',linewidth = 1.5)

plt.text(4000,-0.75e6, r'Alloy Al-6wt%Cu-4wt%Si-0.11wt%Fe' , {'color': 'black','weight':'bold', 'fontsize': 12})
ax25.axvline(GL_EXP,0,1,color='k', linestyle='--',lw=1.5)

ax25.set_ylabel(r'Gibbs Configurational Free-Energy, $\mathrm{ \Delta G_{B} }$  $\mathrm{[J.m^{-3}]}$',{'color': 'black', 'fontsize': 18})
#ax1.set_ylabel(r'Chemical potential, $\mu\,\, \mathrm{[J.mol^{-1}]}$',{'color': 'black', 'fontsize': 16})
#ax4.set_xlabel(r'$\mathrm{1 - \frac {r}{r_{ref}}}$',{'color': 'black', 'fontsize': 18})
ax25.set_xlabel(r'Macroscopic Thermal Gradient,  $\mathrm{\nabla T \,\, [K.m^{-1}] }$',{'color': 'black', 'fontsize': 18})
plt.setp(ax25.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax25.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Gibbs Configurational Free-Energy, $\mathrm{\Delta G_{B}}$', r'Mean Thermal Gradient, $\overline{\nabla T}_{EXP}$ = %g $\mathrm{[K.m^{-1}]}$' % (GL_EXP) ),frameon=True,edgecolor='k',shadow=True, loc= 'upper right', fontsize = 14 )
plt.axis([0,10000,-7,0])






plt.show()



    #print(f'r_hom_v[{i}] = {r_hom_v[i]:g},r_het_v[{i}] = {r_het_v[i]:g},r_hom_2nd_v[{i}] = {r_hom_2nd_v[i]:g}, sigma_v[{i}] = {sigma_v[i]:g},sigma_het_v[{i}] = {sigma_het_v[i]:g}, surface_stress_v[{i}]={surface_stress_v[i]:g}, surface_stress_het_v[{i}]={surface_stress_het_v[i]:g}, gam_hom_v[{i}]= {gam_hom_v[i]:g},gam_het_v[{i}]= {gam_het_v[i]:g}, dfgamdr_num_v[{i}] = {dfgamdr_num_v[i]:g}, dfgamdr_ana_v[{i}] = {dfgamdr_ana_v[i]:g},dfgamdr_ana_het_v[{i}]={dfgamdr_ana_het_v[i]:g},Delta_Sv[{i}] = {Delta_Sv_hom:g},  Delta_Sv_het[{i}]={Delta_Sv_het_v[i]:g},dTdr_v[{i}] = {dTdr_v[i]:g}, DT_v[{i}] = {DT_v[i]:g}, GT_hom_v[{i}] = {GT_hom_v[i]:g},GT_hom_2nd_v[{i}] = {GT_hom_2nd_v[i]:g},GT_het_v[{i}] = {GT_het_v[i]:g},GT_het_2nd_v[{i}] = {GT_het_2nd_v[i]:g}, DGc_hom_expr_v[{i}] = {DGc_hom_expr_v[i]:g},dDT_dr_v[{i}] = {dDT_dr_v[i]:g}, dDeltaSvhomdr_v[{i}] = {dDeltaSvhomdr_v[i]:g},  DGc_hom_2nd_expr_v[{i}]= {DGc_hom_2nd_expr_v[i]:g},theta_v[{i}]={theta_v[i]:g} \n')

    
"""   
    
    
    GT  = gibbs_thomson_hom_r(gam, Delta_Sv_hom, DT, dfgamdr_ana) 
    GT_hom_v[i]  = GT
 #   DGc_hom_expr_v[i] = DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)
    DGc_hom_expr_v[i] = ( pi * gam_hom_v[i] * gam_hom_v[i] * r_hom_v[i] ** 2 + pi/3 * Delta_Sv_hom * DT_v[i] * r_hom_v[i] ** 3 )  * 4   #DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)



## Second order-calculations
    if i > 0:
       dDT_dr_v[i] = (DT_v[i]-DT_v[i-1])/(r_hom_v[i]-r_hom_v[i-1]) #(DT_v[i]-DT_v[i-1])/dreq
       dDeltaSvhomdr = -3* (Delta_Sv_hom * DT_v[i] + dfgamdr_ana)/(2*r_hom_v[i]*DT) - Delta_Sv_hom  /DT_v[i] * dDT_dr_v[i]
       dDeltaSvhomDTdr_v[i] = dDeltaSvhomdr
 #      dDeltaSvhomDTdr_v[i] = -3/(2*r_hom_v[i]) * (Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i] ) #dDeltaSvhomdr
       r_hom_2nd_v[i] = r_hom_v[i] #(3/4) * r_hom_v[i] * ((Delta_Sv_hom*DT_v[i] + dfgamdr_ana_v[i]) ** 2) / (gam_hom_v[i]*(dDeltaSvhomdr * DT_v[i] + Delta_Sv_hom * dDT_dr_v[i]) )
#       DGc_hom_2nd_expr_v[i] = DGc_hom_2nd_expr(gam_hom_v[i], Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr_v[i], dDT_dr_v[i])
       DGc_hom_2nd_expr_v[i] = ( pi * gam_hom_v[i] * gam_hom_v[i] * r_hom_2nd_v[i] ** 2 + pi/3 * Delta_Sv_hom * DT_v[i] * r_hom_2nd_v[i] ** 3 )  * 4   #DGc_hom_expr(gam,dfgamdr_ana, Delta_Sv_hom, DT)

#       GT_hom_2nd_v[i] = DT_v[i] * r_hom_2nd_v[i] / 2 #gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana, dDeltaSvhomdr, dDT_dr_v[i])  
#       GT_hom_2nd_v[i] =  -3/4 * (Delta_Sv_hom * DT_v[i] + dfgamdr_ana_v[i]) / dDeltaSvhomDTdr_v[i] ##gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr, dDT_dr_v[i])  
       GT_hom_2nd_v[i] =  gibbs_thomson_hom_r_2nd (Delta_Sv_hom, DT_v[i], dfgamdr_ana_v[i], dDeltaSvhomdr, dDT_dr_v[i])  


 
 """