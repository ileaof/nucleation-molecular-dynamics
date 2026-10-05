from math import sin

r0 = 5.98498e-6
r1 = 5.95348e-6

gamma0 = 0.169 #0.311285
gamma1 = 0.326032

#theta0 = 3.05895
#theta1 = 3.05830

theta0 = 3.04652
theta1 = 3.04646
r0 = 4.56748E-06
r1 = 4.53598E-06


gamma = gamma0 

DP = 4 * gamma * ( r0**2 / r1 ** 3 + r0 ** 2 * sin(theta0)**2 / (r0 ** 3 * sin(theta1)**2) ) 

print(f'DP = {DP}')

DP2 = 4 * gamma * ( r0**3 / r1 ** 3 + r0 ** 3 * sin(theta0)**3 / (r0 ** 3 * sin(theta1)**3) ) 

print(f'DP2 = {DP2}')


DP3 = 4 * gamma * ( r0**2 / r1 ** 2 + r0 ** 2 * sin(theta0)**2 / (r0 ** 2 * sin(theta1)**2) ) 

print(f'DP3 = {DP3}')
