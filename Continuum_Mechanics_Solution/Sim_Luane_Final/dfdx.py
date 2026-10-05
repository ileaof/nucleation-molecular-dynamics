import matplotlib.pyplot as plt
from sympy import diff, Symbol, Function, sin
from math import pi
from derivative_shim import derivative
from numpy import zeros

nv  = 21
r_v = zeros(nv)
yline_ana_v = zeros(nv)
yline_num_v = zeros(nv)


r0 = Symbol('r0')
r  = Symbol('r')
theta = Symbol('theta')
theta0 = Symbol('theta0')

metrics = Function('metrics')

metrics = 1 + r ** 2 / r0 ** 2 + r**2 * sin(theta) ** 2 / (r0**2 * sin(theta0) ** 2) 

f = 3 / metrics 

y = diff(f,r)

print(f"f'= {y}")

df = 3*(-2*r*sin(theta)**2/(r0**2*sin(theta0)**2) - 2*r/r0**2)/(r**2*sin(theta)**2/(r0**2*sin(theta0)**2) + r**2/r0**2 + 1)**2


r0 = 1.0
r  = 0.99
theta  = pi/2
theta0 = pi/2


def y(r, r0, theta, theta0):
    from math import sin

    metrics = 1 + r ** 2 / r0 ** 2 + r**2 * sin(theta) ** 2 / (r0**2 * sin(theta0) ** 2) 
    
    y_ = 3 / metrics 

    return y_


def yline(r, r0, theta, theta0):
    from math import sin

    yline_ = 3*(-2*r*sin(theta)**2/(r0**2*sin(theta0)**2) - 2*r/r0**2)/(r**2*sin(theta)**2/(r0**2*sin(theta0)**2) + r**2/r0**2 + 1)**2
    
    return yline_

print(f'yline(0.96) = {yline(r,r0,theta,theta0)}')

args = (r0, theta, theta0)
   
yline_num = derivative(y, r, dx=1e-10,args=args)

print(f'yline_numeric({r:g}) = {yline_num}')

for i in range(20):
    dr = 0.0025 * r0
    r  = r0 - i * dr
    yline_ana = yline(r, r0, theta, theta0)
    args = (r0, theta, theta0)   
    yline_num = derivative(y, r, dx=1e-10,args=args)
    r_v[i] = r
    yline_ana_v[i] = yline_ana
    yline_num_v[i] = yline_num

    print(f'{i},  {r},     {yline_ana:.6f},      {yline_num:6f}')



fig1,ax1 = plt.subplots(1)
ax1.set_title(r'Derivative of $\mathrm{f}\, =\, \frac {tr \left ( \delta g_{ij}\otimes \delta g^{ij} \right )}{tr \left ( g_{ij}\otimes \delta g^{ij} \right ) }$',{'color': 'black', 'fontsize': 14})
ax1.plot(r_v[0:20],yline_ana_v[0:20],'b-',linewidth = 1.5)
ax1.plot(r_v[0:20],yline_num_v[0:20],'ro',linewidth = 1.5)
ax1.set_ylabel(r'Derivatives of y',{'color': 'black', 'fontsize': 18})
ax1.set_xlabel(r'Radius r',{'color': 'black', 'fontsize': 18})
plt.setp(ax1.get_xticklabels(), fontsize=16, fontweight="normal")
plt.setp(ax1.get_yticklabels(), fontsize=16, fontweight="normal")
plt.legend((r'Analytical',r'Numerical'),frameon=True,edgecolor='k',shadow=True, loc= 'upper left', fontsize = 16 )
plt.axis([0.95,1.0005,-1.50,-1.30])

plt.show()

