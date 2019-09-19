# -*- coding: utf-8 -*-
import numpy as np
import scipy.optimize
from tanvir import beta_model, polarize, diffuse
scale = 1.0
bohr = 0.529177
eps0 = (2*6.9+3.8)/3
a_sub = scale * (3.29/bohr)
c_sub = 3.2/bohr
n_sub = (1 + 2*6)/(a_sub ** 2 * c_sub * np.cos(np.radians(30)))

rs = (3/4 * np.pi * n_sub) ** (1/3)
d0_perp = 0.02*rs**2 - 0.27*rs + 2.06
w_p_bar = np.sqrt(4*np.pi*n_sub)
w_sub = np.sqrt(w_p_bar/2)
Ef = 0.5 * (3 * n_sub * np.pi **2)**(2/3)

def penn_model(x):
    return eps0-1-((2*w_p_bar **2)*((1+x**2/16*Ef**2) - (x/(4*Ef))))/(3* x**2)

# inital guess was that the solution is between 0.1 to 0.9
wg = scipy.optimize.brentq(penn_model,0.1,0.9)
a_0_1 = 9.945
a_0_2 = a_0_1 ** (5/3)
a_graphene = scale * (2.46/bohr)
c_graphene = 3.4
# n_graphene is m in the paper
n_graphene = 4*2/(a_graphene **2 * c_graphene * np.cos(np.radians(30)))

wp = np.sqrt(4 * np.pi * n_graphene)
w1 = wp/np.sqrt(3)
w2 = wp * np.sqrt(2/5)
kf = np.sqrt(2*Ef)
delta = np.sqrt(3*kf**2/5)
big_delta = wg/(4*Ef)
y = 1/big_delta
P = np.sqrt(1 + y**2)
# game begins Vaan

C3 = np.zeros(10)
C4d = np.zeros(10)
C5d = np.zeros(10)
C4nl = np.zeros(10)
C5nl = np.zeros(10)
C5q = np.zeros(10)

#image loop
N = 0.1
u = np.arange(N,10,N)
#for n from 1 to 10 images
for n in range(1,11):
    for jj in range(len(u)):
         a_u = np.zeros(2)
         beta = beta_model(u[jj],w_p_bar,wg,y,delta,P)
         a_u = polarize(u[jj],w1,w2,a_0_1,a_0_2)
         # C3 term
         C3[n-1] += N/(4*np.pi) * a_u[0]*beta**(2*n-1)
         
         g0 = beta * (1 + u[jj] ** 2/w_sub ** 2)
         lamb, du_perp, fw = diffuse(u[jj],g0,w_sub,w1, d0_perp)
          #C4, C5 diffusion terms
         C4d[n-1] += N * a_0_1*g0*w1 ** 2*w_sub ** 2*0.5*fw*du_perp*beta**(2*n-2)
         C5d[n-1] += N * a_0_1*g0*w1 ** 2*w_sub ** 4/(2*(w1 ** 2+ w_sub ** 2))*fw*du_perp ** 2*beta **(2*n-2)
         
         C5q[n-1] += N/(4*np.pi) * a_u[1] * beta ** (2*n-1)
         zz = delta ** 2 /(w_p_bar ** 2 + u[jj] ** 2)
         
         # non-local terms
         #C4nl[n-1] += N * (-3/(8*np.pi))* e1* a_u[0]*np.sqrt(zz)*beta ** (2*n)
         #C5nl[n-1] += (3*N/4*np.pi) * a_u[0]*e1*((e1-1)^2*zz/(e1+1)^3)*beta ** (2*n-2)
print('C3 = ', C3)
print('C4d = ',C4d)
print('C5d = ',C5d)
print('C5d = ',C5q)