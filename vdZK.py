#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:22:25 2019

@author: Tanvir Chowdhury
"""

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
# importing my own library
from tanvir import beta_model, polarize, diffuse, damping
print('***** Input some parameters for the substrate *****')
scale = float(input('Enter the scale factor: '))
bohr = 0.529177
eps_perp = float(input('Enter the perpendicular component of static dielectric function: '))
eps_para = float(input('Enter the parallel component of static dielectric function: '))
eps0 = (2*eps_perp + eps_para)/3       
a_vector = float(input('Enter substrate unit vector a value, (vector b=a), (Angstroms): '))
a_sub = scale * (a_vector/bohr)
c_vector = float(input('Enter substrate unit vector c value (Angstroms): '))
c_sub = c_vector/bohr
angle_sub = float(input('And what is the angle between vectors a and b in degrees? '))
valence_1 = eval(input('Enter the number of valence electron of Mo or W: '))
valence_2 = eval(input('Enter the number of valence electron of S or Se: '))
no_of_layers = eval(input('Enter the number of layers: '))

n_sub = no_of_layers * (valence_1 + 2*valence_2)/(a_sub ** 2 * c_sub * np.sin(np.radians(angle_sub)))

rs = (3/(4 * np.pi * n_sub)) ** (1/3)
d0_perp = 0.02*rs**2 - 0.27*rs + 2.06
w_p_bar = np.sqrt(4*np.pi*n_sub)
w_sub = w_p_bar/np.sqrt(2)
Ef = 0.5 * (3 * n_sub * np.pi **2)**(2/3)

# penn model to compute wg
def penn_model(x):
       return (eps0-1)**2*(9/4)*x**4/(w_p_bar**4) + (eps0-1)*(3/4)*x**3/(w_p_bar**2*Ef) - 1

# inital guess was that the solution is between 0.1 to 0.9
wg = scipy.optimize.brentq(penn_model,0.1,0.7)

print('***** Thanks! I have modeled the substrate. Not input some parameters for the graphene layer *****')
a_0_1 = float(input('Enter the static dipole polarizability (in a.u.): '))
a_0_2 = a_0_1 ** (5/3)
a_graphene = scale * (float(input('Enter graphene unit vector a value, (vector b=a), (Angstroms): '))/bohr)
c_graphene = float(input('Enter graphene unit vector c value (in Bohrs): '))
angle_graphene = float(input('And what is the angle between vectors a and b in degrees? '))

print('\n Thank you! I have got all the information I need. Now, please wait while I am computing...')

# n_graphene is m in the paper
n_graphene = 4*2/(a_graphene **2 * c_graphene * np.sin(np.radians(angle_graphene)))

wp = np.sqrt(4 * np.pi * n_graphene)
w1 = wp/np.sqrt(3)
w2 = wp * np.sqrt(2/5)
kf = np.sqrt(2*Ef)
delta = np.sqrt(3*(kf**2)/5)
big_delta = wg/(4*Ef)
y = 1/big_delta
P = np.sqrt(1 + y**2)

# game begins (computing the coefficients) C3, C4, C5

C3 = np.zeros(10)
C4d = np.zeros(10)
C5d = np.zeros(10)
C4nl = np.zeros(10)
C5nl = np.zeros(10)
C5q = np.zeros(10)

#image loop

N = 0.001 # number of increment
u = np.arange(N,500,N)
#for n from 1 to 10 images
for n in range(1,11):
    for jj in range(len(u)):
         a_u = np.zeros(2)
         # calling function
         beta, e1 = beta_model(u[jj],w_p_bar,wg,y,big_delta,P)
         # calling another function
         a_u = polarize(u[jj],w1,w2,a_0_1,a_0_2)
         
         # C3 term
         C3[n-1] += N/(4*np.pi) * a_u[0]*beta**(2*n-1)
         
         g0 = beta * (1 + u[jj] ** 2/w_sub ** 2)
         # calling another function
         lamb, du_perp, fw = diffuse(u[jj],g0,w_sub,w1, d0_perp)
         
         #C4, C5 diffusion terms
         C4d[n-1] += N * a_0_1*g0*w1 ** 2*w_sub ** 2*0.5*fw*du_perp*beta**(2*n-2)
         C5d[n-1] += N * a_0_1*g0*w1 ** 2*w_sub ** 4/(2*(w1 ** 2+ w_sub ** 2))*fw*du_perp ** 2*beta **(2*n-2)
         
         # quadruple term
         C5q[n-1] += N/(4*np.pi) * a_u[1] * beta ** (2*n-1)
         zz = delta ** 2 /(w_p_bar ** 2 + u[jj] ** 2)
         
         # non-local terms
         C4nl[n-1] += N * (-3/(8*np.pi))* e1* a_u[0]* (e1-1)/(e1+1)**2* np.sqrt(zz)*beta ** (2*n-2)
         C5nl[n-1] += (3*N/(4*np.pi)) * a_u[0]*e1*((e1-1) ** 2*zz/(e1+1) ** 3)*beta ** (2*n-2)

# for loop (image loop) ends here
C5 = C5d + C5q + C5nl
C4 = C4d + C4nl

# Computing van der Waals Energy
b_bar = 4.5 # cutoff parameter in Bohrs
c = c_sub * bohr
# load data
din, Edft = np.loadtxt('ws2_2L.txt',skiprows=1,unpack=True)
#compute damping factor fd begins
g = 2 * b_bar ** 2 * (C3/C5)
h = 10 * b_bar ** 4 * (C3/C5) ** 2
NN = len(din)
x = np.zeros((10,NN),dtype=float)
term_1 = np.zeros((10,NN),dtype=float)
term_2 = np.zeros((10,NN),dtype=float)
fd = np.zeros((10,NN),dtype=float)

for n in range(1,11):
    term_1[n-1,:] = (din + ((n-1) * c))/bohr
    term_2[n-1,:] = bohr/((din + n * c))
    x[n-1,:] = term_1[n-1,:]/b_bar

# calling the damping function
fd = damping(g,h,x,NN)

# Equ. (10) first term computing
Evdw1 = -np.dot(C3, fd * (1/term_1)**3 ) - np.dot(C4, fd * (1/term_1)**4) - np.dot(C5, fd * (1/term_1)**5)

# Equ. (10) second term computing

Evdw2 = np.dot(C3, fd * term_2 ** 3) + np.dot(C4,fd * term_2 ** 4) + np.dot(C5,fd * term_2 ** 5)

# converting energy into meV
Evdw = 27.2113966 * 1000 * (Evdw1 + Evdw2)

E_total = Edft + Evdw
print('van der Waals Energy = ',Evdw)
print('Total Energy= ',E_total)
# plotting din vs E_total
plt.plot(din,E_total)
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Total Energy (Edft + Evdw)')
plt.show()
#game over