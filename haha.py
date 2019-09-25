#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 15:22:25 2019

@author: Tanvir Chowdhury
"""

import numpy as np
import scipy.optimize
import matplotlib.pyplot as plt
from tanvir import beta_model, polarize, diffuse, damping
print('***** Input some parameters for the substrate *****')
scale = float(input('Enter the scale factor: '))
bohr = 0.529177
eps_perp = float(input('Enter the perpendicular component of static dielectric function: '))
eps_para = float(input('Enter the parallel component of static dielectric function: '))
#   eps0 = (eps_perp + 2*eps_para)/3       
#  a_vector = float(input('Enter substrate unit vector a value, (vector b=a), (Angstroms): '))
#  a_sub = scale * (a_vector/bohr)
#  c_vector = float(input('Enter substrate unit vector c value (Angstroms): '))
#  c_sub = c_vector/bohr
#  angle_sub = float(input('And what is the angle between vectors a and b in degrees? '))
valence_1 = eval(input('Enter the number of valence electron of Mo or W: '))
valence_2 = eval(input('Enter the number of valence electron of S or Se: '))

n_sub = (valence_1 + 2*valence_2)/(a_sub ** 2 * c_sub * np.cos(np.radians(angle_sub)))

#  rs = (3/(4 * np.pi * n_sub)) ** (1/3)
d0_perp = 0.02*rs**2 - 0.27*rs + 2.06
w_p_bar = np.sqrt(4*np.pi*n_sub)
w_sub = np.sqrt(w_p_bar/2)
Ef = 0.5 * (3 * n_sub * np.pi **2)**(2/3)

def penn_model(x):
#       return eps0-1-((2*w_p_bar **2)*((1+x**2/16*Ef**2) - (x/(4*Ef))))/(3* x**2)
       return (eps0-1)**2*(9/4)*x**4/(w_p_bar**4) + (eps0-1)*(3/4)*x**3/(w_p_bar**2*Ef) - 1
# inital guess was that the solution is between 0.1 to 0.9
wg = scipy.optimize.brentq(penn_model,0.1,0.9)

print('***** Thanks! I have modeled the substrate. Not input some parameters for the graphene layer *****')
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
N = 0.001
u = np.arange(N,1000,N)
#for n from 1 to 10 images
for n in range(1,11):
    for jj in range(len(u)):
         a_u = np.zeros(2)
         beta, e1 = beta_model(u[jj],w_p_bar,wg,y,delta,P)
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
         C4nl[n-1] += N * (-3/(8*np.pi))* e1* a_u[0]*np.sqrt(zz)*beta ** (2*n)
         C5nl[n-1] += (3*N/4*np.pi) * a_u[0]*e1*((e1-1) ** 2*zz/(e1+1) ** 3)*beta ** (2*n-2)

'''print('C3 = ', C3)
print('C4d = ',C4d)
print('C5d = ',C5d)
print('C5q = ',C5q)
print('C4nl = ',C4nl)
print('C5nl = ',C5nl)'''
C5 = C5d + C5q + C5nl
C4 = C4d + C4nl
'''print('C4 = ',C4)
print('C5 = ',C5)'''

#E_vdw
b_bar = 4.5
c = c_sub * bohr
# load distance data
din, Edft = np.loadtxt('test_data.txt',skiprows=1,max_rows=21,unpack=True)
#compute damping factor fd begins
g = 2 * b_bar ** 2 * (C3/C5)
h = 10 * b_bar ** 4 * (C3/C5) ** 2
NN = len(din)
x = np.zeros((10,NN),dtype=float)
term_1 = np.zeros((10,NN),dtype=float)
term_2 = np.zeros((10,NN),dtype=float)
fd = np.zeros((10,NN),dtype=float)

#will try to avoid this loop
for n in range(1,11):
    term_1[n-1,:] = din + ((n-1) * c)/bohr
    term_2[n-1,:] = 1/(bohr * (din + n * c))
    x[n-1,:] = term_1[n-1,:]/(bohr * b_bar)

fd = damping(g,h,x,NN)
#fd = damping(g,h,x)
#print('g= ',g.shape)
#print('h= ',h.shape)
#print('fd= ',fd.shape)

# Equ. (10) for Vaan
Vaan = np.dot(C3,-1/term_1 ** 3) + np.dot(C4,-1/term_1 ** 4) + np.dot(C5,-1/term_1 ** 5)
fdT = fd.transpose()
Evdw1 = np.dot(Vaan,fdT)
Evdw1 = np.sum(Evdw1)
#print(Evdw1.shape)
# Equ. (10) for Vbbn
Vbbn = np.dot(C3,term_2 ** 3) + np.dot(C4,term_2 ** 4) + np.dot(C5,term_2 ** 5)
Evdw2 = np.dot(Vbbn, fdT)
Evdw2 = np.sum(Evdw2)

Evdw = 27.2113966 * 1000 * (Evdw1 + Evdw2)
#print(Evdw.shape)

E_total = Edft + Evdw
print('Total Energy= ',E_total)
# plot din vs E_total
plt.plot(din,E_total)
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Total Energy (Edft + Evdw)')
plt.show()
#game over
