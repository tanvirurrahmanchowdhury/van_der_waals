#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 01:28:31 2020

@author: tanvir
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d,CubicSpline,Akima1DInterpolator

d = np.loadtxt('data.txt',usecols=0,unpack=True)
Evdw = np.loadtxt('data.txt',usecols=1,unpack=True)
#f1 = interp1d(d,Evdw,kind='cubic')
f2 = Akima1DInterpolator(d,Evdw)
d_new = np.arange(3,50,.01)
'''ln_D = np.log(D)
ln_Evdw = np.log(-Evdw)
power_law = (ln_Evdw[1:]-ln_Evdw[:-1])/(ln_D[1:]-ln_D[:-1])
plt.plot(d,Evdw,'o',d_new,f2(d_new),'--',linewidth=2)
plt.xlabel('Adsorption Distance (D/Angstrom)',fontsize=26)
plt.ylabel('vdw Power Law',fontsize=26)
plt.rc('xtick',labelsize=26)
plt.rc('ytick',labelsize=26)
plt.show()'''

#writing data into a file
info = 'Power Law Data'
info += '\nDate: 4-May-2020'
info += '\nAuthor: Tanvir Chowdhury'
info += '\nAkima Interpolation\n'
info += 'D   Evdw'
np.savetxt('power_law.txt',list(zip(d_new,f2(d_new))), header=info,fmt='%20.8f')
