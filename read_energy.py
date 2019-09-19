"""
Created on Thu Aug 29 17:30:08 2019
@author: Tanvir Chowdhury
"""
import numpy as np
from scipy.interpolate import interp1d
#import matplotlib.pyplot as plt
# ask for input
no_of_C = input('Enter the number of C atom in your Graphene layer: ')

# cell_size = scaling * float(input('Enter the size of the cell: '))

d0 = float(input('Enter the equillibrium distance: '))

distance = np.zeros(7)

#loop to create the distance array
for n in range(-3,4):
    distance[n + 3] = d0 + n * 0.25
    
#read file, skip the first line, read 2nd column
raw_energy = np.loadtxt('energy_data.txt',skiprows=1,usecols=1,max_rows=9,unpack=True)
monomer_energy = np.sum([raw_energy[-1], raw_energy[-2]]) #last two rows, 2nd column
graphene_sub_energy = raw_energy[0:7] #first 6 rows (0 to 6 I mean n-1)

#calculate binding energy
binding_nrg = graphene_sub_energy - monomer_energy
print(binding_nrg)

#divide them by the number of carbon atoms
nrg_per_C = binding_nrg/eval(no_of_C)
print(nrg_per_C)
# interpolation
f1 = interp1d(distance,nrg_per_C)
f2 = interp1d(distance,nrg_per_C,kind='cubic')
distance_new = np.linspace(distance[0],distance[6],50)
print('The minimum value of binding energy per C atom is ',1000* (f2(distance_new).min()),'meV')
ii = np.where(f2(distance_new) == f2(distance_new).min())
print('At the minimum distance of ',distance_new[ii[0]], 'Angstroms')

'''
plt.plot(distance,1000*nrg_per_C,'o',distance_new,1000*f1(distance_new),'-',distance_new,1000*f2(distance_new),'--')
plt.legend(['original data','linear interpolation','cubic interpolation'], loc = 'best')
plt.xlabel('Distance (Angstroms)')
plt.ylabel('Binding Energy per C atom (meV)')
plt.show()
#writing data into a file
info = 'Binding Energy data'
info += '\nDate: 3-September-2019'
info += '\nAuthor: Tanvir Chowdhury'
info += '\n for Graphene-hBN system\n\n'
info += 'Binding Energy (eV)   Binding Energy per C (eV)  Distance (Angstrom)'
np.savetxt('binding_energy.txt',list(zip(binding_nrg,nrg_per_C, distance)), header=info,fmt='%20.8f')'''
