import numpy as np

scaling_factor = float(input('Enter scaling factor: '))
cell_size = float(input('Enter cell size (z value): '))
c_atom = eval(input('Enter the number of C atoms: '))
all_atoms = eval(input('Enter the total number of atoms: '))
#read the CONTCAR file
z_data = np.loadtxt('CONTCAR',skiprows=8,usecols=2,max_rows=all_atoms,unpack=True)
#seperate graphene layer z data
graphene_layer = z_data[0:c_atom]
# find the average position
graphene_layer_center = graphene_layer.mean()

sub_layer = z_data[c_atom:all_atoms]
# sort the array
sub_layer = -np.sort(-sub_layer)
ii = 0
while ii != 1:
    flag = input('Is the first layer of Se atoms? Press y or n: ')
    if flag=='y' or flag=='Y':
        print('Okay! So it is Se')
        upper_sub_layer_center = sub_layer[0:9]
        ii = 1
    elif flag=='n' or flag=='N':
        print('Okay! So it is S')
        upper_sub_layer_center = sub_layer[0:16]
        ii = 1
    else:
        print('Invalid Character. Please enter either y or n.')
del ii
mean_layer = upper_sub_layer_center.mean()
#print(upper_se_layer_center)
'''direct distance between Graphern layer and the
first layer of atoms in the layered material'''
direct_distance = scaling_factor * cell_size * (graphene_layer_center - mean_layer)
print('d_0 = ', direct_distance)
