import numpy as np
import sys
import make_maps_v3_master as mmap

# read-in
numslices = int(sys.argv[2])
sliceind = int(sys.argv[1])

# basic box and projection region parameters (splitting into (sub)slices is done later)
simulation = 'eagle'
simnum = 'L0050N0752'
snapnum = 28
var='REFERENCE' 
npix_x=1024
npix_y=1024

L_x = 50.
L_y = 50.
L_z = 50.
centre = [25.,25.,25.]
LsinMpc = True
axis = 'z'

periodic = True
velcut = False
kernel = 'C2'

parttype = '0'

ptypeW = 'basic'
ionW = 'None'
abundsW = 'Pt'
quantityW = 'Mass'
excludeSFRW = 'T4'

ptypeQ = None
ionQ = None
abundsQ = 'auto'
quantityQ = 'Temperature'
excludeSFRQ = 'T4'


saveres = True
log = True
misc = None
ompproj = True

       
         

# box slicing using input parameters (no need to adjust per run)
if axis == 'z':
    L_z = L_z/np.float(numslices) # split into numslices slices along projection axis
    centre[2] = centre[2] - (numslices+1.)*L_z/2. + sliceind*L_z  
if axis == 'x':
    L_x = L_x/np.float(numslices) # split into numslices slices along projection axis
    centre[0] = centre[0] - (numslices+1.)*L_x/2. + sliceind*L_x  
if axis == 'y':
    L_y = L_y/np.float(numslices) # split into numslices slices along projection axis
    centre[1] = centre[1] - (numslices+1.)*L_y/2. + sliceind*L_y      

# document input parameters:(processed input is also printed by make_map)
print('\n')

print('Overview of function input parameters: [cMpc] where applicable \n')
print('simnum: \t %s' %simnum)
print('snapnum: \t %i' % snapnum)
print('centre: \t %s' %str(centre))
print('L_x, L_y, L_z: \t %f, %f, %f \n' %(L_x, L_y, L_z))

print('kernel: \t %s' %kernel)
print('axis: \t %s' %axis)
print('periodic: \t %s' %str(periodic))
print('npix_x,npix_y: \t %i, %i \n' %(npix_x, npix_y))

print('saveres: \t %s' %str(saveres))
print('log: \t %s' %str(True))

print('\n')

# function call
resultW, resultQ = mmap.make_map(simnum, snapnum, centre, L_x, L_y, L_z, npix_x, npix_y, \
         ptypeW,\
         ionW = ionW, abundsW = abundsW, quantityW = quantityW,\
         ionQ = ionQ, abundsQ = abundsQ, quantityQ = quantityQ, ptypeQ = ptypeQ,\
         excludeSFRW = excludeSFRW, excludeSFRQ = excludeSFRQ, parttype = parttype,\
         theta=0.0, phi=0.0, psi=0.0, \
         var = var, axis = axis, log = log, velcut = velcut,\
         periodic = periodic, kernel = kernel, saveres = saveres,\
         simulation = simulation, LsinMpc = LsinMpc,\
         select = None, misc = None, ompproj = ompproj)



#print resultW

#print resultQ




