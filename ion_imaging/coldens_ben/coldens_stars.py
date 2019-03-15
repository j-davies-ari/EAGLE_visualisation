#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

sim= '.'
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
snip = sys.argv[3]
xcoord = float(sys.argv[4])
ycoord = float(sys.argv[5])
zcoord = float(sys.argv[6])
lgrid = float(sys.argv[7])
ngrid = int(sys.argv[8])

center = np.array([xcoord, ycoord, zcoord])

pixelpc = lgrid*1.e+03/ngrid

lgrid = lgrid/1e+03
lgridz = lgrid*4 # Four times longer than Lgrid in zdirection.

if(snip=='0'):
    path = "snapshot_%s/snap_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNAP"
else:
    path = "snipshot_%s/snip_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNIP"

aex = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/HubbleParam")
center = center/hubble_param*aex
print "center= ", center

coords = eagle.readArray(sniptag, sim, input_filename_base, "/PartType4/Coordinates",numThreads=1)

stellar_mass = eagle.readArray(sniptag, sim, input_filename_base, "/PartType4/Mass",numThreads=1)

hsmooth = eagle.readArray(sniptag, sim, input_filename_base, "/PartType4/SmoothingLength",numThreads=1)

aexform = eagle.readArray(sniptag, sim, input_filename_base, "/PartType4/StellarFormationTime",numThreads=1)

index_young = np.where(aexform/aex>0.95) # Last 5% of time.  

cm_per_pc = 3.086e+18

nhydr_per_msol = 1.989e+33/1.673e-24

conversion = cm_per_pc**2/nhydr_per_msol

#mass_to_density = pixelpc**-2/npartcm2tomsolpc2
hsmooth = hsmooth/10.

print hsmooth

stellar_mass *= 1.e+10 #/hubble_param

result_stars_z = coldens.main(coords, hsmooth, stellar_mass*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars.z.l%5.3f.png'%(snapname,lgrid),Vmin=-1, Vmax=4,ion='Msol/pc^-2',npix=ngrid)

result_stars_y = coldens.main(coords, hsmooth, stellar_mass*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars.y.l%5.3f.png'%(snapname,lgrid),Vmin=-1, Vmax=4,ion='Msol/pc^-2',theta=90,npix=ngrid)

result_stars_x = coldens.main(coords, hsmooth, stellar_mass*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars.x.l%5.3f.png'%(snapname,lgrid),Vmin=-1, Vmax=4,ion='Msol/pc^-2',phi=90,npix=ngrid)

result_stars_young_z = coldens.main(coords[index_young], hsmooth[index_young], stellar_mass[index_young]*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars_young.z.l%5.3f.png'%(snapname,lgrid),Vmin=-2, Vmax=3,ion='Msol/pc^-2',npix=ngrid)

result_stars_young_y = coldens.main(coords[index_young], hsmooth[index_young], stellar_mass[index_young]*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars_young.y.l%5.3f.png'%(snapname,lgrid),Vmin=-2, Vmax=3,ion='Msol/pc^-2',theta=90,npix=ngrid)

result_stars_young_x = coldens.main(coords[index_young], hsmooth[index_young], stellar_mass[index_young]*conversion, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.stars_young.x.l%5.3f.png'%(snapname,lgrid),Vmin=-2, Vmax=3,ion='Msol/pc^-2',phi=90,npix=ngrid)
