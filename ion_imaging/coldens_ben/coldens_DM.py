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

lgrid = lgrid/1e+03
lgridz = lgrid*4 # Four times longer than Lgrid in zdirection.

if(snip=='0'):
    path = "snapshot_%s/snap_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNAP"
else:
    path = "snipshot_%s/snip_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNIP"

redshift = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/Redshift")
aex = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/HubbleParam")
center = center/hubble_param*aex
print "center= ", center
coords = eagle.readArray(sniptag, sim, input_filename_base, "/PartType1/Coordinates",numThreads=1)
mass_table = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/MassTable")

DM_mass = eagle.readArray(sniptag, sim, input_filename_base, "/PartType1/GroupNumber",numThreads=1)*0.0 + mass_table[1]

SofteningHaloMaxPhys = eagle.readAttribute(sniptag, sim, input_filename_base, "/RuntimePars/SofteningHaloMaxPhys")

hsmooth = DM_mass*0.0 + SofteningHaloMaxPhys

print hsmooth

DM_mass *= 1.e+10 #/hubble_param

WIMP_mass = 106.58

result_dm_z = coldens.main(coords, hsmooth, DM_mass/WIMP_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.DM.z.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='DM(100 GeV)',extralabel='DM',redshift=redshift)

result_dm_y = coldens.main(coords, hsmooth, DM_mass/WIMP_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.DM.y.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='DM(100 GeV)',theta=90,extralabel='DM',redshift=redshift)

result_dm_x = coldens.main(coords, hsmooth, DM_mass/WIMP_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.DM.x.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='DM(100 GeV)',phi=90,extralabel='DM',redshift=redshift)
