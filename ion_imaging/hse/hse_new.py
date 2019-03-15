import sys
import tables
import eagle
import h5py
import glob as glob
import os
import numpy as np
import scipy.stats as scist
import galmanip.readGalaxy as rg
import galmanip.rotateGalaxy as rotg
import galmanip.writeGalaxy as wg
import galmanip.binRadial as br
import galmanip.binSpherical as bs
import coldens_ben.coldens as coldens
from astropy import constants as const
import matplotlib.pyplot as plt
import hse_spherical 
import gc

G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24


snap_intag = sys.argv[1]
galidstr = sys.argv[2]
    #snap_outtag = sys.argv[2]
    #sim_name = sys.argv[2]  # e.g. m1e11f30/data
x_gal = float(sys.argv[3])
y_gal = float(sys.argv[4])
z_gal = float(sys.argv[5])
v_x_gal = float(sys.argv[6])
v_y_gal = float(sys.argv[7])
v_z_gal = float(sys.argv[8])

chem = 0
dofigs = 0
cosmoowls = 0
###acceleration = 1

runlabel = 'halo'

unit_mass_in_cgs = 1.989e33 * 1.0e10 
unit_length_in_cgs = 3.0857e24 
proton_mass_cgs = 1.67e-24

unit_Density_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3/proton_mass_cgs 




if (runlabel=='halo'):
    haloname = galidstr.split("halo")[1]
    halo_id = haloname.split("_")[1]
    M200 = float(haloname.split("_")[2])
    ms = float(haloname.split("_")[3])
    sfr = float(haloname.split("_")[4])
    haloinfostr = 'lg M$_{200}=%4.1f$, lg M$_{*}=%3.1f$, SFR$=%4.2f$'%(M200,ms,sfr)
    z_hold = galidstr.split("_z")[1]
    zname = z_hold.split("_")[0]
    zint = zname.split("p")[0]
    zfrac = zname.split("p")[1]
    redshift = float(zint)*1.0 + float(zfrac)/1000.

else:
    haloinfostr = ''


data_dir = "" #% sim_name
sniptag = "SNAP"
sim = "."

#h5file_inbase = "%ssnapshot_%s/snap_%s" % (data_dir, snap_intag, snap_intag)

hubbleparam,aex,redshift,boxsize,mass_table = rg.header(sniptag,sim,snap_intag)

omegaM = 0.307
omegaL = 1-omegaM
omegaratio = (omegaM+omegaL/(1+redshift)**3)
R200 = 1.63e-5*(10**float(M200)*hubbleparam)**0.333/omegaratio**0.333/(1+redshift)/hubbleparam
v200 = np.sqrt(G_Grav*10**M200*M_Solar/(R200*cmpermpc))/1e+05

boxsize_proper = boxsize*aex/hubbleparam
print "boxsize_proper = ", boxsize_proper

gal_coords = ((float(x_gal*aex/hubbleparam),float(y_gal*aex/hubbleparam),float(z_gal*aex/hubbleparam)))
gal_vels = ((float(v_x_gal),float(v_y_gal),float(v_z_gal)))

box_center = ((float(boxsize_proper/2.),float(boxsize_proper/2.),float(boxsize_proper/2.)))

coords_gas,vels_gas = rg.coords(sniptag,sim,snap_intag,0)

gas_indexes = rg.distanceindexes(coords_gas,gal_coords,1.0,boxsize_proper)

coords_gas = coords_gas[gas_indexes]
vels_gas = vels_gas[gas_indexes]

gc.collect()
phi=0
theta=0
psi=0
coords_gas_rot, vels_gas_rot = rotg.rotategalaxy(coords_gas,vels_gas,gal_coords,gal_vels,boxsize_proper,phi,theta,psi)

del(cooords
gc.collect()



mass_gas = rg.mass_indexes(sniptag,sim,snap_intag,0,gas_indexes)*1e+10

gc.collect()



if(cosmoowls):
    nH, T = rg.gasprops_cosmoowls(sniptag,sim,snap_intag)
else:
    nH, T, Z, SFR, hsmooth = rg.gasprops(sniptag,sim,snap_intag)

for i in range(len(nH)):
    if(i%1000000==0): print i, nH[i]

coords_DM,vels_DM = rg.coords(sniptag,sim,snap_intag,1)

coords_stars,vels_stars = rg.coords(sniptag,sim,snap_intag,4)
mass_stars = rg.mass(sniptag,sim,snap_intag,4)*1e+10

