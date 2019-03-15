#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

#sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
input_filename_base_pos = sys.argv[1]
input_filename_base_neg = sys.argv[2]
snapname = sys.argv[3]
xcoord = float(sys.argv[4])
ycoord = float(sys.argv[5])
zcoord = float(sys.argv[6])
lgrid = float(sys.argv[7])
ngrid = int(sys.argv[8])
ion = sys.argv[9]
direction = sys.argv[10]
#'047_z000p000.ioneq'
#tag='ioneq_025'
#center = np.array([6.98,5.21,6.55])
#center= np.array([17.5995,14.08347,15.8329])  #snapshot 31
#center= np.array([16.3672,13.0542,14.6979])  #z=0.271
#center = np.array([15.2534, 10.9404,  9.0412])
#center = np.array([15.2691, 10.934, 9.03164])
#center = np.array([15.2974, 10.9540,  9.0412])
center = np.array([xcoord, ycoord, zcoord])

lgrid = lgrid/1e+03
lgridz = lgrid*2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

redshift = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/Redshift")
aex = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/HubbleParam")
boxsize = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/BoxSize")

boxsize = boxsize/hubble_param*aex
print "boxsize=", boxsize 
center = center/hubble_param*aex
print "center= ", center

coords = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Coordinates",numThreads=1)
mass = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Mass",numThreads=1)
mass *= 1.e+10
print "mass= ", mass
hsmooth = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/SmoothingLength",numThreads=1)
hydrogen = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
chem_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
h1_pos = chem_pos[:,1]
c4_pos = chem_pos[:,10]
o6_pos = chem_pos[:,28]
o7_pos = chem_pos[:,29]
o8_pos = chem_pos[:,30]
mg2_pos = chem_pos[:,45] 

chem_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
h1_neg = chem_neg[:,1]
c4_neg = chem_neg[:,10]
o6_neg = chem_neg[:,28]
o7_neg = chem_neg[:,29]
o8_neg = chem_neg[:,30]
mg2_neg = chem_neg[:,45] 

if(ion=="h1"):
    numerator = mass*hydrogen*h1_pos
    denominator = mass*hydrogen*h1_neg
    ionname = 'HI'
    colmin = 13
    colmax = 21
if(ion=="c4"):
    numerator = mass*hydrogen*c4_pos
    denominator = mass*hydrogen*c4_neg
    ionname = 'CIV'
    colmin = 11
    colmax = 15
if(ion=="o6"):
    numerator = mass*hydrogen*o6_pos
    denominator = mass*hydrogen*o6_neg
    ionname = 'OVI'
    colmin = 11
    colmax = 15
if(ion=="o7"):
    numerator = mass*hydrogen*o7_pos
    denominator = mass*hydrogen*o7_neg
    ionname = 'OVII'
    colmin = 13.5
    colmax = 16.5
if(ion=="o8"):
    numerator = mass*hydrogen*o8_pos
    denominator = mass*hydrogen*o8_neg
    ionname = 'OVIII'
    colmin = 13.5
    colmax = 16.5


if(direction=="z"):
    thetaangle = 0
    phiangle = 0
if(direction=="y"):
    thetaangle = 90
    phiangle = 0
if(direction=="x"):
    thetaangle = 0
    phiangle = 90


result_den = coldens.main(coords, hsmooth, denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_diff.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)

result_num = coldens.main(coords, hsmooth, numerator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='num.png',Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
coldens.make_colourmap('coldens_diff.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -0.5, 0.5, ionname, ngrid, redshift=redshift, extralabel=ionname'$_{\mathrm{NEQ-Equ}}')
print "result_num= ", result_num
print "result_den= ", result_den

