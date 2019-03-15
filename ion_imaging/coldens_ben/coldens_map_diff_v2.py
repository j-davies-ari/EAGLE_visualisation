#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
from pylab import get_cmap
from matplotlib.colors import LogNorm
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib import rc

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
runlabel = sys.argv[9]
xvel = float(sys.argv[10])
yvel = float(sys.argv[11])
zvel = float(sys.argv[12])

cmgray = plt.get_cmap('gray_r')

unit_mass_in_cgs = 1.989e33 * 1.0e10 
unit_length_in_cgs = 3.0857e24 
proton_mass_cgs = 1.67e-24
unitDensity_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3/proton_mass_cgs

if (runlabel=='halo'):
    haloname = snapname.split("halo")[1]
    halo_id = haloname.split("_")[1]
    mh = haloname.split("_")[2]
    ms = haloname.split("_")[3]
    sfr = haloname.split("_")[4]
    haloinfostr = 'lg M$_{\mathrm{halo}}=' + mh + '$, lg M$_{*}=' + ms + '$, SFR$=' + sfr + '$'
else:
    haloinfostr = runlabel


center = np.array([xcoord, ycoord, zcoord])

lgrid = lgrid/1e+03
##lgridz = lgrid*0.01 # *2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

redshift = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/Redshift")
aex = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/HubbleParam")
boxsize = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/BoxSize")

boxsize = boxsize/hubble_param*aex
print "boxsize=", boxsize 
center = center/hubble_param*aex
print "center= ", center

coords_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Coordinates",numThreads=1)
velocity = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Velocity",numThreads=1)
mass_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Mass",numThreads=1)*1e+10
print "mass= ", mass_pos

coords_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/Coordinates",numThreads=1)
mass_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/Mass",numThreads=1)*1e+10

gal_coords = ([float(xcoord)*aex/hubble_param,float(ycoord)*aex/hubble_param,float(zcoord)*aex/hubble_param])
gal_vels = ([float(xvel),float(yvel),float(zvel)])

#shiftcoords = np.abs(coords - gal_coords)
#shiftcoords = np.where(shiftcoords > 0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param - shiftcoords, shiftcoords)

shiftcoords = coords_pos - gal_coords
shiftcoords = np.where(shiftcoords > 0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param - shiftcoords, shiftcoords)
shiftcoords = np.where(shiftcoords < -0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param + shiftcoords, shiftcoords)

shiftvels = velocity - gal_vels

vx_mass = shiftvels[:,0]*mass_pos
vy_mass = shiftvels[:,1]*mass_pos
vz_mass = shiftvels[:,2]*mass_pos

vlimit = 600

dist = np.sqrt((shiftcoords[:,0])**2+(shiftcoords[:,1])**2+(shiftcoords[:,2])**2)
vdotroverr = shiftcoords[:,0]*velocity[:,0]+shiftcoords[:,1]*velocity[:,1]+shiftcoords[:,2]*velocity[:,2]
vdotroverr /= dist
vdotroverr_mass = vdotroverr*mass_pos
print "shiftcoords[:,0]= ", shiftcoords[:,0]
print "xcoord= ", xcoord
print "dist= ", dist, min(dist), max(dist)
print "shiftvels[:,0]= ", shiftvels[:,0]
print "vdotroverr= ", vdotroverr, min(vdotroverr), max(vdotroverr)

temperature = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Temperature",physicalUnits=1,noH=0,numThreads=1)
temp_mass = temperature*mass_pos #/nh**(2/3.)*internalenergy/(entropy*3/2.)
print "temp= ", min(temperature), max(temperature), temperature

print "temp_mass= ", min(temp_mass), max(temp_mass), temp_mass
#entropy = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Entropy",numThreads=1)
#internalenergy = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/InternalEnergy",numThreads=1)
nH = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Density",numThreads=1)
nH =  nH*unitDensity_in_cgs
nH_mass = nH*mass_pos

hsmooth_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/SmoothingLength",numThreads=1)
hsmooth_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/SmoothingLength",numThreads=1)
hydrogen_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
hydrogen_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
chem_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
h1_pos = chem_pos[:,1]
he2_pos = chem_pos[:,5]
c4_pos = chem_pos[:,10]
n5_pos = chem_pos[:,19]
o1_pos = chem_pos[:,23]
o2_pos = chem_pos[:,24]
o3_pos = chem_pos[:,25]
o4_pos = chem_pos[:,26]
o5_pos = chem_pos[:,27]
o6_pos = chem_pos[:,28]
o7_pos = chem_pos[:,29]
o8_pos = chem_pos[:,30]
o9_pos = chem_pos[:,31]
ne8_pos = chem_pos[:,40]
mg2_pos = chem_pos[:,45] 
mg10_pos = chem_pos[:,53] 
si3_pos = chem_pos[:,59] 
si4_pos = chem_pos[:,60] 
si12_pos = chem_pos[:,68] 
fe17_pos = chem_pos[:,126] 
fe20_pos = chem_pos[:,129] 
fe23_pos = chem_pos[:,132] 
fe25_pos = chem_pos[:,134] 

chem_pos = np.delete(chem_pos, np.s_[::1])
#del chem_pos
#chem_pos = chem_pos[:1]

chem_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
h1_neg = chem_neg[:,1]
he2_neg = chem_neg[:,5]
c4_neg = chem_neg[:,10]
n5_neg = chem_neg[:,19]
o1_neg = chem_neg[:,23]
o2_neg = chem_neg[:,24]
o3_neg = chem_neg[:,25]
o4_neg = chem_neg[:,26]
o5_neg = chem_neg[:,27]
o6_neg = chem_neg[:,28]
o7_neg = chem_neg[:,29]
o8_neg = chem_neg[:,30]
o9_neg = chem_neg[:,31]
ne8_neg = chem_neg[:,40]
mg2_neg = chem_neg[:,45] 
mg10_neg = chem_neg[:,53]
si3_neg = chem_neg[:,59] 
si4_neg = chem_neg[:,60] 
si12_neg = chem_neg[:,68] 
fe17_neg = chem_neg[:,126] 
fe20_neg = chem_neg[:,129] 
fe23_neg = chem_neg[:,132] 
fe25_neg = chem_neg[:,134] 



chem_neg = np.delete(chem_neg, np.s_[::1])
#del chem_neg
#chem_neg = chem_neg[:1]


ndim = 3
nspecies = 16

particle_IDs_gas = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ParticleIDs",numThreads=1)
metallicity_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Metallicity",numThreads=1)
SFR_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/StarFormationRate",numThreads=1)
oxygen_pos = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
SFR_neg = eagle.readArray("SNAP", sim, input_filename_base_neg, "/PartType0/StarFormationRate",numThreads=1)

#OnEquationOfState = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/OnEquationOfState", numThreads=1)


#fpart = file('particles_outflowspast.%s.l%3.1f.dat'%(snapname,lgrid), 'w')
#for i in xrange(len(particle_IDs_gas)):
#    if((dist[i]<0.3) & (vdotroverr[i] > 200.0) & ((particle_IDs_gas[i]+1)%2==0)):
#        fpart.write("%17d % 5.3f %5.3f %6.4f %5.3e %7.4f % 7.1f % 7.5f % 5.2f"%(particle_IDs_gas[i], np.log10(nh[i]), np.log10(temperature[i]), metallicity[i], np.log10(SFR[i]+1e-20),dist[i], vdotroverr[i], OnEquationOfState[i], np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]+o5_pos[i]+o6_pos[i]+o7_pos[i]+o8_pos[i]+o9_pos[i])))
#        fpart.write("  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f\n"%(np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]+1.1e-10), np.log10(o1_neg[i]+o2_neg[i]+o3_neg[i]+o4_neg[i]+1.1e-10), np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]) - np.log10(o1_neg[i]+o2_neg[i]+o3_neg[i]+o4_neg[i]), np.log10(o5_pos[i]+1.1e-10), np.log10(o5_neg[i]+1.1e-10), np.log10(o5_pos[i]) - np.log10(o5_neg[i]), np.log10(o6_pos[i]+1.1e-10), np.log10(o6_neg[i]+1.1e-10), np.log10(o6_pos[i]) - np.log10(o6_neg[i]), np.log10(o7_pos[i]+1.1e-10), np.log10(o7_neg[i]+1.1e-10), np.log10(o7_pos[i]) - np.log10(o7_neg[i]), np.log10(o8_pos[i]+1.1e-10), np.log10(o8_neg[i]+1.1e-10), np.log10(o8_pos[i]) - np.log10(o8_neg[i]), np.log10(o9_pos[i]+1.1e-10), np.log10(o9_neg[i]+1.1e-10), np.log10(o9_pos[i]) - np.log10(o9_neg[i])))

#fpart.close()


#fpart = file('particles_odiff.%s.l%3.1f.dat'%(snapname,lgrid), 'w')

#for i in xrange(len(o6_neg)):
#    if( ((particle_IDs[i]+1)%20==0) & (np.log10(o6_pos[i]+o7_pos[i]+o8_pos[i])>-3.61) & ((np.fabs(np.log10(o6_pos[i]) - np.log10(o6_neg[i])) > 0.3) or (np.fabs(np.log10(o7_pos[i]) - np.log10(o7_neg[i])) > 0.3) or  (np.fabs(np.log10(o8_pos[i]) - np.log10(o8_neg[i])) > 0.3)) ):
#        fpart.write("%17d % 5.3f %5.3f %6.4f %5.3e %7.4f % 7.1f % 7.5f % 5.2f"%(particle_IDs[i], np.log10(nh[i]), np.log10(temperature[i]), metallicity[i], np.log10(SFR[i]+1e-20),dist[i], vdotroverr[i], OnEquationOfState[i], np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]+o5_pos[i]+o6_pos[i]+o7_pos[i]+o8_pos[i]+o9_pos[i])))
#        fpart.write("  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f  % 5.2f % 5.2f % 5.2f\n"%(np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]+1.1e-10), np.log10(o1_neg[i]+o2_neg[i]+o3_neg[i]+o4_neg[i]+1.1e-10), np.log10(o1_pos[i]+o2_pos[i]+o3_pos[i]+o4_pos[i]) - np.log10(o1_neg[i]+o2_neg[i]+o3_neg[i]+o4_neg[i]), np.log10(o5_pos[i]+1.1e-10), np.log10(o5_neg[i]+1.1e-10), np.log10(o5_pos[i]) - np.log10(o5_neg[i]), np.log10(o6_pos[i]+1.1e-10), np.log10(o6_neg[i]+1.1e-10), np.log10(o6_pos[i]) - np.log10(o6_neg[i]), np.log10(o7_pos[i]+1.1e-10), np.log10(o7_neg[i]+1.1e-10), np.log10(o7_pos[i]) - np.log10(o7_neg[i]), np.log10(o8_pos[i]+1.1e-10), np.log10(o8_neg[i]+1.1e-10), np.log10(o8_pos[i]) - np.log10(o8_neg[i]), np.log10(o9_pos[i]+1.1e-10), np.log10(o9_neg[i]+1.1e-10), np.log10(o9_pos[i]) - np.log10(o9_neg[i])))
#fpart.close()


for i in xrange(ndim):
    if (i==0):
        direction="x"
        thetaangle = 0
        phiangle = 90
    if (i==1):
        direction="y"
        thetaangle = 90
        phiangle = 0
    if (i==2):
        direction="z"
        thetaangle = 0
        phiangle = 0

    lgridz = lgrid*0.01 # *2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

    result_Nh = coldens.main(coords_pos, hsmooth_pos, mass_pos*hydrogen_pos, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.hydrogen.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=18.0+np.log10(lgridz/2./lgrid), Vmax=22.0+np.log10(lgridz/2./lgrid),ion='H',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='H')
    print "result_Nh= ", result_Nh
    result_mass = coldens.main(coords_pos, hsmooth_pos, mass_pos, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.mass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=18.0+np.log10(lgridz/2./lgrid), Vmax=22.0+np.log10(lgridz/2./lgrid),ion='mass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='mass')
    result_tempmass = coldens.main(coords_pos, hsmooth_pos, temp_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.Tmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=22.0, Vmax=29.0,ion='Tmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='T')
    coldens.make_colourmap('coldens_num.%s.T.%s.l%3.1f.png'%(snapname,direction,lgrid), result_tempmass-result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, 4.0,7.0, "T", ngrid, redshift=redshift, extralabel='T',haloinfostr=haloinfostr)
    print "result_tempmass= ", result_tempmass
    result_nHmass = coldens.main(coords_pos, hsmooth_pos, nH_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.nHmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=14.0, Vmax=22.0,ion='nHmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    coldens.make_colourmap('coldens_num.%s.nH.%s.l%3.1f.png'%(snapname,direction,lgrid), result_nHmass-result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -7.0,-1.0, "nH", ngrid, redshift=redshift, extralabel='n$_{\mathrm{H}}$',haloinfostr=haloinfostr)
    print "result_nHmass= ", result_nHmass

    result_vdotroverrmass = coldens.main(coords_pos, hsmooth_pos, vdotroverr_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vdotroverrmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vdotroverrmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    print "result_vdotroverrmass= ", result_vdotroverrmass
    print "10**result_mass= ", 10**result_mass
    coldens.make_colourmap('coldens_num.%s.vdotroverr.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vdotroverrmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vdotroverr", ngrid, redshift=redshift, extralabel='v$_\mathrm{rad}$',haloinfostr=haloinfostr)
    print "result_vdotroverrmass/10**result_mass= ", result_vdotroverrmass/10**result_mass


    #result_vxmass = coldens.main(coords_pos, hsmooth_pos, vx_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vxmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vxmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    #print "result_vxmass= ", result_vxmass
    #print "10**result_mass= ", 10**result_mass
    #coldens.make_colourmap('coldens_num.%s.vx.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vxmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vx", ngrid, redshift=redshift, extralabel='v$_{x}$',haloinfostr=haloinfostr)
    #print "result_vxmass/10**result_mass= ", result_vxmass/10**result_mass

    #result_vymass = coldens.main(coords_pos, hsmooth_pos, vy_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vymass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vymass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    #coldens.make_colourmap('coldens_num.%s.vy.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vymass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vy", ngrid, redshift=redshift, extralabel='v$_{y}$',haloinfostr=haloinfostr)

    #result_vzmass = coldens.main(coords_pos, hsmooth_pos, vz_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vzmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vzmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    #coldens.make_colourmap('coldens_num.%s.vz.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vzmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vz", ngrid, redshift=redshift, extralabel='v$_{z}$',haloinfostr=haloinfostr)

    lgridz = lgrid*2 #*2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

    for j in xrange(nspecies):
      if (j>1):
        indexes_pos = np.where(SFR_pos==0)
        indexes_neg = np.where(SFR_neg==0)
        
        if(j==0):
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*h1_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*h1_neg[indexes_neg]
            ionname = 'HI'
            ion = 'h1'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 18.0+np.log10(lgridz/2./lgrid)
#        if(j==1):
#            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*he2_pos[indexes_pos]
#            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*he2_neg[indexes_neg]
#            ionname = 'HeII'
#            ion = 'he2'
#            colmin = 13.0+np.log10(lgridz/2./lgrid)
#            colmax = 18.0+np.log10(lgridz/2./lgrid)
        if(j==1):
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*(o1_pos[indexes_pos]+o2_pos[indexes_pos]+o3_pos[indexes_pos]+o4_pos[indexes_pos]+o5_pos[indexes_pos]+o6_pos[indexes_pos]+o7_pos[indexes_pos]+o8_pos[indexes_pos]+o9_pos[indexes_pos])
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*(o1_neg[indexes_neg]+o2_neg[indexes_neg]+o3_neg[indexes_neg]+o4_neg[indexes_neg]+o5_neg[indexes_neg]+o6_neg[indexes_neg]+o7_neg[indexes_neg]+o8_neg[indexes_neg]+o9_neg[indexes_neg])
            ionname = 'O'
            ion = 'o'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 18.0+np.log10(lgridz/2./lgrid)
        if(j==2):
            indexes_pos = np.where((SFR_pos==0) & (c4_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (c4_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*c4_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*c4_neg[indexes_neg]
            ionname = 'CIV'
            ion = 'c4'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==3):
            indexes_pos = np.where((SFR_pos==0) & (n5_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (n5_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*n5_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*n5_neg[indexes_neg]
            ionname = 'NV'
            ion = 'n5'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==4):
            indexes_pos = np.where((SFR_pos==0) & (o6_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (o6_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*o6_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*o6_neg[indexes_neg]
            ionname = 'OVI'
            ion = 'o6'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==5):
            indexes_pos = np.where((SFR_pos==0) & (o7_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (o7_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*o7_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*o7_neg[indexes_neg]
            ionname = 'OVII'
            ion = 'o7'
            colmin = 12.5+np.log10(lgridz/2./lgrid)
            colmax = 16.5+np.log10(lgridz/2./lgrid)
        if(j==6):
            indexes_pos = np.where((SFR_pos==0) & (o8_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (o8_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*o8_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*o8_neg[indexes_neg]
            ionname = 'OVIII'
            ion = 'o8'
            colmin = 12.5+np.log10(lgridz/2./lgrid)
            colmax = 16.5+np.log10(lgridz/2./lgrid)
        if(j==7):
            indexes_pos = np.where((SFR_pos==0) & (o9_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (o9_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*o9_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*o9_neg[indexes_neg]
            ionname = 'OIX'
            ion = 'o9'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 17.0+np.log10(lgridz/2./lgrid)
        if(j==8):
            indexes_pos = np.where((SFR_pos==0) & (ne8_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (ne8_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*ne8_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*ne8_neg[indexes_neg]
            ionname = 'NeVIII'
            ion = 'ne8'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==9):
            indexes_pos = np.where((SFR_pos==0) & (mg2_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (mg2_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*mg2_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*mg2_neg[indexes_neg]
            ionname = 'MgII'
            ion = 'mg2'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==10):
            indexes_pos = np.where((SFR_pos==0) & (mg10_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (mg10_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*mg10_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*mg10_neg[indexes_neg]
            ionname = 'MgX'
            ion = 'mg10'            
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==11):
            indexes_pos = np.where((SFR_pos==0) & (si3_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (si3_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*si3_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*si3_neg[indexes_neg]
            ionname = 'SiIII'
            ion = 'si3'
            colmin = 12.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==12):
            indexes_pos = np.where((SFR_pos==0) & (si4_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (si4_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*si4_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*si4_neg[indexes_neg]
            ionname = 'SiIV'
            ion = 'si4'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==13):
            indexes_pos = np.where((SFR_pos==0) & (si12_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (si12_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*si12_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*si12_neg[indexes_neg]
            ionname = 'SiXII'
            ion = 'si12'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==14):
            indexes_pos = np.where((SFR_pos==0) & (fe17_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (fe17_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*fe17_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*fe17_neg[indexes_neg]
            ionname = 'FeXVII'
            ion = 'fe17'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==15):
            indexes_pos = np.where((SFR_pos==0) & (fe20_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (fe20_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*fe20_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*fe20_neg[indexes_neg]
            ionname = 'FeXX'
            ion = 'fe20'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==16):
            indexes_pos = np.where((SFR_pos==0) & (fe23_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (fe23_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*fe23_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*fe23_neg[indexes_neg]
            ionname = 'FeXXIII'
            ion = 'fe23'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==17):
            indexes_pos = np.where((SFR_pos==0) & (fe25_pos>1e-30))
            indexes_neg = np.where((SFR_neg==0) & (fe25_neg>1e-30))
            numerator = mass_pos[indexes_pos]*hydrogen_pos[indexes_pos]*fe25_pos[indexes_pos]
            denominator = mass_neg[indexes_neg]*hydrogen_neg[indexes_neg]*fe25_neg[indexes_neg]
            ionname = 'FeXXV'
            ion = 'fe25'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)

        numerator_T = numerator*temperature[indexes_pos]
        numerator_nH = numerator*nH[indexes_pos]
        numerator_vdotroverr = numerator*vdotroverr[indexes_pos]


        result_den = coldens.main(coords_neg[indexes_neg], hsmooth_neg[indexes_neg], denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_den.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
        result_num = coldens.main(coords_pos[indexes_pos], hsmooth_pos[indexes_pos], numerator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
        coldens.make_colourmap('coldens_diff.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -0.5, 0.5, ionname, ngrid, redshift=redshift, extralabel=ionname + '$_\mathrm{NEQ-Equ}$',haloinfostr=haloinfostr)

        result_num_vdotroverr = coldens.main(coords_pos[indexes_pos], hsmooth_pos[indexes_pos], numerator_vdotroverr, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.%s_vdotroverr.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
        print "result_num_vdotroverr= ", result_num_vdotroverr
        print "10**result_num= ", 10**result_num
        coldens.make_colourmap('coldens_num.%s.%s_vdotroverr.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num_vdotroverr/10**result_num, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit, vlimit, "vdotroverr", ngrid, redshift=redshift, extralabel='v$_\mathrm{rad,' + ionname + '}$',haloinfostr=haloinfostr)
        print "result_num_vdotroverr/10**result_num= ", result_num_vdotroverr/10**result_num

    #result_vdotroverrmass = coldens.main(coords_pos, hsmooth_pos, vdotroverr_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vdotroverrmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vdotroverrmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    #print "result_vdotroverrmass= ", result_vdotroverrmass
    #print "10**result_mass= ", 10**result_mass
    #coldens.make_colourmap('coldens_num.%s.vdotroverr.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vdotroverrmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vdotroverr", ngrid, redshift=redshift, extralabel='v$_\mathrm{rad}$',haloinfostr=haloinfostr)
    #print "result_vdotroverrmass/10**result_mass= ", result_vdotroverrmass/10**result_mass



        result_num_T = coldens.main(coords_pos[indexes_pos], hsmooth_pos[indexes_pos], numerator_T, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.%s_T.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
        coldens.make_colourmap('coldens_num.%s.%s_T.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num_T-result_num, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, 4.0, 7.0, "T", ngrid, redshift=redshift, extralabel='T$_\mathrm{' + ionname + '}$',haloinfostr=haloinfostr)

        result_num_nH = coldens.main(coords_pos[indexes_pos], hsmooth_pos[indexes_pos], numerator_nH, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.%s_nH.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
        coldens.make_colourmap('coldens_num.%s.%s_nH.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num_nH-result_num, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -7.0, -1.0, "nH", ngrid, redshift=redshift, extralabel='n$_\mathrm{H,' + ionname + '}$',haloinfostr=haloinfostr)


        f = file('colion_map_diff.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid), 'w')
        for ii in xrange(ngrid):
            for jj in xrange(ngrid):
                f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 6.1f % 6.1f\n'%(ii,jj,result_num[ii,jj],result_den[ii,jj],result_Nh[ii,jj],result_nHmass[ii,jj]-result_mass[ii,jj],result_tempmass[ii,jj]-result_mass[ii,jj],result_num_nH[ii,jj]-result_num[ii,jj],result_num_T[ii,jj]-result_num[ii,jj],result_vdotroverrmass[ii,jj]/10**result_mass[ii,jj], result_num_vdotroverr[ii,jj]/10**result_num[ii,jj]))
        f.close()

        fig = plt.figure(figsize = (5, 5))
        fig.clf()

        Tlow = 4.0
        Thi = 7.0
        nhlow = -6.0
        nhhi = -1.0
    
        ax1 = fig.add_subplot(111)
        
        hex1 = ax1.hexbin(np.reshape(result_num_nH-result_num,ngrid*ngrid),np.reshape(result_num_T-result_num,ngrid*ngrid), gridsize=100, C=np.reshape(result_num-result_den,ngrid*ngrid),vmin=-0.5,vmax=0.5,cmap=get_cmap("RdBu"))
        ax1.axis([nhlow,nhhi,Tlow,Thi])
        ax1.set_xlabel(r'log $n_{\mathrm{H}}$ [cm$^{-3}$]',fontsize=14)
        ax1.set_ylabel(r'log $T$ [K]',fontsize=14)

        ax2 = ax1.twinx()
        ax2.yaxis.set_major_formatter( NullFormatter() )
        nm,binsm,patches = ax2.hist(np.reshape(result_num_nH-result_num,ngrid*ngrid),40, normed=0,range=[nhlow,nhhi], weights=np.reshape(10**result_num,ngrid*ngrid), histtype='stepfilled',alpha=0.25,color='black')
        ax2.set_xlim(nhlow,nhhi)
        ax2.set_ylim(0,np.max(nm)*5.5/2.)
        ax2.set_xticks([-6,-5,-4,-3,-2,-1])

        ax3 = ax1.twiny()
        ax3.xaxis.set_major_formatter( NullFormatter() )
        nm,binsm,patches = ax3.hist(np.reshape(result_num_T-result_num,ngrid*ngrid),40, normed=0,range=[Tlow,Thi], weights=np.reshape(10**result_num,ngrid*ngrid), histtype='stepfilled',alpha=0.25,orientation='horizontal',color='black')
        ax3.axis([0,2,Tlow,Thi])
        ax3.set_xlim(0, np.max(nm)*5.5/2.)
        ax3.set_yticks([4,5,6,7])


        fig.subplots_adjust(left=0.135, bottom=0.135)
        plt.savefig('rhoT.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid))

        fig = plt.figure(figsize = (5, 5))
        fig.clf()
        ax_N_diff = fig.add_subplot(111)


        indexes_circ = np.where(np.sqrt((ii-ngrid/2.)**2+(jj-ngrid/2.)**2)<ngrid/2.)
        #result_num_circ = result_num[indexes_circ] 

        ax_N_diff.hist2d(np.reshape(result_num,ngrid*ngrid), np.reshape(result_num-result_den,ngrid*ngrid), bins=100,cmap=cmgray, range=[[colmin,colmax],[-1.0,1.0]], norm=mpl.colors.LogNorm()) 
        #hexdiff = ax_N_diff.hexbin(np.reshape(result_num,ngrid*ngrid), np.reshape(result_num-result_den,ngrid*ngrid), gridsize=hexspacing, norm=mpl.colors.LogNorm(), cmap=cmgray, reduce_C_function=np.nansum)


        ax2 = ax_N_diff.twinx()
        ax2.yaxis.set_major_formatter( NullFormatter() )
        nm,binsm,patches = ax2.hist(np.reshape(result_num,ngrid*ngrid),40, normed=0,range=[colmin,colmax], histtype='stepfilled',alpha=0.25,color='red')
        ax2.set_xlim(colmin,colmax)
        ax2.set_ylim(0,np.max(nm)*5.5/2.)

        ax3 = ax_N_diff.twiny()
        ax3.xaxis.set_major_formatter( NullFormatter() )
        nm,binsm,patches = ax3.hist(np.reshape(result_num-result_den,ngrid*ngrid),40, normed=0,range=[-1.0,1.0], histtype='stepfilled',alpha=0.25,orientation='horizontal',color='red')
        #ax3.axis([0,2,Tlow,Thi])
        ax3.set_xlim(0, np.max(nm)*5.5/2.)
        ax3.set_ylim(-1.0,1.0)
        
        ax_N_diff.axis([colmin,colmax,-1.0,1.0])
        ax_N_diff.set_xlabel('log N$_{\mathrm{' + ionname + '}}$',fontsize=14)
        ax_N_diff.set_ylabel('$\Delta$log N$_{\mathrm{NEQ-Equ}}$',fontsize=14)

        fig.subplots_adjust(left=0.16, bottom=0.135)
        plt.savefig('N_diff.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid))

        fig.clf()
        
        fig = plt.figure(figsize = (5, 5))
        ax_N_hist = fig.add_subplot(111)

        #N_bins_num, N_ave_num = np.histogram(np.reshape(result_num,ngrid*ngrid),20, range=[colmin,colmax])
        ax_N_hist.hist(np.reshape(result_den,ngrid*ngrid),20, range=[colmin,colmax], log=1, color='blue', histtype='step', label=r'Equ %s'%ionname)
        ax_N_hist.hist(np.reshape(result_num,ngrid*ngrid),20, range=[colmin,colmax], log=1, color='red', histtype='step', label=r'NEQ %s'%ionname)
        #ax_N_hist.axis([colmin,colmax,-1.0,1.0])
        ax_N_hist.set_xlabel('log N$_{\mathrm{' + ionname + '}}$',fontsize=14)
        ax_N_hist.set_ylabel('f',fontsize=14)
        ax_N_hist.legend(loc='lower left', fontsize=10)

        fig.subplots_adjust(left=0.135, bottom=0.135)
        plt.savefig('N_hist.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid))
