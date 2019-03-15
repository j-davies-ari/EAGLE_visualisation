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

#sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
input_filename_base_pos = sys.argv[1]
#input_filename_base_neg = sys.argv[2]
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
lgridz = lgrid*0.01 # *2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

redshift = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/Redshift")
aex = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/HubbleParam")
boxsize = eagle.readAttribute("SNAP", sim, input_filename_base_pos, "/Header/BoxSize")

boxsize = boxsize/hubble_param*aex
print "boxsize=", boxsize 
center = center/hubble_param*aex
print "center= ", center

coords = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Coordinates",numThreads=1)
velocity = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Velocity",numThreads=1)
mass = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Mass",numThreads=1)*1e+10
print "mass= ", mass
gal_coords = ([float(xcoord)*aex/hubble_param,float(ycoord)*aex/hubble_param,float(zcoord)*aex/hubble_param])
gal_vels = ([float(xvel),float(yvel),float(zvel)])

#shiftcoords = np.abs(coords - gal_coords)
#shiftcoords = np.where(shiftcoords > 0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param - shiftcoords, shiftcoords)

shiftcoords = coords - gal_coords
shiftcoords = np.where(shiftcoords > 0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param - shiftcoords, shiftcoords)
shiftcoords = np.where(shiftcoords < -0.5*boxsize*aex/hubble_param, boxsize*aex/hubble_param + shiftcoords, shiftcoords)

shiftvels = velocity - gal_vels

vx_mass = shiftvels[:,0]*mass
vy_mass = shiftvels[:,1]*mass
vz_mass = shiftvels[:,2]*mass

vlimit = 1000

dist = np.sqrt((shiftcoords[:,0])**2+(shiftcoords[:,1])**2+(shiftcoords[:,2])**2)
vdotroverr = shiftcoords[:,0]*velocity[:,0]+shiftcoords[:,1]*velocity[:,1]+shiftcoords[:,2]*velocity[:,2]
vdotroverr /= dist
vdotroverr_mass = vdotroverr*mass
print "shiftcoords[:,0]= ", shiftcoords[:,0]
print "xcoord= ", xcoord
print "dist= ", dist, min(dist), max(dist)
print "shiftvels[:,0]= ", shiftvels[:,0]
print "vdotroverr= ", vdotroverr, min(vdotroverr), max(vdotroverr)

temperature = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Temperature",numThreads=1)
temp_mass = temperature*mass #/nh**(2/3.)*internalenergy/(entropy*3/2.)
print "temp= ", min(temperature), max(temperature), temperature
print "mass= ", min(mass), max(mass), mass

print "temp_mass= ", min(temp_mass), max(temp_mass), temp_mass
#entropy = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Entropy",numThreads=1)
#internalenergy = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/InternalEnergy",numThreads=1)
nh = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Density",numThreads=1)
nh =  nh*unitDensity_in_cgs
nh_mass = nh*mass

hsmooth = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/SmoothingLength",numThreads=1)
hydrogen = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
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
si12_pos = chem_pos[:,68] 
fe17_pos = chem_pos[:,126] 
fe20_pos = chem_pos[:,129] 
fe23_pos = chem_pos[:,132] 
fe25_pos = chem_pos[:,134] 

chem_pos = np.delete(chem_pos, np.s_[::1])
#del chem_pos
#chem_pos = chem_pos[:1]



ndim = 3
nspecies = 16

particle_IDs = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ParticleIDs",numThreads=1)
metallicity = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/Metallicity",numThreads=1)
SFR = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/StarFormationRate",numThreads=1)
oxygen = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/ElementAbundance/Oxygen",numThreads=1)

OnEquationOfState = eagle.readArray("SNAP", sim, input_filename_base_pos, "/PartType0/OnEquationOfState", numThreads=1)

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

    result_Nh = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.hydrogen.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=18.0+np.log10(lgridz/2./lgrid), Vmax=22.0+np.log10(lgridz/2./lgrid),ion='H',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='H')
    print "result_Nh= ", result_Nh
    result_mass = coldens.main(coords, hsmooth, mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.mass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=18.0+np.log10(lgridz/2./lgrid), Vmax=22.0+np.log10(lgridz/2./lgrid),ion='mass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='mass')
    result_tempmass = coldens.main(coords, hsmooth, temp_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.Tmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=22.0, Vmax=29.0,ion='Tmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='T')
    coldens.make_colourmap('coldens_num.%s.T.%s.l%3.1f.png'%(snapname,direction,lgrid), result_tempmass-result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, 4.0,7.0, "T", ngrid, redshift=redshift, extralabel='T',haloinfostr=haloinfostr)
    print "result_tempmass= ", result_tempmass
    result_nhmass = coldens.main(coords, hsmooth, nh_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.nHmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=14.0, Vmax=22.0,ion='nhmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    coldens.make_colourmap('coldens_num.%s.nh.%s.l%3.1f.png'%(snapname,direction,lgrid), result_nhmass-result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -7.0,-1.0, "nH", ngrid, redshift=redshift, extralabel='n$_{\mathrm{H}}$',haloinfostr=haloinfostr)
    print "result_nhmass= ", result_nhmass

    result_vdotroverrmass = coldens.main(coords, hsmooth, vdotroverr_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vdotroverrmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vdotroverrmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    print "result_vdotroverrmass= ", result_vdotroverrmass
    print "10**result_mass= ", 10**result_mass
    coldens.make_colourmap('coldens_num.%s.vdotroverr.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vdotroverrmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vdotroverr", ngrid, redshift=redshift, extralabel='v',haloinfostr=haloinfostr)
    print "result_vdotroverrmass/10**result_mass= ", result_vdotroverrmass/10**result_mass


    result_vxmass = coldens.main(coords, hsmooth, vx_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vxmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vxmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    print "result_vxmass= ", result_vxmass
    print "10**result_mass= ", 10**result_mass
    coldens.make_colourmap('coldens_num.%s.vx.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vxmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vx", ngrid, redshift=redshift, extralabel='v$_{x}$',haloinfostr=haloinfostr)
    print "result_vxmass/10**result_mass= ", result_vxmass/10**result_mass

    result_vymass = coldens.main(coords, hsmooth, vy_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vymass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vymass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    coldens.make_colourmap('coldens_num.%s.vy.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vymass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vy", ngrid, redshift=redshift, extralabel='v$_{y}$',haloinfostr=haloinfostr)

    result_vzmass = coldens.main(coords, hsmooth, vz_mass, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.vzmass.%s.l%3.1f.png'%(snapname,direction,lgrid),Vmin=-vlimit, Vmax=vlimit,ion='vzmass',npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel='')
    coldens.make_colourmap('coldens_num.%s.vz.%s.l%3.1f.png'%(snapname,direction,lgrid), result_vzmass/10**result_mass, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -vlimit,vlimit, "vz", ngrid, redshift=redshift, extralabel='v$_{z}$',haloinfostr=haloinfostr)

    lgridz = lgrid*0.01 #*2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

    for j in xrange(nspecies):
        indexes = np.where(SFR==0)
        if(j==0):
            numerator = mass[indexes]*hydrogen[indexes]*h1_pos[indexes]
            ionname = 'HI'
            ion = 'h1'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 18.0+np.log10(lgridz/2./lgrid)
        if(j==1):
            numerator = mass[indexes]*hydrogen[indexes]*he2_pos[indexes]
            ionname = 'HeII'
            ion = 'he2'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 18.0+np.log10(lgridz/2./lgrid)
#        if(j==1):
#            numerator = mass*hydrogen*(o1_pos+o2_pos+o3_pos+o4_pos+o5_pos+o6_pos+o7_pos+o8_pos+o9_pos)
#            ionname = 'O'
#            ion = 'o'
#            colmin = 13.0
#            colmax = 18.0
        if(j==2):
            indexes = np.where((SFR==0) & (c4_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*c4_pos[indexes]
            ionname = 'CIV'
            ion = 'c4'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==3):
            indexes = np.where((SFR==0) & (n5_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*n5_pos[indexes]
            ionname = 'NV'
            ion = 'n5'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==4):
            indexes = np.where((SFR==0) & (o6_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*o6_pos[indexes]
            ionname = 'OVI'
            ion = 'o6'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==5):
            indexes = np.where((SFR==0) & (o7_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*o7_pos[indexes]
            ionname = 'OVII'
            ion = 'o7'
            colmin = 12.5+np.log10(lgridz/2./lgrid)
            colmax = 16.5+np.log10(lgridz/2./lgrid)
        if(j==6):
            indexes = np.where((SFR==0) & (o8_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*o8_pos[indexes]
            ionname = 'OVIII'
            ion = 'o8'
            colmin = 12.5+np.log10(lgridz/2./lgrid)
            colmax = 16.5+np.log10(lgridz/2./lgrid)
        if(j==7):
            indexes = np.where((SFR==0) & (o9_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*o9_pos[indexes]
            ionname = 'OIX'
            ion = 'o9'
            colmin = 13.0+np.log10(lgridz/2./lgrid)
            colmax = 17.0+np.log10(lgridz/2./lgrid)
        if(j==8):
            indexes = np.where((SFR==0) & (ne8_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*ne8_pos[indexes]
            ionname = 'NeVIII'
            ion = 'ne8'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==9):
            indexes = np.where((SFR==0) & (mg2_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*mg2_pos[indexes]
            ionname = 'MgII'
            ion = 'mg2'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)
        if(j==10):
            indexes = np.where((SFR==0) & (mg10_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*mg10_pos[indexes]
            ionname = 'MgX'
            ion = 'mg10'            
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==11):
            indexes = np.where((SFR==0) & (si12_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*si12_pos[indexes]
            ionname = 'SiXII'
            ion = 'si12'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 15.0+np.log10(lgridz/2./lgrid)
        if(j==12):
            indexes = np.where((SFR==0) & (fe17_pos>1e-10))
            numerator = mass[indexes]*hydrogen[indexes]*fe17_pos[indexes]
            ionname = 'FeXVII'
            ion = 'fe17'
            colmin = 11.0+np.log10(lgridz/2./lgrid)
            colmax = 16.0+np.log10(lgridz/2./lgrid)

        result_num = coldens.main(coords[indexes], hsmooth[indexes], numerator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_num.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)

        f = file('colion_map_out.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid), 'w')
        for ii in xrange(ngrid):
            for jj in xrange(ngrid):
                f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 6.1f\n'%(ii,jj,result_num[ii,jj],result_Nh[ii,jj],result_nhmass[ii,jj]-result_mass[ii,jj],result_tempmass[ii,jj]-result_mass[ii,jj],result_vdotroverrmass[ii,jj]/10**result_mass[ii,jj]))
        f.close()

