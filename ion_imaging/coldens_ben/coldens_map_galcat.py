#import numpy
#import eagle
import sys
#sys.path.insert(1,"/home/oppenheimer/.local/lib/python2.7/site-packages/")
import eagle
#import os
import coldens
import numpy as np
import h5py
import tables


#sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
#tag= sys.argv[1]  # Changed on 11/13/14.  
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
snip = sys.argv[3]
xcoord = float(sys.argv[4])
ycoord = float(sys.argv[5])
zcoord = float(sys.argv[6])
lgrid = float(sys.argv[7])
ngrid = int(sys.argv[8])
runlabel = sys.argv[9]
docbar = sys.argv[10]

if (runlabel=='halo'):
    haloname = snapname.split("halo")[1]
    halo_id = haloname.split("_")[1]
    mh = float(haloname.split("_")[2])
    ms = float(haloname.split("_")[3])
    sfr = float(haloname.split("_")[4])
    haloinfostr = 'lg M$_{200}=%4.1f$, lg M$_{*}=%3.1f$, SFR$=%4.2f$'%(mh,ms,sfr)
    #haloinfostr = 'lg M$_{200}=' + mh + '$, lg M$_{*}=' + ms + '$, SFR$=' + sfr + '$'
else:
    haloinfostr = runlabel

center = np.array([xcoord, ycoord, zcoord])

lgrid = lgrid/1e+03
#lgrid_half = lgrid/2. 
lgridz = lgrid*4 # Four times longer than Lgrid in zdirection.

if(snip=='0'):
    path = "snapshot_%s/snap_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNAP"
else:
    path = "snipshot_%s/snip_%s.0.hdf5"%(input_filename_base,input_filename_base)
    sniptag = "SNIP"
#filein = h5py.File(path)
#redshift = filein['Header'].attrs['Redshift']
#aex = 1/(1+redshift)
#center = center*aex

aex = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/HubbleParam")
redshift = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/Redshift")
boxsize = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/BoxSize")

boxsize = boxsize/hubble_param*aex
print "boxsize=", boxsize 
center = center/hubble_param*aex
print "center= ", center
coords = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Coordinates",numThreads=1)
mass = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Mass",numThreads=1)

if(snip=='1'): # Have to back out hsmooth from density
    density = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/Density",numThreads=1)
    print "density= ", density
    hsmooth = (density/mass)**(-0.3333)*2.39  #2.39 conversion factor for 58 neighbors?  
    print "hsmooth= ",hsmooth
else:
    hsmooth = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/SmoothingLength",numThreads=1)

mass *= 1.e+10
print "mass= ", mass


temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1,physicalUnits=0,noH=0)
SFR = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/StarFormationRate",numThreads=1)
indexes_cool = np.where((temperature<1e+05) & (SFR == 0))

hydrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
oxygen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
silicon = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Silicon",numThreads=1)

snapname0 = sim + "/snapshot_" + input_filename_base + "/snap_" + input_filename_base + ".0.hdf5" 
print "snapname0= ", snapname0
h50file = h5py.File(snapname0,"r")

ChemicalAbundances = "/PartType0/ChemicalAbundances/HydrogenI" in h50file
print "ChemicalAbundances= ", ChemicalAbundances


if(ChemicalAbundances): #(snip==1) | (snip==0)):
    h1 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1,physicalUnits=1,noH=0)
    print 'h1= ', h1
    c2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonII",numThreads=1,physicalUnits=1,noH=0)
    c3 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIII",numThreads=1,physicalUnits=1,noH=0)
    c4 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1,physicalUnits=1,noH=0)
    print 'c4= ', c4
    o6 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1,physicalUnits=1,noH=0)
    o7 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1,physicalUnits=1,noH=0)
    o8 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1,physicalUnits=1,noH=0)
    mg2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1,physicalUnits=1,noH=0)
    si2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconII",numThreads=1,physicalUnits=1,noH=0)
    si3 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIII",numThreads=1,physicalUnits=1,noH=0)
    si4 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIV",numThreads=1,physicalUnits=1,noH=0)
    ne8 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/NeonVIII",numThreads=1,physicalUnits=1,noH=0)
    mg10 = si12 = n5 = o1 = o2 = o3 = o4 = o5 = c6 = o8*0.0
else:
    print "snip = ", snip
    if(snip=="0"):
        print "reading ChemistryAbundanecs Array"
        chem = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
        h1 = chem[:,1]
        c4 = chem[:,10]
        o6 = chem[:,28]
        o7 = chem[:,29]
        o8 = chem[:,30]
        mg2 = chem[:,45] 
        ne8 = chem[:,40]
        mg10 = chem[:,53]
        si12 = chem[:,68]
        n5 = chem[:,19]
        o1 = chem[:,23]
        o2 = chem[:,24]
        o3 = chem[:,25]
        o4 = chem[:,26]
        o5 = chem[:,27]
        c2 = chem[:,8]
        c3 = chem[:,9]
        c6 = chem[:,12]
        si2 = chem[:,58]
        si3 = chem[:,59]
        si4 = chem[:,60]
    else:        
        h1 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
        print 'h1= ', h1
        c4 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
        print 'c4= ', c4
        o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
        o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
        o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
        mg2 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)
        #ne8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/NeonVIII",numThreads=1)
        mg10 = si12 = n5 = o1 = o2 = o3 = o4 = o5 = o8*0.0
        ne8 = o8*0.0

print "sim = ", sim
print "input_filename_base = ", input_filename_base
group_filename_base = input_filename_base.split("_")[1] + "_" + input_filename_base.split("_")[2]
print "group_filename_base = ", group_filename_base


M_200 = eagle.readArray("SUBFIND_GROUP", sim, group_filename_base, "/FOF/Group_M_Crit200",numThreads=1) * 1e10
R_200 = eagle.readArray("SUBFIND_GROUP", sim, group_filename_base, "/FOF/Group_R_Crit200",numThreads=1) 

print "M_200 = ", M_200
print "R_200 = ", R_200

M_halo = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/Mass",numThreads=1) * 1e+10
CenterofPotential = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/CentreOfPotential",numThreads=1)
Velocity = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/Velocity",numThreads=1) 
GrpID = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/GroupNumber",numThreads=1)
SubID = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/SubGroupNumber",numThreads=1)
M_Star = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/Stars/Mass",numThreads=1) * 1e10
M_SF = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/SF/Mass",numThreads=1) * 1e10
M_NSF = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/NSF/Mass",numThreads=1) * 1e10
SFR = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/StarFormationRate",numThreads=1)
M_star_30kpc = eagle.readArray("SUBFIND", sim, group_filename_base, "/Subhalo/ApertureMeasurements/Mass/030kpc",numThreads=1)[:,4] * 1e10
SFR_30kpc = eagle.readArray("SUBFIND", sim, group_filename_base, "/Subhalo/ApertureMeasurements/SFR/030kpc",numThreads=1)

MassStarLimit = 1e+08

h5file_out = tables.openFile("galcat.%s.l%3.1f.hdf5"%(snapname,lgrid), "w")
fout = open("galcat.%s.l%3.1f.dat"%(snapname,lgrid),"w")

xcoord_ = []
ycoord_ = []
zcoord_ = []
xvel_ = []
yvel_ = []
zvel_ = []
xpix_ = []
ypix_ = []
zpix_ = []
GrpID_ = []
lM200_ = []
lMstar_ = []
SFR_ = []
lMhalo_ = []

for i in range(len(M_halo)):
    if(M_star_30kpc[i] > MassStarLimit):
        GrpID_ = np.append(GrpID_,GrpID[i])
        lM200_ = np.append(lM200_,np.log10(M_200[GrpID[i]-1]))
        lMstar_ = np.append(lMstar_,np.log10(M_star_30kpc[i]))
        SFR_ = np.append(SFR_,SFR_30kpc[i])
        lMhalo_ = np.append(lMhalo_, np.log10(M_halo[i]))
        xcoord_ = np.append(xcoord_,CenterofPotential[i,0]*hubble_param/aex)
        ycoord_ = np.append(ycoord_,CenterofPotential[i,1]*hubble_param/aex)
        zcoord_ = np.append(zcoord_,CenterofPotential[i,2]*hubble_param/aex)
        xvel_ = np.append(xvel_,Velocity[i,0]*hubble_param/aex)
        yvel_ = np.append(yvel_,Velocity[i,1]*hubble_param/aex)
        zvel_ = np.append(zvel_,Velocity[i,2]*hubble_param/aex)
        xpix_ = np.append(xpix_,ngrid*((CenterofPotential[i,0]*hubble_param/aex-xcoord)/aex)/lgrid+ngrid*0.5)
        ypix_ = np.append(ypix_,ngrid*((CenterofPotential[i,1]*hubble_param/aex-ycoord)/aex)/lgrid+ngrid*0.5)
        zpix_ = np.append(zpix_,ngrid*((CenterofPotential[i,2]*hubble_param/aex-zcoord)/aex)/lgrid+ngrid*0.5)
        ###xpix_ = np.append(xpix_,ngrid*((CenterofPotential[i,0]*hubble_param/aex-xcoord)/aex+1)/lgrid)
        ###ypix_ = np.append(ypix_,ngrid*((CenterofPotential[i,1]*hubble_param/aex-ycoord)/aex+1)/lgrid)
        ###zpix_ = np.append(zpix_,ngrid*((CenterofPotential[i,2]*hubble_param/aex-zcoord)/aex+1)/lgrid)

        #xpix_ = ngrid*((CenterofPotential[i,0]*hubble_param/aex-xcoord)/aex+1)/lgrid
        #ypix_ = ngrid*((CenterofPotential[i,1]*hubble_param/aex-ycoord)/aex+1)/lgrid
        #zpix_ = ngrid*((CenterofPotential[i,2]*hubble_param/aex-zcoord)/aex+1)/lgrid
        fout.write(" %3d %5.2f %5.2f %7.3f %5.2f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f % 10.6f %6.1f %6.1f %6.1f\n"%(GrpID[i], np.log10(M_200[GrpID[i]-1]), np.log10(M_star_30kpc[i]), SFR_30kpc[i], np.log10(M_halo[i]), CenterofPotential[i,0]*hubble_param/aex, CenterofPotential[i,1]*hubble_param/aex, CenterofPotential[i,2]*hubble_param/aex, xcoord, ycoord, zcoord, (CenterofPotential[i,0]*hubble_param/aex-xcoord)/aex/lgrid*ngrid+ngrid*0.5, (CenterofPotential[i,1]*hubble_param/aex-ycoord)/aex/lgrid*ngrid+ngrid*0.5, (CenterofPotential[i,2]*hubble_param/aex-zcoord)/aex/lgrid*ngrid+ngrid*0.5))

print "xpix_ = ", xpix_

GalCat_group = h5file_out.create_group(h5file_out.root, "GalCat", "Galaxy Catalogue")
#h5file.create_array(gcolumns, 'pressure', array(pressure), "Pressure column selection")

h5file_out.createArray(h5file_out.root.GalCat, "GrpID", np.array(GrpID_, dtype = np.int))
h5file_out.createArray(h5file_out.root.GalCat, "lM_200", np.array(lM200_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "lMstar", np.array(lMstar_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "SFR", np.array(SFR_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "lMhalo", np.array(lMhalo_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "x_pixel", np.array(xpix_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "y_pixel", np.array(ypix_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "z_pixel", np.array(zpix_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "x_coord", np.array(xcoord_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "y_coord", np.array(ycoord_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "z_coord", np.array(zcoord_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "x_vel", np.array(xvel_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "y_vel", np.array(yvel_, dtype = np.float32))
h5file_out.createArray(h5file_out.root.GalCat, "z_vel", np.array(zvel_, dtype = np.float32))



print "lgrid, ngrid = ", lgrid, ngrid

fout.close()

dist = np.zeros((ngrid,ngrid))
xval = np.zeros((ngrid,ngrid))
yval = np.zeros((ngrid,ngrid))
zval = np.zeros((ngrid,ngrid))


ndim = 3
nions = 11

for k in range(ndim):
    if(k==2):
        thetaangle = 0
        phiangle = 0
        for i in xrange(ngrid):
            for j in xrange(ngrid):
                xval[i,j] = i
                yval[i,j] = j
                dist[i,j] = np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03
                
        Zmap_group = h5file_out.create_group(h5file_out.root, "Zmap", "Mapping in z direction")
            
        atom = tables.Atom.from_dtype(dist.dtype)
        colmap = h5file_out.createCArray(h5file_out.root.Zmap, 'Distance', atom, dist.shape)
        colmap[:] = dist
        colmap = h5file_out.createCArray(h5file_out.root.Zmap, 'x-pixel', atom, dist.shape)
        colmap[:] = xval
        colmap = h5file_out.createCArray(h5file_out.root.Zmap, 'y-pixel', atom, dist.shape)
        colmap[:] = yval
        direction = 'z'

    if(k==1):
        thetaangle = -90
        phiangle = 0
        for i in xrange(ngrid):
            for j in xrange(ngrid):
                xval[i,j] = i
                zval[i,j] = j
                dist[i,j] = np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03
                
        Ymap_group = h5file_out.create_group(h5file_out.root, "Ymap", "Mapping in y direction")
            
        atom = tables.Atom.from_dtype(dist.dtype)
        colmap = h5file_out.createCArray(h5file_out.root.Ymap, 'Distance', atom, dist.shape)
        colmap[:] = dist
        colmap = h5file_out.createCArray(h5file_out.root.Ymap, 'x-pixel', atom, dist.shape)
        colmap[:] = xval
        colmap = h5file_out.createCArray(h5file_out.root.Ymap, 'z-pixel', atom, dist.shape)
        colmap[:] = zval
        direction = 'y'

    if(k==0):
        thetaangle = 0
        phiangle = 90
        for i in xrange(ngrid):
            for j in xrange(ngrid):
                zval[i,j] = i
                yval[i,j] = j
                dist[i,j] = np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03
                
        Xmap_group = h5file_out.create_group(h5file_out.root, "Xmap", "Mapping in x direction")
            
        atom = tables.Atom.from_dtype(dist.dtype)
        colmap = h5file_out.createCArray(h5file_out.root.Xmap, 'Distance', atom, dist.shape)
        colmap[:] = dist
        colmap = h5file_out.createCArray(h5file_out.root.Xmap, 'z-pixel', atom, dist.shape)
        colmap[:] = zval
        colmap = h5file_out.createCArray(h5file_out.root.Xmap, 'y-pixel', atom, dist.shape)
        colmap[:] = yval
        direction = 'x'


    for l in range(nions):

        if(l==0):        
            indexes = np.where(mass>0.0)
            ionarray = mass[indexes]*hydrogen[indexes]
            ionname = 'hydrogen'
            ionlabel = 'H'
            colmin = 18.5
            colmax = 21.5
        if(l==1):
            indexes = np.where(oxygen>1e-08)
            ionarray = mass[indexes]*oxygen[indexes]/16.0
            ionname = 'oxygen'
            ionlabel = 'O'
            colmin = 14
            colmax = 18
        if(l==2):        
            indexes = np.where(h1>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*h1[indexes]
            ionname = 'h1'
            ionlabel = 'HI'
            colmin = 13
            colmax = 21
        if(l==3):        
            indexes = np.where(c4>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*c4[indexes]
            ionname = 'c4'
            ionlabel = 'CIV'
            colmin = 11.5
            colmax = 15.5
        if(l==4):        
            indexes = np.where(o6>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*o6[indexes]
            ionname = 'o6'
            ionlabel = 'OVI'
            colmin = 11.5
            colmax = 15.5
        if(l==5):        
            indexes = np.where(o7>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*o7[indexes]
            ionname = 'o7'
            ionlabel = 'OVII'
            colmin = 13.0
            colmax = 17.0
        if(l==6):        
            indexes = np.where(o8>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*o8[indexes]
            ionname = 'o8'
            ionlabel = 'OVIII'
            colmin = 13.0
            colmax = 17.0
        if(l==7):        
            indexes = np.where(mg2>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*mg2[indexes]
            ionname = 'mg2'
            ionlabel = 'MgII'
            colmin = 11.0
            colmax = 15.0
        if(l==8):        
            indexes = np.where(si2>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*si2[indexes]
            ionname = 'si2'
            ionlabel = 'SiII'
            colmin = 11.0
            colmax = 15.0
        if(l==9):        
            indexes = np.where(si3>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*si3[indexes]
            ionname = 'si3'
            ionlabel = 'SiIII'
            colmin = 11.0
            colmax = 15.0
        if(l==10):        
            indexes = np.where(si4>1e-08)
            ionarray = mass[indexes]*hydrogen[indexes]*si4[indexes]
            ionname = 'si4'
            ionlabel = 'SiIV'
            colmin = 11.0
            colmax = 15.0


        result_ion = coldens.main(coords[indexes], hsmooth[indexes], ionarray, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.%s.l%3.1f.png'%(snapname,ionname,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ionlabel,redshift=redshift,theta=thetaangle,phi=phiangle,extralabel=ionlabel,haloinfostr=haloinfostr,docbar=docbar)


        atom = tables.Atom.from_dtype(result_ion.dtype)
        if(k==0):
            colmap = h5file_out.createCArray(h5file_out.root.Xmap, ionlabel, atom, result_ion.shape)
        if(k==1):
            colmap = h5file_out.createCArray(h5file_out.root.Ymap, ionlabel, atom, result_ion.shape)
        if(k==2):
            colmap = h5file_out.createCArray(h5file_out.root.Zmap, ionlabel, atom, result_ion.shape)
        colmap[:] = result_ion


h5file_out.close()

