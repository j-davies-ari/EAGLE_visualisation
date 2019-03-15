#import numpy
import eagle 
import sys
#import os
import coldens
import numpy as np
import h5py
import tables
import ctypes as ct
import math

def calc_radial(ion_array):

    nbins = 100
    dist = ion_array*0.0

    for i in range(0,ngrid-1):
        for j in range(0,ngrid-1):
            xcoord = abs(ngrid/2.-0.5-i)/ngrid*lgrid*1e+03
            ycoord = abs(ngrid/2.-0.5-j)/ngrid*lgrid*1e+03
            dist[i,j] = np.sqrt(xcoord**2+ycoord**2)


    npart,bins = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.])
    #print npart
    #print bins
    
    hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.], weights=ion_array)[0]
    dist_bin = hist*0.0
    dbin = hist*0.0
    for i in xrange(nbins):
        hist[i] = hist[i]/npart[i]
        dist_bin[i] = (bins[i] + bins[i+1])/2.
        dbin[i] = bins[i+1]-bins[i]
        
    return(hist,dist_bin,dbin)
    

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

unit_mass_in_cgs = 1.989e33 * 1.0e10 
unit_length_in_cgs = 3.0857e24 
proton_mass_cgs = 1.67e-24
unitDensity_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3


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

pixsize_kpc = lgrid/ngrid

lgrid = lgrid/1e+03
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

solarabundance_carbon = eagle.readAttribute(sniptag, sim, input_filename_base, "/Parameters/ChemicalElements/SolarAbundance_Carbon")
solarabundance_oxygen = eagle.readAttribute(sniptag, sim, input_filename_base, "/Parameters/ChemicalElements/SolarAbundance_Oxygen")
solarabundance_hydrogen = eagle.readAttribute(sniptag, sim, input_filename_base, "/Parameters/ChemicalElements/SolarAbundance_Hydrogen")
solarabundance_magnesium = eagle.readAttribute(sniptag, sim, input_filename_base, "/Parameters/ChemicalElements/SolarAbundance_Magnesium")


if(redshift < 0.045):
    redshift_str = "0.0000"
    cm2arcsec2 = (44.6*1e+03)**2*4*np.pi/0.212**2 # 1.80e-12
    kpcarcsec = 0.212
    redshift = 0.01 #SHIFT TO z=0.01.
else:
    if(redshift < 0.075):
        redshift_str = "0.0000"
        cm2arcsec2 = (229.5*1e+03)**2*4*np.pi/1.009**2 #1.38e-12
        kpcarcsec = 1.009
    else:
        if(redshift < 0.185):
            redshift_str = "0.1006"
            cm2arcsec2 = (457*1e+03)**2*4*np.pi/1.904**2 #1.38e-12
            kpcarcsec = 1.904
        else:
            if(redshift < 0.335):
                redshift_str = "0.2709"
                cm2arcsec2 = (1300*1e+03)**2*4*np.pi/4.304**2 #8.72e-13 
                kpcarcsec = 4.304
            else:
                if((redshift > 0.9) | (redshift<1.2)):
                    redshift_str = "0.9567"
                    cm2arcsec2 = (6796.4*1e+03)**2*4*np.pi/8.238**2
                    kpcarcsec = 8.238
                else:                                        
                    redshift_str = "0.4675"

pixsize_arcsec = pixsize_kpc/kpcarcsec
npix_5arcmin = int(300./pixsize_arcsec)
print "npix_5arcmin= ", npix_5arcmin
            
boxsize = boxsize/hubble_param*aex
print "boxsize=", boxsize 
center = center/hubble_param*aex
print "center= ", center
coords = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Coordinates",numThreads=1)
mass = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Mass",numThreads=1)*unit_mass_in_cgs
density = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Density",numThreads=1)*unitDensity_in_cgs
metallicity = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Metallicity",numThreads=1)

if(snip=='1'): # Have to back out hsmooth from density
    print "density= ", density
    hsmooth = (density/mass)**(-0.3333)*2.39  #2.39 conversion factor for 58 neighbors?  
    print "hsmooth= ",hsmooth
else:
    hsmooth = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/SmoothingLength",numThreads=1)

#mass *= unit_mass_in_cgs
print "mass= ", mass


temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1,physicalUnits=0,noH=0)
SFR = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/StarFormationRate",numThreads=1)
indexes_cool = np.where((temperature<1e+05) & (SFR == 0))

hydrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
helium = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Helium",numThreads=1)
carbon = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Carbon",numThreads=1)
nitrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Nitrogen",numThreads=1)
oxygen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
neon = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Neon",numThreads=1)
magnesium = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Magnesium",numThreads=1)
silicon = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Silicon",numThreads=1)
iron = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Iron",numThreads=1)

nh = density*hydrogen/proton_mass_cgs

calcEmissionLib = ct.cdll.LoadLibrary("/home/oppenheimer/pythonCode/calcEmission/lookupEmission/lookupEmission.so")
#calcEmissionLib = ct.cdll.LoadLibrary("/home/beop5934/pythonCode/calcEmission/lookupEmission/lookupEmission.so")
calcEmissionLib.return_lineEmission.restype = ct.c_double
calcEmissionLib.return_contEmission.restype = ct.c_double


#if(snip=="1"):
#    h1 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
#    c4 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
#    o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
#    o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
#    o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
#    mg2 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)
#    mg10 = si12 = n5 = o1 = o2 = o3 = o4 = o5 = o8*0.0
#    ne8 = o8*0.0
#else:
#    chem = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
#    h1 = chem[:,1]
#    c4 = chem[:,10]
#    o6 = chem[:,28]
#    o7 = chem[:,29]
#    o8 = chem[:,30]
#    mg2 = chem[:,45] 
#    ne8 = chem[:,40]
#    mg10 = chem[:,53]
#    si12 = chem[:,68]
#    n5 = chem[:,19]
#    o1 = chem[:,23]
#    o2 = chem[:,24]
#    o3 = chem[:,25]
#    o4 = chem[:,26]
#    o5 = chem[:,27]
#    c2 = chem[:,8]
#    c3 = chem[:,9]
#    c6 = chem[:,12]
#    si2 = chem[:,58]
#    si3 = chem[:,59]
#    si4 = chem[:,60]
    

print "sim = ", sim
print "input_filename_base = ", input_filename_base
if (input_filename_base.split("_")[0] == "noneq"): 
    group_filename_base = input_filename_base.split("_")[1] + "_" + input_filename_base.split("_")[2]
else:
    group_filename_base = input_filename_base.split("_")[0] + "_" + input_filename_base.split("_")[1]

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
SFR_sub = eagle.readArray("SUBFIND", sim, group_filename_base, "Subhalo/StarFormationRate",numThreads=1)
M_star_30kpc = eagle.readArray("SUBFIND", sim, group_filename_base, "/Subhalo/ApertureMeasurements/Mass/030kpc",numThreads=1)[:,4] * 1e10
SFR_30kpc = eagle.readArray("SUBFIND", sim, group_filename_base, "/Subhalo/ApertureMeasurements/SFR/030kpc",numThreads=1)

MassStarLimit = 1e+08

h5file_out = tables.openFile("emission.%s.l%3.1f.hdf5"%(snapname,lgrid), "w")

xcoord_ = []
ycoord_ = []
zcoord_ = []
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
        xpix_ = np.append(xpix_,ngrid*((CenterofPotential[i,0]*hubble_param/aex-xcoord)/aex+1)/lgrid)
        ypix_ = np.append(ypix_,ngrid*((CenterofPotential[i,1]*hubble_param/aex-ycoord)/aex+1)/lgrid)
        zpix_ = np.append(zpix_,ngrid*((CenterofPotential[i,2]*hubble_param/aex-zcoord)/aex+1)/lgrid)

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


FiveArcMin_group = h5file_out.create_group(h5file_out.root, "FiveArcMinSum", "Sums within 5 arcmin pixel")
Impactbin_group = h5file_out.create_group(h5file_out.root, "ImpactBins", "Impact parameter binning in kpc bins")

dist = np.zeros((ngrid,ngrid))
xval = np.zeros((ngrid,ngrid))
yval = np.zeros((ngrid,ngrid))
zval = np.zeros((ngrid,ngrid))

ndim = 3 #1
nlines = 11 #11

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


    for l in range(nlines):
        emmin = -26.0
        emmax = -16.0
        if (l==0):
            indexes = np.where((mass>0.0) & (SFR==0) & (temperature>=10**6.5)) #5.3))
            calcEmissionLib.load_contEmissionTable_redshift("bolometric") 
            emissionarray = mass[indexes]*hydrogen[indexes]*0.0+1e-95
            linelabel = "Bolometric" 
            linename = "Bolometric" 
            #emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs
            emission_multiple = nh[indexes]**2*mass[indexes]/density[indexes]*(hydrogen[indexes]+2./4.*(1-hydrogen[indexes]))/hydrogen[indexes]
        if (l==1):
            indexes = np.where((mass>0.0) & (SFR==0) & (temperature>=10**5.3))
            print "START loading calcEmissionLib.load_contEmissionTable_redshift"
            calcEmissionLib.load_contEmissionTable_redshift("0.5_2.0keV") 
            emissionarray = mass[indexes]*hydrogen[indexes]*0.0+1e-95
            linelabel = "SoftXray" #"0.5-2.0keV"
            linename = "SoftXray"
            #emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs
            emission_multiple = nh[indexes]**2*mass[indexes]/density[indexes]*(hydrogen[indexes]+2./4.*(1-hydrogen[indexes]))/hydrogen[indexes]
        if (l==2):
            indexes = np.where((mass>0.0) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"hydrogen",35,1) 
            emissionarray = mass[indexes]*hydrogen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "Lya" #"Ly$\alpha$"
            linename = "Lya"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*hydrogen[indexes]/solarabundance_hydrogen
        if (l==3):
            indexes = np.where((mass>0.0) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"hydrogen",35,9) 
            emissionarray = mass[indexes]*hydrogen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "Ha" #"H$\alpha$"
            linename = "Ha"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*hydrogen[indexes]/solarabundance_hydrogen
        if (l==4):
            indexes = np.where((oxygen>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"oxygen",151,119) #OVI1032
            emissionarray = mass[indexes]*oxygen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "OVI_1032" 
            linename = "OVI_1032"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*oxygen[indexes]/solarabundance_oxygen
        if (l==5):
            indexes = np.where((oxygen>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"oxygen",151,136) #OVII22.1
            emissionarray = mass[indexes]*oxygen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "OVII_21.60" 
            linename = "OVII_21.60"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*oxygen[indexes]/solarabundance_oxygen
        if (l==6):
            indexes = np.where((oxygen>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"oxygen",151,138) #OVII22.1
            emissionarray = mass[indexes]*oxygen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "OVII_21.81" 
            linename = "OVII_21.81"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*oxygen[indexes]/solarabundance_oxygen
        if (l==7):
            indexes = np.where((oxygen>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"oxygen",151,139) #OVII22.1
            emissionarray = mass[indexes]*oxygen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "OVII_22.10" 
            linename = "OVII_22.10"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*oxygen[indexes]/solarabundance_oxygen
        if (l==8):
            indexes = np.where((oxygen>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"oxygen",151,149) #OVIII18.97
            emissionarray = mass[indexes]*oxygen[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "OVIII_18.98" 
            linename = "OVIII_18.98"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*oxygen[indexes]/solarabundance_oxygen
        if (l==9):
            indexes = np.where((carbon>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"carbon",57,31) #CIII1907
            emissionarray = mass[indexes]*carbon[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "CIII" 
            linename = "CIII"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*carbon[indexes]/solarabundance_carbon
        if (l==10):
            indexes = np.where((magnesium>1e-08) & (SFR==0))
            calcEmissionLib.load_lineEmissionTable_redshift(redshift_str,"magnesium",184,1) #MgII-2796
            emissionarray = mass[indexes]*carbon[indexes]*0.0+1e-95 # Fix. Msolar to cm^-2
            linelabel = "MgII" 
            linename = "MgII"
            emission_multiple = nh[indexes]*mass[indexes]*1.00794/proton_mass_cgs*magnesium[indexes]/solarabundance_magnesium


        H = hydrogen[indexes]
        He = helium[indexes]
        C = carbon[indexes]
        N = nitrogen[indexes]
        O = oxygen[indexes]
        Ne = neon[indexes]
        Mg = magnesium[indexes]
        Si = silicon[indexes]
        Fe = iron[indexes]         

        print "Indexes for metals done."
   

        dens = np.log10(nh[indexes])
        T = np.log10(temperature[indexes])

        emission_ergs = np.log10(nh[indexes])*0.0
        emission_sum = np.float128(0.0)
        emission_ergs = np.ndarray((len(nh[indexes]),),np.float128)

        for i in range(len(dens)):
            if((l==0) | (l==1)):
                emissionarray[i] = calcEmissionLib.return_contEmission(ct.c_double(T[i]),ct.c_double(H[i]),ct.c_double(He[i]),ct.c_double(C[i]),ct.c_double(N[i]),ct.c_double(O[i]),ct.c_double(Ne[i]),ct.c_double(Mg[i]),ct.c_double(Si[i]),ct.c_double(Fe[i]))
            else:
                emissionarray[i] = calcEmissionLib.return_lineEmission(ct.c_float(dens[i]),ct.c_float(T[i]),ct.c_float(redshift))
                #if(emissionarray[i]>0): emissionarray[i] = -95.0 # problem disaster not yet understood!!! 

            emission_ergs[i] = emissionarray[i] + np.log10(emission_multiple[i])
            emission_ergs[i] = 10**emission_ergs[i]

            if((i%100000==0) | (math.isinf(emission_ergs[i]))): print "em = ", emissionarray[i], 10**emissionarray[i] *emission_multiple[i],T[i],dens[i],i,len(dens),emission_multiple[i], emission_ergs[i], emission_sum #* dens[i] * mass[i] *1.00794 /proton_mass_cgs*hydrogen[i]/solarabundance_hydrogen, T[i], dens[i], mass[i], hydrogen[i], i, len(dens)
            if(math.isinf(emission_ergs[i])): # Used to be a problem.  
                print "em inf= ",emission_ergs[i],i,T[i],dens[i], " em_mult= ", emission_multiple[i], " em_array= ", emissionarray[i], "emission_erg,1,2,3= ", emissionarray[i] + np.log10(emission_multiple[i]), 10**emissionarray[i]*emission_multiple[i], 10**(emissionarray[i] + np.log10(emission_multiple[i]))
            emission_sum += np.float128(emission_ergs[i])
        
        print emissionarray, emission_ergs
        print "emission_ergs ", linelabel, " min,max,med= ", np.min(emission_ergs),np.max(emission_ergs), np.median(emission_ergs)
        print "emission_ergs ", linelabel, "sum= ",emission_sum #np.sum(emission_ergs[np.where(math.isinf(10**emissionarray))])

        emission_ergscmarcsec = emission_ergs/cm2arcsec2
        
        #m_H = 1.00794
        #n_H = gasrho * gasX / protonmass

            #lineem = 10^lineem * n_H^2 * gasmass * m_H / (gasrho * gasX) ;erg/s per particle

            #lineem =  10^lineem * gasrho^2 * gasX^2 /(gasrho * gasX) * m_H/protonmass^2 *gasmass
            #lineem = 10^lineem * gasrho * gasX / protonmass * m_H/protonmass * gasmass
            #lineem = 10^lineem * nH * 1.00794* gasmass /protonmass
            #lineem = lineem * gasAbund / Element_solar
            #; What to do with starforming gas?
            #wheresfr = where(gassfr NE 0.)
            #lineem[wheresfr] = min(lineem)/1d10;0.1*1.65*1d42 * gassfr[wheresfr] ;erg/s per particle

            #print, alog10(minmax(lineem))
            #IF NOT keyword_set(ergs) THEN BEGIN

            #rcomoving = comovingdistance(rshift, omega0, omegalambda0, hubble0)
            #lumdistance = double(rcomoving * (1. + rshift) * cm_per_mpc)

            #lineem = lineem / (4. * !pi * lumdistance^2.) ;erg/cm^2/s

            #lineem = lineem * 1216e-8 / (1. + rshift) / (planck * lightspeed) ; photons/cm^2./s per particle
            
            
        result_ion = coldens.main(coords[indexes], hsmooth[indexes], emission_ergscmarcsec, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='emission.%s.%s.%s.l%3.1f.png'%(snapname,linename,direction,lgrid),Vmin=emmin, Vmax=emmax,ion=linelabel,redshift=redshift,theta=thetaangle,phi=phiangle,extralabel=linelabel,haloinfostr=haloinfostr,docbar=docbar,emission=1)

        ion_ave, rbins, dbin = calc_radial(result_ion)
        print "ion_ave = ",ion_ave
        print "rbins = ", rbins
        #print "dbin = ", dbin
        if(k==0):
            if(l==0):
                h5file_out.createArray(h5file_out.root.ImpactBins, "bbins", np.array(rbins,dtype=np.float32))
                h5file_out.createArray(h5file_out.root.ImpactBins, "deltabins", np.array(dbin,dtype=np.float32))
            h5file_out.createArray(h5file_out.root.ImpactBins, "X_bins_" + linelabel, np.array(ion_ave,dtype=np.float32))
        if(k==1):
            h5file_out.createArray(h5file_out.root.ImpactBins, "Y_bins_" + linelabel, np.array(ion_ave,dtype=np.float32))
        if(k==2):
            h5file_out.createArray(h5file_out.root.ImpactBins, "Z_bins_" + linelabel, np.array(ion_ave,dtype=np.float32))
            
            
        sum_5arcmin = 0.
        Nsum_5arcmin = 0
        for i in xrange(ngrid):
            for j in xrange(ngrid):
                if((ngrid/2.-0.5-j<npix_5arcmin/2.) & (j-(ngrid/2.-0.5)<npix_5arcmin/2.)):
                    if((ngrid/2.-0.5-i<npix_5arcmin/2.) & (i-(ngrid/2.-0.5)<npix_5arcmin/2.)):
                        sum_5arcmin += 10**(result_ion[i,j])
                        Nsum_5arcmin += 1 

        print "5arcmin_sum= ", sum_5arcmin, Nsum_5arcmin, sum_5arcmin/Nsum_5arcmin

        atom = tables.Atom.from_dtype(result_ion.dtype)
        if(k==0):
            colmap = h5file_out.createCArray(h5file_out.root.Xmap, linelabel, atom, result_ion.shape)
            if(l==0):
                h5file_out.createArray(h5file_out.root.FiveArcMinSum, "npix", np.array(Nsum_5arcmin, dtype=np.int))
            h5file_out.createArray(h5file_out.root.FiveArcMinSum, "Xave_" + linelabel, np.array(sum_5arcmin/Nsum_5arcmin, dtype=np.float32))
        if(k==1):
            colmap = h5file_out.createCArray(h5file_out.root.Ymap, linelabel, atom, result_ion.shape)
            h5file_out.createArray(h5file_out.root.FiveArcMinSum, "Yave_" + linelabel, np.array(sum_5arcmin/Nsum_5arcmin, dtype=np.float32))
        if(k==2):
            colmap = h5file_out.createCArray(h5file_out.root.Zmap, linelabel, atom, result_ion.shape)
            h5file_out.createArray(h5file_out.root.FiveArcMinSum, "Zave_" + linelabel, np.array(sum_5arcmin/Nsum_5arcmin, dtype=np.float32))
        colmap[:] = result_ion


        
h5file_out.close()

