#import numpy
import socket
if(socket.gethostname()=="quasar.strw.leidenuniv.nl"): 
    import eagle_ben as eagle
else:
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
#tag= sys.argv[1]  # Changed on 11/13/14.  
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
snip = int(sys.argv[3])
ionshort = sys.argv[4] # "h1_o6"
envirocat = sys.argv[5]
ngrid = int(sys.argv[6])
nslice = int(sys.argv[7])
print "ngrid= ", ngrid
#runlabel = sys.argv[8]
docbar = sys.argv[8]

ionshort = "h1_o6"

do_metals = 1

#center = np.array([6.98,5.21,6.55])
#center= np.array([17.5995,14.08347,15.8329])  #snapshot 31
#center= np.array([16.3672,13.0542,14.6979])  #z=0.271
#center = np.array([15.2534, 10.9404,  9.0412])
#center = np.array([15.2691, 10.934, 9.03164])
#center = np.array([15.2974, 10.9540,  9.0412])

if(snip!=1):
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
center = np.array([boxsize/2.,boxsize/2.,boxsize/2.])
print "center= ", center

lgrid = boxsize
lgridz = lgrid*1./nslice

print "lgrid, lgridz= ", lgrid, lgridz


#lgrid = lgrid/1e+03
#lgridz = lgrid/8. ###lgrid*4 # Four times longer than Lgrid in zdirection.


coords = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Coordinates",numThreads=1)
mass = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Mass",numThreads=1)


if(snip==1): # Have to back out hsmooth from density
    density = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/Density",numThreads=1)
    print "density= ", density
    hsmooth = (density/mass)**(-0.3333)*2.39  #2.39 conversion factor for 58 neighbors?  
    print "hsmooth= ",hsmooth
else:
    hsmooth = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/SmoothingLength",numThreads=1)

mass *= 1.e+10
print "mass= ", mass


#temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1,physicalUnits=0,noH=0)
#SFR = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/StarFormationRate",numThreads=1)

hydrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
#oxygen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
#silicon = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Silicon",numThreads=1)

NEQ = 0

if(NEQ):
    if(snip==1):
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
    else:
        if(snip==0):
            print "chem = 0 " 
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

print "done reading "

enviro = True
if(envirocat == "All"):
    enviro = False

M_200 = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/Enviro/M_200",numThreads=1,noH=0,physicalUnits=1)
fRvir = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/Enviro/fRvir_min",numThreads=1,noH=0,physicalUnits=1)

if(enviro):
    print "Enviro maxes/mins= ",np.max(M_200), np.min(M_200), np.max(fRvir), np.min(fRvir)
    fRvir_low = 0.0  

    if(envirocat == "M110"):
        fRvir_high = 1.0
        M_200_low = 10**10.5
        M_200_high = 10**11.5        
        enviro_label = "M105_M115"
    if(envirocat == "M120"):
        fRvir_high = 1.0
        M_200_low = 10**11.5
        M_200_high = 10**12.5        
        enviro_label = "M115_M125"
    if(envirocat == "M130"):
        fRvir_high = 1.0
        M_200_low = 10**12.5
        M_200_high = 10**15.5        
        enviro_label = "M125_M155"
    if(envirocat == "IGM"):
        fRvir_high = 1000.
        fRvir_low = 1.0  
        M_200_low = 10**10.5
        M_200_high = 10**15.5        
        enviro_label = "IGM"

    if(envirocat == "M110_2R"):
        fRvir_high = 2.0
        M_200_low = 10**10.5
        M_200_high = 10**11.5        
        enviro_label = "M105_M115_2R"
    if(envirocat == "M120_2R"):
        fRvir_high = 2.0
        M_200_low = 10**11.5
        M_200_high = 10**12.5        
        enviro_label = "M115_M125_2R"
    if(envirocat == "M130_2R"):
        fRvir_high = 2.0
        M_200_low = 10**12.5
        M_200_high = 10**15.5        
        enviro_label = "M125_M155_2R"
    if(envirocat == "IGM_2R"):
        fRvir_high = 1000.
        fRvir_low = 2.0  
        M_200_low = 10**10.5
        M_200_high = 10**15.5        
        enviro_label = "IGM_2R"


else:
    fRvir_high = 1e+06
    fRvir_low = 0
    M_200_low = 0
    M_200_high = 1e+20 
    enviro_label = "All"

if(ionshort=="h1"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'HI'
    colmin = 12
    colmax = 18
    coldens_limit = 11.0

if(ionshort=="si2"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'SiII'
    colmin = 11
    colmax = 15
    coldens_limit = 10.0

if(ionshort=="si3"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'SiIII'
    colmin = 11
    colmax = 15
    coldens_limit = 10.0

if(ionshort=="si4"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIV",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'SiIV'
    colmin = 11
    colmax = 15
    coldens_limit = 10.0

if(ionshort=="c2"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'CII'
    colmin = 11
    colmax = 15
    coldens_limit = 11.0

if(ionshort=="c4"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'CIV'
    colmin = 11
    colmax = 15
    coldens_limit = 10.0

if(ionshort=="mg2"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'MgII'
    colmin = 11
    colmax = 15
    coldens_limit = 11.0

if(ionshort=="o6"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'OVI'
    colmin = 11
    colmax = 15
    coldens_limit = 11.0

if(ionshort=="h1_o6"):
    ion_o6 = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1,noH=0,physicalUnits=1)
    ion_h1 = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1,noH=0,physicalUnits=1)
    ionmass_o6 = mass*hydrogen*ion_o6
    ionmass_h1 = mass*hydrogen*ion_h1
    ionname = 'HI_OVI'
    colmin = 11
    colmax = 15
    coldens_limit = 11.0

if(ionshort=="o7"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'OVII'
    colmin = 12
    colmax = 16
    coldens_limit = 13.0

if(ionshort=="o8"):
    ion = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1,noH=0,physicalUnits=1)
    ionmass = mass*hydrogen*ion
    ionname = 'OVIII'
    colmin = 12
    colmax = 16
    coldens_limit = 13.0

#indexes = np.where((ion>1e-08) & (fRvir>= fRvir_low) & (fRvir< fRvir_high) & (M_200>=M_200_low) & (M_200<M_200_high))
indexes = np.where((fRvir>= fRvir_low) & (fRvir< fRvir_high) & (M_200>=M_200_low) & (M_200<M_200_high))

for i in range(nslice):


    sliceinfostr = 'slice-%d'%i
    
    zcoord = boxsize/nslice*(i+0.5)
    center = np.array([boxsize/2.,boxsize/2.,zcoord])

    print "center= ", center

    result_h1 = coldens.main(coords[indexes], hsmooth[indexes], ionmass_h1[indexes], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.%s.l%3.1f.%s.png'%(snapname,ionname,sliceinfostr,lgrid,enviro_label),Vmin=colmin, Vmax=colmax,ion=ionname,redshift=redshift,extralabel=ionname,haloinfostr=sliceinfostr,docbar=docbar)  

    result_o6 = coldens.main(coords[indexes], hsmooth[indexes], ionmass_o6[indexes], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.%s.l%3.1f.%s.png'%(snapname,ionname,sliceinfostr,lgrid,enviro_label),Vmin=colmin, Vmax=colmax,ion=ionname,redshift=redshift,extralabel=ionname,haloinfostr=sliceinfostr,docbar=docbar)  

    
    if(do_metals):
        metallicity = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Metallicity",numThreads=1)

        result_Z = coldens.main(coords[indexes], hsmooth[indexes], ionmass_h1[indexes]*metallicity[indexes], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens_Z.%s.%s.%s.l%3.1f.%s.png'%(snapname,ionname,sliceinfostr,lgrid,enviro_label),Vmin=colmin, Vmax=colmax,ion=ionname,redshift=redshift,extralabel=ionname,haloinfostr=sliceinfostr,docbar=docbar)  

    sumslice = np.sum(result_h1)

    sumslice = np.sum(10**(result_h1[np.where(np.isfinite(result_h1))]))
    sumslice /= (ngrid*ngrid)
    if(do_metals):
        f = file('colion_map_Z.%s.%s.%s.l%3.1f.%s.dat'%(snapname,ionname,sliceinfostr,lgrid,enviro_label), 'w')
    else:
        f = file('colion_map.%s.%s.%s.l%3.1f.%s.dat'%(snapname,ionname,sliceinfostr,lgrid,enviro_label), 'w')
    for i in xrange(ngrid):
        for j in xrange(ngrid):
            if(result_h1[i,j] >= coldens_limit):
                if(do_metals):
                    f.write('%5d %5d % 5.2f %5.3e % 5.2f\n'%(i,j,result_h1[i,j],10**result_Z[i,j]/10**result_h1[i,j],result_o6[i,j]))
                else:
                    f.write('%3d %3d % 5.2f\n'%(i,j,result_h1[i,j]))
                    
    #f.write("#sumslice= %5.3e\n"%sumslice)
    print "sumslice = ", sumslice
    f.close()
                        
