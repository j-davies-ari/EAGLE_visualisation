#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
tag='047_z000p000.ioneq'
#tag='ioneq_025'
#center = np.array([6.98,5.21,6.55])
#center= np.array([17.5995,14.08347,15.8329])  #snapshot 31
#center= np.array([16.3672,13.0542,14.6979])  #z=0.271
center = np.array([15.2466, 10.9404,  9.0412])
center = center/0.6777
ngrid = 400
lgrid = 600 
lgrid = lgrid/1e+03
lgrid2 = 4000 
lgrid2 = lgrid2/1e+03
lgridlong = 12500
lgridlong = lgridlong/1e+03

coords = eagle.readArray("SNAP", sim, tag, "/PartType0/Coordinates",numThreads=1)
hsmooth = eagle.readArray("SNAP", sim, tag, "/PartType0/SmoothingLength",numThreads=1)
mass = eagle.readArray("SNAP", sim, tag, "/PartType0/Mass",numThreads=1)
mass *= 1.e+10

density = eagle.readArray("SNAP", sim, tag, "/PartType0/Density",numThreads=1)

temperature = eagle.readArray("SNAP", sim, tag, "/PartType0/Temperature",numThreads=1)

metals = eagle.readArray("SNAP", sim, tag, "/PartType0/Metallicity",numThreads=1)

hydrogen = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)

carbon = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Carbon",numThreads=1)

oxygen = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Oxygen",numThreads=1)

chem = eagle.readArray("SNAP", sim, tag, "/PartType0/ChemistryAbundances",numThreads=1)

#chem = eagle.readArray("SNAP", sim, tag, "/PartType0/ChemistryAbundances"[:,1],numThreads=1)

#center = [8 8 8]

print center
print coords[0]

orig_stdout = sys.stdout
f = file('/net/galaxy/data2/oppenheimer/snap31.txt', 'w')
#f = file('/net/galaxy/data2/oppenheimer/snap25.txt', 'w')
sys.stdout = f

for i in range(0,len(mass)-1):
    print i,mass[i],density[i],temperature[i],hydrogen[i],oxygen[i],chem[i,0],chem[i,28]

sys.stdout = orig_stdout
f.close()

print "mass whole"
result = coldens.main(coords, hsmooth, mass, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.mass.whole.png'%tag,Vmin=16, Vmax=21)
print "Z whole"
result = coldens.main(coords, hsmooth, mass*metals, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.metals.whole.png'%tag,Vmin=16, Vmax=21)
print "H whole"
result = coldens.main(coords, hsmooth, mass*hydrogen, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.hydrogen.whole.png'%tag,Vmin=16, Vmax=21)
print "O whole"
result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.oxygen.whole.png'%tag,Vmin=12, Vmax=18)
print "OI whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,23], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o1.whole.png'%tag,Vmin=11, Vmax=15)
print "OII whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,24], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o2.whole.png'%tag,Vmin=11, Vmax=15)
print "OIII whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,25], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o3.whole.png'%tag,Vmin=11, Vmax=15)
print "OIV whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,26], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o4.whole.png'%tag,Vmin=11, Vmax=15)
print "OV whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,27], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o5.whole.png'%tag,Vmin=11, Vmax=15)
print "OVI whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o6.whole.png'%tag,Vmin=11, Vmax=15)
print "OVII whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o7.whole.png'%tag,Vmin=11, Vmax=15)
print "OVIII whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,30], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o8.whole.png'%tag,Vmin=11, Vmax=15)
print "OVIX whole"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,31], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='%s.o9.whole.png'%tag,Vmin=11, Vmax=15)
