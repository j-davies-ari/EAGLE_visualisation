#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
tag= sys.argv[1]  
#'047_z000p000.ioneq'
#tag='ioneq_025'
#center = np.array([6.98,5.21,6.55])
#center= np.array([17.5995,14.08347,15.8329])  #snapshot 31
#center= np.array([16.3672,13.0542,14.6979])  #z=0.271
#center = np.array([15.2534, 10.9404,  9.0412])
#center = np.array([15.2691, 10.934, 9.03164])
center = np.array([15.2974, 10.9540,  9.0412])

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

result = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.png'%tag,Vmin=16, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogenrot.png'%tag,Vmin=16, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*carbon/12.0, center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.carbonrot.png'%tag,Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.png'%tag,Vmin=12, Vmax=18)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygenrot.png'%tag,Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.png'%tag,Vmin=13, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1rot.png'%tag,Vmin=13, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4rot.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6rot.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7rot.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid, lgrid, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2rot.png'%tag,Vmin=11, Vmax=15,theta=90)





result = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.big.png'%tag,Vmin=16, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogenrot.big.png'%tag,Vmin=16, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*carbon/12.0, center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.carbonrot.big.png'%tag,Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.big.png'%tag,Vmin=12, Vmax=18)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygenrot.big.png'%tag,Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.big.png'%tag,Vmin=13, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1rot.big.png'%tag,Vmin=13, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.big.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.big.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.big.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.big.png'%tag,Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4rot.big.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6rot.big.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7rot.big.png'%tag,Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid2, lgrid2, lgridlong, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2rot.big.png'%tag,Vmin=11, Vmax=15,theta=90)

print "H all"
result = coldens.main(coords, hsmooth, mass*hydrogen, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.all.png'%tag,Vmin=16, Vmax=21)
print "O all"
result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.all.png'%tag,Vmin=12, Vmax=18)
print "OVI all"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.all.png'%tag,Vmin=11, Vmax=15)
print "OVII all"
result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, 12.5, 12.5, 12.5, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.all.png'%tag,Vmin=11, Vmax=15)
