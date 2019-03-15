#import numpy
import eagle
#import sys
#import os
import coldens
import numpy as np
sim='/net/galaxy/data2/oppenheimer/noneqhalozoomtest_snaprestart/data/'
tag='ioneq_026_z000p183'


coords = eagle.readArray("SNAP", sim, tag, "/PartType0/Coordinates",numThreads=1)
hsmooth = eagle.readArray("SNAP", sim, tag, "/PartType0/SmoothingLength",numThreads=1)
mass = eagle.readArray("SNAP", sim, tag, "/PartType0/Mass",numThreads=1)
mass *= 1.e+10

hydrogen = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)

carbon = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Carbon",numThreads=1)

oxygen = eagle.readArray("SNAP", sim, tag, "/PartType0/ElementAbundance/Oxygen",numThreads=1)

chem = eagle.readArray("SNAP", sim, tag, "/PartType0/ChemistryAbundances",numThreads=1)


#chem = eagle.readArray("SNAP", sim, tag, "/PartType0/ChemistryAbundances"[:,1],numThreads=1)

#center = [8 8 8]
center = np.array([14.26,10.65,13.48])
center = np.array([6.98,5.21,6.55])
center= np.array([17.5995,14.08347,15.8329]) 
print center
print coords[0]

result = coldens.main(coords, hsmooth, mass*hydrogen, center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_hydrogen_z02.pdf',Vmin=16, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen, center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_hydrogen_z02rot.pdf',Vmin=16, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*carbon/12.0, center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_carbon_z02rot.pdf',Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_oxygen_z02.pdf',Vmin=12, Vmax=18)

result = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_oxygen_z02rot.pdf',Vmin=12, Vmax=18,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_h1_z02.pdf',Vmin=13, Vmax=21)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_h1_z02rot.pdf',Vmin=13, Vmax=21,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_c4_z02.pdf',Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_o6_z02.pdf',Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_o7_z02.pdf',Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_mg2_z02.pdf',Vmin=11, Vmax=15)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_c4_z02rot.pdf',Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_o6_z02rot.pdf',Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_o7_z02rot.pdf',Vmin=11, Vmax=15,theta=90)

result = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, 0.5, 0.5, 0.5, 400, 400, 58, 25,fig_name='test_mg2_z02rot.pdf',Vmin=11, Vmax=15,theta=90)


