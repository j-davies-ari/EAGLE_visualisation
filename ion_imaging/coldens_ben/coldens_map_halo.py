#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,harr,oarr,h1arr,c4arr,o6arr,o7arr,o8arr,mg2arr,ngrid,lgrid,lin):
#def write_radial_column(filename,harr,ngrid,lgrid):

    nbins= 20
    dist = harr*0.0

    for i in range(0,ngrid-1):
        for j in range(0,ngrid-1):
            xcoord = abs(ngrid/2.-0.5-i)/ngrid*lgrid*1e+03
            ycoord = abs(ngrid/2.-0.5-j)/ngrid*lgrid*1e+03
            dist[i,j] = np.sqrt(xcoord**2+ycoord**2)
            #if(dist[i,j]<15.0): print dist[i,j],i,j,xcoord,ycoord

    npart,bins = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.])
    print npart
    print bins

    hhist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=harr)[0]
    ohist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=oarr)[0]
    h1hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=h1arr)[0]
    c4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c4arr)[0]
    o6hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o6arr)[0]
    o7hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o7arr)[0]
    o8hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o8arr)[0]
    mg2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=mg2arr)[0]

#    hhist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**((harr)))[0]
#    ohist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**((oarr)))[0]
#    h1hist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**(h1arr))[0]
#    c4hist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**(c4arr))[0]
#    o6hist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**(o6arr))[0]
#    o7hist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**(o7arr))[0]
#    mg2hist_lin = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=10**(mg2arr))[0]


#    print 'hhist_lin=', hhist_lin
#    print 'h1hist_lin=', h1hist_lin 

    f = file(filename, 'w')
    f.write('#kpclo kpchi Hcol Ocol  HIcol C4col O6col O7col O8col Mg2col n\n')
    for i in xrange(nbins):
        if lin > 0:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(h1hist[i]/npart[i]),np.log10(c4hist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(o7hist[i]/npart[i]),np.log10(o8hist[i]/npart[i]),np.log10(mg2hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],h1hist[i]/npart[i],c4hist[i]/npart[i],o6hist[i]/npart[i],o7hist[i]/npart[i],o8hist[i]/npart[i],mg2hist[i]/npart[i],npart[i]))
    f.close()

#sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
#tag= sys.argv[1]  # Changed on 11/13/14.  
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
xcoord = float(sys.argv[3])
ycoord = float(sys.argv[4])
zcoord = float(sys.argv[5])
lgrid = float(sys.argv[6])
ngrid = int(sys.argv[7])
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
lgridz = lgrid*4 # Four times longer than Lgrid in zdirection.

path = "snapshot_%s/snap_%s.0.hdf5"%(input_filename_base,input_filename_base)
print path

#filein = h5py.File(path)
#redshift = filein['Header'].attrs['Redshift']
#aex = 1/(1+redshift)
#center = center*aex

aex = eagle.readAttribute("SNAP", sim, input_filename_base, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute("SNAP", sim, input_filename_base, "/Header/HubbleParam")

center = center/hubble_param*aex

coords = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/Coordinates",numThreads=1)
hsmooth = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/SmoothingLength",numThreads=1)
mass = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/Mass",numThreads=1)
mass *= 1.e+10



hydrogen = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
#carbon = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ElementAbundance/Carbon",numThreads=1)
oxygen = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
chem = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)


#orig_stdout = sys.stdout
#f = file('/net/galaxy/data2/oppenheimer/snap31.txt', 'w')
#f = file('/net/galaxy/data2/oppenheimer/snap25.txt', 'w')
#sys.stdout = f

#for i in range(0,len(mass)-1):
#print i,mass[i],density[i],temperature[i],hydrogen[i],oxygen[i],chem[i,0],chem[i,28]

#sys.stdout = orig_stdout
#f.close()

result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H')
result_o_z = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O')
result_h1_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI')
result_c4_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV')
result_o6_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI')
result_o7_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII')
result_o8_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,30], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o8.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII')
result_mg2_z = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII')

write_radial_column('columnave.%s.z.%3.1f.dat'%(snapname,lgrid),result_h_z,result_o_z,result_h1_z,result_c4_z,result_o6_z,result_o7_z,result_o8_z,result_mg2_z,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.z.%3.1f.dat'%(snapname,lgrid),10**result_h_z,10**result_o_z,10**result_h1_z,10**result_c4_z,10**result_o6_z,10**result_o7_z,10**result_o8_z,10**result_mg2_z,ngrid,lgrid,1)

result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H',theta=90)
result_o_y = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O',theta=90)
result_h1_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',theta=90)
result_c4_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV',theta=90)
result_o6_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',theta=90)
result_o7_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII',theta=90)
result_o8_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,30], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o8.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII',theta=90)
result_mg2_y = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',theta=90)

write_radial_column('columnave.%s.y.%3.1f.dat'%(snapname,lgrid),result_h_y,result_o_y,result_h1_y,result_c4_y,result_o6_y,result_o7_y,result_o8_y,result_mg2_y,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.y.%3.1f.dat'%(snapname,lgrid),10**result_h_y,10**result_o_y,10**result_h1_y,10**result_c4_y,10**result_o6_y,10**result_o7_y,10**result_o8_y,10**result_mg2_y,ngrid,lgrid,1)

result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H',phi=90)
result_o_x = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O',phi=90)
result_h1_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,1], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',phi=90)
result_c4_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,10], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV',phi=90)
result_o6_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,28], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',phi=90)
result_o7_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,29], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o7.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII',phi=90)
result_o8_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,30], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o8.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII',phi=90)
result_mg2_x = coldens.main(coords, hsmooth, mass*hydrogen*chem[:,45], center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',phi=90)

write_radial_column('columnave.%s.x.%3.1f.dat'%(snapname,lgrid),result_h_x,result_o_x,result_h1_x,result_c4_x,result_o6_x,result_o7_x,result_o8_x,result_mg2_x,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.x.%3.1f.dat'%(snapname,lgrid),10**result_h_x,10**result_o_x,10**result_h1_x,10**result_c4_x,10**result_o6_x,10**result_o7_x,10**result_o8_x,10**result_mg2_x,ngrid,lgrid,1)
