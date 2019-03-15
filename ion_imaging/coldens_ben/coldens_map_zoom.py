#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,harr,oarr,h1arr,c4arr,o6arr,o7arr,o8arr,mg2arr,ne8arr,mg10arr,si12arr,n5arr,o1arr,o2arr,o3arr,o4arr,o5arr,c2arr,c3arr,si2arr,si3arr,si4arr,silowarr,ngrid,lgrid,lin):
#def write_radial_column(filename,harr,ngrid,lgrid):

    nbins= 20 #50
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
    ne8hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=ne8arr)[0]
    mg10hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=mg10arr)[0]
    si12hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si12arr)[0]
    n5hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=n5arr)[0]
    o1hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o1arr)[0]
    o2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o2arr)[0]
    o3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o3arr)[0]
    o4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o4arr)[0]
    o5hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o5arr)[0]
    c2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c2arr)[0]
    c3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c3arr)[0]
    si2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si2arr)[0]
    si3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si3arr)[0]
    si4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si4arr)[0]
    silowhist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=silowarr)[0]

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
    f.write('#kpclo kpchi Hcol Ocol  HIcol C4col O6col O7col O8col Mg2col Ne8col Mg10col Si12col N5col n O1col O2col O3col O4col o5col c2col c3col si2col si3col si4col silowcol n\n')
    for i in xrange(nbins):
        if lin > 0:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d\n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(h1hist[i]/npart[i]),np.log10(c4hist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(o7hist[i]/npart[i]),np.log10(o8hist[i]/npart[i]),np.log10(mg2hist[i]/npart[i]),np.log10(ne8hist[i]/npart[i]),np.log10(mg10hist[i]/npart[i]),np.log10(si12hist[i]/npart[i]),np.log10(n5hist[i]/npart[i]),npart[i],np.log10(o1hist[i]/npart[i]),np.log10(o2hist[i]/npart[i]),np.log10(o3hist[i]/npart[i]),np.log10(o4hist[i]/npart[i]),np.log10(o5hist[i]/npart[i]),np.log10(c2hist[i]/npart[i]),np.log10(c3hist[i]/npart[i]),np.log10(si2hist[i]/npart[i]),np.log10(si3hist[i]/npart[i]),np.log10(si4hist[i]/npart[i]),np.log10(silowhist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d\n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],h1hist[i]/npart[i],c4hist[i]/npart[i],o6hist[i]/npart[i],o7hist[i]/npart[i],o8hist[i]/npart[i],mg2hist[i]/npart[i],ne8hist[i]/npart[i],mg10hist[i]/npart[i],si12hist[i]/npart[i],n5hist[i]/npart[i],npart[i],o1hist[i]/npart[i],o2hist[i]/npart[i],o3hist[i]/npart[i],o4hist[i]/npart[i],o5hist[i]/npart[i],c2hist[i]/npart[i],c3hist[i]/npart[i],si2hist[i]/npart[i],si3hist[i]/npart[i],si4hist[i]/npart[i],silowhist[i]/npart[i],npart[i]))
    f.close()

#sim='/net/galaxy/data2/oppenheimer/noneqhalozoom_HM01/data/'
#sim='/net/galaxy/data2/oppenheimer/halozoomtest_janus/data/'
#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
#tag= sys.argv[1]  # Changed on 11/13/14.  
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
snip = int(sys.argv[3])
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
lgridz = lgrid ###lgrid*4 # Four times longer than Lgrid in zdirection.

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
center = center/hubble_param*aex
print "center= ", center
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

    c2 = np.where((np.isinf(c2)) | (np.isnan(c2)), 0.0, c2)
    c3 = np.where((np.isinf(c3)) | (np.isnan(c3)), 0.0, c3)
    c4 = np.where((np.isinf(c4)) | (np.isnan(c4)), 0.0, c4)
    o6 = np.where((np.isinf(o6)) | (np.isnan(o6)), 0.0, o6)
    o7 = np.where((np.isinf(o7)) | (np.isnan(o7)), 0.0, o7)
    o8 = np.where((np.isinf(o8)) | (np.isnan(o8)), 0.0, o8)
    si2 = np.where((np.isinf(si2)) | (np.isnan(si2)), 0.0, si2)
    si3 = np.where((np.isinf(si3)) | (np.isnan(si3)), 0.0, si3)
    si4 = np.where((np.isinf(si4)) | (np.isnan(si4)), 0.0, si4)


    #ne8 = o8*0.0
else:
    if(snip==0):
        print "chem = 0 " 
        chem = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
        ###chem = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChimesAbundances",physicalUnits=1,noH=0,numThreads=1)
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

print "snapname = ", snapname        
result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)
result_o_z = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.z.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',redshift=redshift,extralabel='O',haloinfostr=haloinfostr,docbar=docbar)
result_olow_z = coldens.main(coords[indexes_cool], hsmooth[indexes_cool], mass[indexes_cool]*oxygen[indexes_cool]/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.olow.z.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O-cool',redshift=redshift,extralabel='O-cool',haloinfostr=haloinfostr,docbar=docbar)
if(snip<=1):
    result_h1_z = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',redshift=redshift,extralabel='HI',haloinfostr=haloinfostr,docbar=docbar)
    result_c4_z = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIV',redshift=redshift,extralabel='CIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o6_z = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='OVI',redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar)
    result_o7_z = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.z.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVII',redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar)
    result_o8_z = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.z.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVIII',redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar)
    result_mg2_z = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',redshift=redshift,extralabel='MgII',haloinfostr=haloinfostr,docbar=docbar)
    result_ne8_z = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NeVIII',redshift=redshift,extralabel='NeVIII',haloinfostr=haloinfostr,docbar=docbar)
if(snip==0):
    result_mg10_z = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgX',redshift=redshift,extralabel='MgX',haloinfostr=haloinfostr,docbar=docbar)
    result_si12_z = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiXII',redshift=redshift,extralabel='SiXII',haloinfostr=haloinfostr,docbar=docbar)
    result_n5_z = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NV',redshift=redshift,extralabel='NV',haloinfostr=haloinfostr,docbar=docbar)
    result_o1_z = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o1.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OI',redshift=redshift,extralabel='OI',haloinfostr=haloinfostr,docbar=docbar)
    result_o2_z = coldens.main(coords, hsmooth, mass*hydrogen*o2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OII',redshift=redshift,extralabel='OII',haloinfostr=haloinfostr,docbar=docbar)
    result_o3_z = coldens.main(coords, hsmooth, mass*hydrogen*o3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIII',redshift=redshift,extralabel='OIII',haloinfostr=haloinfostr,docbar=docbar)
    result_o4_z = coldens.main(coords, hsmooth, mass*hydrogen*o4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIV',redshift=redshift,extralabel='OIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o5_z = coldens.main(coords, hsmooth, mass*hydrogen*o5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o5.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OV',redshift=redshift,extralabel='OV',haloinfostr=haloinfostr,docbar=docbar)
    result_c2_z = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CII',redshift=redshift,extralabel='CII',haloinfostr=haloinfostr,docbar=docbar)
    result_c3_z = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIII',redshift=redshift,extralabel='CIII',haloinfostr=haloinfostr,docbar=docbar)
    result_c6_z = coldens.main(coords, hsmooth, mass*hydrogen*c6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c6.z.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='CVI',redshift=redshift,extralabel='CVI',haloinfostr=haloinfostr,docbar=docbar)
    result_si2_z = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiII',redshift=redshift,extralabel='SiII',haloinfostr=haloinfostr,docbar=docbar)
    result_si3_z = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIII',redshift=redshift,extralabel='SiIII',haloinfostr=haloinfostr,docbar=docbar)
    result_si4_z = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIV',redshift=redshift,extralabel='SiIV',haloinfostr=haloinfostr,docbar=docbar)
    result_silow_z = coldens.main(coords[indexes_cool], hsmooth[indexes_cool], mass[indexes_cool]*silicon[indexes_cool]/28.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.silow.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='Si-cool',redshift=redshift,extralabel='Si-cool',haloinfostr=haloinfostr,docbar=docbar)
else:
    result_mg10_z = result_si12_z = result_n5_z = 0*result_h_z
    result_o1_z = result_o2_z = result_o3_z = result_o4_z =result_o5_z = 0*result_h_z

if(snip<=1):
    coldens.make_colourmap('coldens.%s.fo6.z.l%3.1f.png'%(snapname,lgrid), result_o6_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='OVI/O',haloinfostr=haloinfostr,docbar=docbar)
    coldens.make_colourmap('coldens.%s.fo7.z.l%3.1f.png'%(snapname,lgrid), result_o7_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='OVII/O',haloinfostr=haloinfostr,docbar=docbar)
    coldens.make_colourmap('coldens.%s.fo8.z.l%3.1f.png'%(snapname,lgrid), result_o8_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar)

    f = file('colion_map.%s.z.l%3.1f.dat'%(snapname,lgrid), 'w')
    for i in xrange(ngrid):
        for j in xrange(ngrid):
            f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.2f %5.2f\n'%(i,j,result_h_z[i,j],result_o_z[i,j],result_h1_z[i,j],result_c4_z[i,j],result_o6_z[i,j],result_o7_z[i,j],result_o8_z[i,j],result_mg2_z[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_z[i,j],result_mg10_z[i,j],result_si12_z[i,j],result_n5_z[i,j],result_c2_z[i,j],result_c3_z[i,j],result_si2_z[i,j],result_si3_z[i,j],result_si4_z[i,j],result_silow_z[i,j],result_olow_z[i,j]))
    f.close()

    write_radial_column('columnave.%s.z.%3.1f.dat'%(snapname,lgrid),result_h_z,result_o_z,result_h1_z,result_c4_z,result_o6_z,result_o7_z,result_o8_z,result_mg2_z,result_ne8_z,result_mg10_z,result_si12_z,result_n5_z,result_o1_z,result_o2_z,result_o3_z,result_o4_z,result_o5_z,result_c2_z,result_c3_z,result_si2_z,result_si3_z,result_si4_z,result_silow_z,ngrid,lgrid,0)
    write_radial_column('columnave_linsum.%s.z.%3.1f.dat'%(snapname,lgrid),10**result_h_z,10**result_o_z,10**result_h1_z,10**result_c4_z,10**result_o6_z,10**result_o7_z,10**result_o8_z,10**result_mg2_z,10**result_ne8_z,10**result_mg10_z,10**result_si12_z,10**result_n5_z,10**result_o1_z,10**result_o2_z,10**result_o3_z,10**result_o4_z,10**result_o5_z,10**result_c2_z,10**result_c3_z,10**result_si2_z,10**result_si3_z,10**result_si4_z,10**result_silow_z,ngrid,lgrid,1)

result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',theta=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)
result_o_y = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.y.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',theta=90,redshift=redshift,extralabel='O',haloinfostr=haloinfostr, docbar=docbar)
if(snip<=1):
    result_h1_y = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',theta=90,redshift=redshift,extralabel='HI',haloinfostr=haloinfostr,docbar=docbar)
    result_c4_y = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIV',theta=90,redshift=redshift,extralabel='CIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o6_y = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='OVI',theta=90,redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar)
    result_o7_y = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.y.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVII',theta=90,redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar)
    result_o8_y = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.y.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVIII',theta=90,redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar)
    result_mg2_y = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',theta=90,redshift=redshift,extralabel='MgII',haloinfostr=haloinfostr,docbar=docbar)
    result_ne8_y = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NeVIII',theta=90,redshift=redshift,extralabel='NeVIII',haloinfostr=haloinfostr,docbar=docbar)
if(snip==0):
    result_mg10_y = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgX',theta=90,redshift=redshift,extralabel='MgX',haloinfostr=haloinfostr,docbar=docbar)
    result_si12_y = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiXII',theta=90,redshift=redshift,extralabel='SiXII',haloinfostr=haloinfostr,docbar=docbar)
    result_n5_y = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NV',theta=90,redshift=redshift,extralabel='NV',haloinfostr=haloinfostr,docbar=docbar)
    result_o1_y = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o1.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OI',theta=90,redshift=redshift,extralabel='OI',haloinfostr=haloinfostr,docbar=docbar)
    result_o2_y = coldens.main(coords, hsmooth, mass*hydrogen*o2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OII',theta=90,redshift=redshift,extralabel='OII',haloinfostr=haloinfostr,docbar=docbar)
    result_o3_y = coldens.main(coords, hsmooth, mass*hydrogen*o3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIII',theta=90,redshift=redshift,extralabel='OIII',haloinfostr=haloinfostr,docbar=docbar)
    result_o4_y = coldens.main(coords, hsmooth, mass*hydrogen*o4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIV',theta=90,redshift=redshift,extralabel='OIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o5_y = coldens.main(coords, hsmooth, mass*hydrogen*o5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o5.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OV',theta=90,redshift=redshift,extralabel='OV',haloinfostr=haloinfostr,docbar=docbar)
    result_c2_y = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CII',theta=90,redshift=redshift,extralabel='CII',haloinfostr=haloinfostr,docbar=docbar)
    result_c3_y = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIII',theta=90,redshift=redshift,extralabel='CIII',haloinfostr=haloinfostr,docbar=docbar)
    result_c6_y = coldens.main(coords, hsmooth, mass*hydrogen*c6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c6.y.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='CVI',theta=90,redshift=redshift,extralabel='CVI',haloinfostr=haloinfostr,docbar=docbar)
    result_si2_y = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiII',theta=90,redshift=redshift,extralabel='SiII',haloinfostr=haloinfostr,docbar=docbar)
    result_si3_y = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIII',theta=90,redshift=redshift,extralabel='SiIII',haloinfostr=haloinfostr,docbar=docbar)
    result_si4_y = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIV',theta=90,redshift=redshift,extralabel='SiIV',haloinfostr=haloinfostr,docbar=docbar)
    result_silow_y = coldens.main(coords[indexes_cool], hsmooth[indexes_cool], mass[indexes_cool]*silicon[indexes_cool]/28.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.silow.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='Si-cool',theta=90,redshift=redshift,extralabel='Si-cool',haloinfostr=haloinfostr,docbar=docbar)
else:
    result_mg10_y = result_si12_y = result_n5_y = 0*result_h_y
    result_o1_y = result_o2_y = result_o3_y = result_o4_y =result_o5_y = 0*result_h_y

if(snip<=1):    
    coldens.make_colourmap('coldens.%s.fo6.y.l%3.1f.png'%(snapname,lgrid), result_o6_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='O/OVI',haloinfostr=haloinfostr,docbar=docbar)
    coldens.make_colourmap('coldens.%s.fo7.y.l%3.1f.png'%(snapname,lgrid), result_o7_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='O/OVII',haloinfostr=haloinfostr,docbar=docbar)
    coldens.make_colourmap('coldens.%s.fo8.y.l%3.1f.png'%(snapname,lgrid), result_o8_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar)

    f = file('colion_map.%s.y.l%3.1f.dat'%(snapname,lgrid), 'w')
    for i in xrange(ngrid):
        for j in xrange(ngrid):
            f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.2f\n'%(i,j,result_h_y[i,j],result_o_y[i,j],result_h1_y[i,j],result_c4_y[i,j],result_o6_y[i,j],result_o7_y[i,j],result_o8_y[i,j],result_mg2_y[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_y[i,j],result_mg10_y[i,j],result_si12_y[i,j],result_n5_y[i,j],result_c2_y[i,j],result_c3_y[i,j],result_si2_y[i,j],result_si3_y[i,j],result_si4_y[i,j],result_silow_y[i,j]))
    f.close()

    write_radial_column('columnave.%s.y.%3.1f.dat'%(snapname,lgrid),result_h_y,result_o_y,result_h1_y,result_c4_y,result_o6_y,result_o7_y,result_o8_y,result_mg2_y,result_ne8_y,result_mg10_y,result_si12_y,result_n5_y,result_o1_y,result_o2_y,result_o3_y,result_o4_y,result_o5_y,result_c2_y,result_c3_y,result_si2_y,result_si3_y,result_si4_y,result_silow_y,ngrid,lgrid,0)
    write_radial_column('columnave_linsum.%s.y.%3.1f.dat'%(snapname,lgrid),10**result_h_y,10**result_o_y,10**result_h1_y,10**result_c4_y,10**result_o6_y,10**result_o7_y,10**result_o8_y,10**result_mg2_y,10**result_ne8_y,10**result_mg10_y,10**result_si12_y,10**result_n5_y,10**result_o1_y,10**result_o2_y,10**result_o3_y,10**result_o4_y,10**result_o5_y,10**result_c2_y,10**result_c3_y,10**result_si2_y,10**result_si3_y,10**result_si4_y,10**result_silow_y,ngrid,lgrid,1)


result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',phi=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)
result_o_x = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.x.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',phi=90,redshift=redshift,extralabel='O',haloinfostr=haloinfostr,docbar=docbar)
if(snip<=1):
    result_h1_x = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',phi=90,redshift=redshift,extralabel='HI',haloinfostr=haloinfostr,docbar=docbar)
    result_c4_x = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIV',phi=90,redshift=redshift,extralabel='CIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o6_x = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='OVI',phi=90,redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar)
    result_o7_x = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.x.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVII',phi=90,redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar)
    result_o8_x = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.x.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='OVIII',phi=90,redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar)
    result_mg2_x = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',phi=90,redshift=redshift,extralabel='MgII',haloinfostr=haloinfostr,docbar=docbar)
    result_ne8_x = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NeVIII',phi=90,redshift=redshift,extralabel='NeVIII',haloinfostr=haloinfostr,docbar=docbar)
if(snip==0):
    result_mg10_x = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgX',phi=90,redshift=redshift,extralabel='MgX',haloinfostr=haloinfostr,docbar=docbar)
    result_si12_x = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiXII',phi=90,redshift=redshift,extralabel='SiXII',haloinfostr=haloinfostr,docbar=docbar)
    result_n5_x = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NV',phi=90,redshift=redshift,extralabel='NV',haloinfostr=haloinfostr,docbar=docbar)
    result_o1_x = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o1.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OI',phi=90,redshift=redshift,extralabel='OI',haloinfostr=haloinfostr,docbar=docbar)
    result_o2_x = coldens.main(coords, hsmooth, mass*hydrogen*o2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OII',phi=90,redshift=redshift,extralabel='OII',haloinfostr=haloinfostr,docbar=docbar)
    result_o3_x = coldens.main(coords, hsmooth, mass*hydrogen*o3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIII',phi=90,redshift=redshift,extralabel='OIII',haloinfostr=haloinfostr,docbar=docbar)
    result_o4_x = coldens.main(coords, hsmooth, mass*hydrogen*o4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OIV',phi=90,redshift=redshift,extralabel='OIV',haloinfostr=haloinfostr,docbar=docbar)
    result_o5_x = coldens.main(coords, hsmooth, mass*hydrogen*o5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o5.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.0, Vmax=16.0,ion='OV',phi=90,redshift=redshift,extralabel='OV',haloinfostr=haloinfostr,docbar=docbar)
    result_c2_x = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CII',phi=90,redshift=redshift,extralabel='CII',haloinfostr=haloinfostr,docbar=docbar)
    result_c3_x = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='CIII',phi=90,redshift=redshift,extralabel='CIII',haloinfostr=haloinfostr,docbar=docbar)
    result_c6_x = coldens.main(coords, hsmooth, mass*hydrogen*c6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c6.x.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='CVI',phi=90,redshift=redshift,extralabel='CVI',haloinfostr=haloinfostr,docbar=docbar)
    result_si2_x = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiII',phi=90,redshift=redshift,extralabel='SiII',haloinfostr=haloinfostr,docbar=docbar)
    result_si3_x = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIII',phi=90,redshift=redshift,extralabel='SiIII',haloinfostr=haloinfostr,docbar=docbar)
    result_si4_x = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='SiIV',phi=90,redshift=redshift,extralabel='SiIV',haloinfostr=haloinfostr,docbar=docbar)
    result_silow_x = coldens.main(coords[indexes_cool], hsmooth[indexes_cool], mass[indexes_cool]*silicon[indexes_cool]/28.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.silow.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='Si-cool',phi=90,redshift=redshift,extralabel='Si-cool',haloinfostr=haloinfostr,docbar=docbar)    
else:
    result_mg10_x = result_si12_x = result_n5_x = 0*result_h_x
    result_o1_x = result_o2_x = result_o3_x = result_o4_x =result_o5_x = 0*result_h_x

coldens.make_colourmap('coldens.%s.fo6.x.l%3.1f.png'%(snapname,lgrid), result_o6_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='OVI/O',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo7.x.l%3.1f.png'%(snapname,lgrid), result_o7_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='OVII/O',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo8.x.l%3.1f.png'%(snapname,lgrid), result_o8_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar)

if(snip<=1):
    f = file('colion_map.%s.x.l%3.1f.dat'%(snapname,lgrid), 'w')
    for i in xrange(ngrid):
        for j in xrange(ngrid):
            f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.2f\n'%(i,j,result_h_x[i,j],result_o_x[i,j],result_h1_x[i,j],result_c4_x[i,j],result_o6_x[i,j],result_o7_x[i,j],result_o8_x[i,j],result_mg2_x[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_x[i,j],result_mg10_x[i,j],result_si12_x[i,j],result_n5_x[i,j],result_c2_x[i,j],result_c3_x[i,j],result_si2_x[i,j],result_si3_x[i,j],result_si4_x[i,j],result_silow_x[i,j]))
    f.close()

    write_radial_column('columnave.%s.x.%3.1f.dat'%(snapname,lgrid),result_h_x,result_o_x,result_h1_x,result_c4_x,result_o6_x,result_o7_x,result_o8_x,result_mg2_x,result_ne8_x,result_mg10_x,result_si12_x,result_n5_x,result_o1_x,result_o2_x,result_o3_x,result_o4_x,result_o5_x,result_c2_x,result_c3_x,result_si2_x,result_si3_x,result_si4_x,result_silow_x,ngrid,lgrid,0)
    write_radial_column('columnave_linsum.%s.x.%3.1f.dat'%(snapname,lgrid),10**result_h_x,10**result_o_x,10**result_h1_x,10**result_c4_x,10**result_o6_x,10**result_o7_x,10**result_o8_x,10**result_mg2_x,10**result_ne8_x,10**result_mg10_x,10**result_si12_x,10**result_n5_x,10**result_o1_x,10**result_o2_x,10**result_o3_x,10**result_o4_x,10**result_o5_x,10**result_c2_x,10**result_c3_x,10**result_si2_x,10**result_si3_x,10**result_si4_x,10**result_silow_x,ngrid,lgrid,1)
