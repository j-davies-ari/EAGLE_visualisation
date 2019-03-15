#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,harr,oarr,h1arr,c4arr,o6arr,o7arr,o8arr,mg2arr,ne8arr,mg10arr,si12arr,n5arr,o1arr,o2arr,o3arr,o4arr,o5arr,c2arr,c3arr,si2arr,si3arr,si4arr,silowarr,ngrid,lgrid,lin):
#def write_radial_column(filename,harr,ngrid,lgrid):

    nbins= 50
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

coords_orig = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Coordinates_orig",numThreads=1)


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
        
result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)

result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',theta=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)

result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',phi=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)


result_h_z = coldens.main(coords_orig, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.z_orig.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)

result_h_y = coldens.main(coords_orig, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.y_orig.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',theta=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)

result_h_x = coldens.main(coords_orig, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.x_orig.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',phi=90,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)

