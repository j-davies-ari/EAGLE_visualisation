#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,harr,oarr,o6arr,o7arr,o8arr,ngrid,lgrid,lin):
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
    o6hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o6arr)[0]
    o7hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o7arr)[0]
    o8hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o8arr)[0]

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
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(o7hist[i]/npart[i]),np.log10(o8hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],hhist[i]/npart[i],ohist[i]/npart[i],o6hist[i]/npart[i],o7hist[i]/npart[i],o8hist[i]/npart[i],npart[i]))
    f.close()

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
    mh = haloname.split("_")[2]
    ms = haloname.split("_")[3]
    sfr = haloname.split("_")[4]
    haloinfostr = 'lg M$_{200}=' + mh + '$, lg M$_{*}=' + ms + '$, SFR$=' + sfr + '$'
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

hydrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)
oxygen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)

if(snip=="1"):
    o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
    o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
    o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
else:
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

result_o6_z = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar)
result_o7_z = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII',redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar)
result_o8_z = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII',redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar)
result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=17, Vmax=21,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar)
result_o_z = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.z.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',redshift=redshift,extralabel='O',haloinfostr=haloinfostr,docbar=docbar)

coldens.make_colourmap('coldens.%s.fo6.z.l%3.1f.png'%(snapname,lgrid), result_o6_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVI}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo7.z.l%3.1f.png'%(snapname,lgrid), result_o7_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVII}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo8.z.l%3.1f.png'%(snapname,lgrid), result_o8_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVIII}$',haloinfostr=haloinfostr,docbar=docbar)

f = file('colion_map.%s.z.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_z[i,j],result_o_z[i,j],result_h_z[i,j],result_o_z[i,j],result_o6_z[i,j],result_o7_z[i,j],result_o8_z[i,j]))
f.close()

write_radial_column('columnave.%s.z.%3.1f.dat'%(snapname,lgrid),result_h_z,result_o_z,result_o6_z,result_o7_z,result_o8_z,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.z.%3.1f.dat'%(snapname,lgrid),10**result_h_z,10**result_o_z,10**result_o6_z,10**result_o7_z,10**result_o8_z,ngrid,lgrid,1)



result_o6_y = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar,theta=90)
result_o7_y = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII',redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar,theta=90)
result_o8_y = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII',redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar,theta=90)
result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=17, Vmax=21,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar,theta=90)
result_o_y = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.y.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',redshift=redshift,extralabel='O',haloinfostr=haloinfostr,docbar=docbar,theta=90)

coldens.make_colourmap('coldens.%s.fo6.y.l%3.1f.png'%(snapname,lgrid), result_o6_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVI}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo7.y.l%3.1f.png'%(snapname,lgrid), result_o7_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVII}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo8.y.l%3.1f.png'%(snapname,lgrid), result_o8_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVIII}$',haloinfostr=haloinfostr,docbar=docbar)

f = file('colion_map.%s.y.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_y[i,j],result_o_y[i,j],result_h_y[i,j],result_o_y[i,j],result_o6_y[i,j],result_o7_y[i,j],result_o8_y[i,j]))
f.close()

write_radial_column('columnave.%s.y.%3.1f.dat'%(snapname,lgrid),result_h_y,result_o_y,result_o6_y,result_o7_y,result_o8_y,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.y.%3.1f.dat'%(snapname,lgrid),10**result_h_y,10**result_o_y,10**result_o6_y,10**result_o7_y,10**result_o8_y,ngrid,lgrid,1)



result_o6_x = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',redshift=redshift,extralabel='OVI',haloinfostr=haloinfostr,docbar=docbar,phi=90)
result_o7_x = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVII',redshift=redshift,extralabel='OVII',haloinfostr=haloinfostr,docbar=docbar,phi=90)
result_o8_x = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.5, Vmax=16.5,ion='OVIII',redshift=redshift,extralabel='OVIII',haloinfostr=haloinfostr,docbar=docbar,phi=90)
result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=17, Vmax=21,ion='H',redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=docbar,phi=90)
result_o_x = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.x.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',redshift=redshift,extralabel='O',haloinfostr=haloinfostr,docbar=docbar,phi=90)

coldens.make_colourmap('coldens.%s.fo6.x.l%3.1f.png'%(snapname,lgrid), result_o6_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'fOVI', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVI}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo7.x.l%3.1f.png'%(snapname,lgrid), result_o7_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVII}$',haloinfostr=haloinfostr,docbar=docbar)
coldens.make_colourmap('coldens.%s.fo8.x.l%3.1f.png'%(snapname,lgrid), result_o8_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'fOVIII', ngrid, redshift=redshift, extralabel='f$_\mathrm{OVIII}$',haloinfostr=haloinfostr,docbar=docbar)

f = file('colion_map.%s.x.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_x[i,j],result_o_x[i,j],result_h_x[i,j],result_o_x[i,j],result_o6_x[i,j],result_o7_x[i,j],result_o8_x[i,j]))
f.close()

write_radial_column('columnave.%s.x.%3.1f.dat'%(snapname,lgrid),result_h_x,result_o_x,result_o6_x,result_o7_x,result_o8_x,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.x.%3.1f.dat'%(snapname,lgrid),10**result_h_x,10**result_o_x,10**result_o6_x,10**result_o7_x,10**result_o8_x,ngrid,lgrid,1)

