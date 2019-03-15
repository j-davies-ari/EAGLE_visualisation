#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,harr,oarr,h1arr,c2arr,c3arr,c4arr,o1arr,o6arr,mg2arr,si2arr,si3arr,si4arr,ngrid,lgrid,lin):
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
    c2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c2arr)[0]
    c3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c3arr)[0]
    c4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c4arr)[0]
    o1hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o1arr)[0]
    o6hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o6arr)[0]
    mg2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=mg2arr)[0]
    si2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si2arr)[0]
    si3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si3arr)[0]
    si4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si4arr)[0]

    f = file(filename, 'w')
    f.write('#kpclo kpchi Hcol Ocol  HIcol C4col O6col O7col O8col Mg2col n\n')
    for i in xrange(nbins):
        if lin > 0:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.3f %5.3f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(h1hist[i]/npart[i]),np.log10(c2hist[i]/npart[i]),np.log10(c3hist[i]/npart[i]),np.log10(c4hist[i]/npart[i]),np.log10(o1hist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(mg2hist[i]/npart[i]),np.log10(si2hist[i]/npart[i]),np.log10(si3hist[i]/npart[i]),np.log10(si4hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],h1hist[i]/npart[i],c2hist[i]/npart[i],c3hist[i]/npart[i],c4hist[i]/npart[i],o1hist[i]/npart[i],o6hist[i]/npart[i],mg2hist[i]/npart[i],si2hist[i]/npart[i],si3hist[i]/npart[i],si4hist[i]/npart[i],npart[i]))
    f.close()



def radial_bin_statistics(filename,harr,oarr,h1arr,c2arr,c3arr,c4arr,o1arr,o6arr,mg2arr,si2arr,si3arr,si4arr,ngrid,lgrid,lin):
#def write_radial_column(filename,harr,ngrid,lgrid):

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
    c2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c2arr)[0]
    c3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c3arr)[0]
    c4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=c4arr)[0]
    o1hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o1arr)[0]
    o6hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=o6arr)[0]
    mg2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=mg2arr)[0]
    si2hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si2arr)[0]
    si3hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si3arr)[0]
    si4hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si4arr)[0]

    f = file(filename, 'w')
    f.write('#kpclo kpchi Hcol Ocol  HIcol C4col O6col O7col O8col Mg2col n\n')
    for i in xrange(nbins):
        if lin > 0:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.3f %5.3f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(h1hist[i]/npart[i]),np.log10(c2hist[i]/npart[i]),np.log10(c3hist[i]/npart[i]),np.log10(c4hist[i]/npart[i]),np.log10(o1hist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(mg2hist[i]/npart[i]),np.log10(si2hist[i]/npart[i]),np.log10(si3hist[i]/npart[i]),np.log10(si4hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],h1hist[i]/npart[i],c2hist[i]/npart[i],c3hist[i]/npart[i],c4hist[i]/npart[i],o1hist[i]/npart[i],o6hist[i]/npart[i],mg2hist[i]/npart[i],si2hist[i]/npart[i],si3hist[i]/npart[i],si4hist[i]/npart[i],npart[i]))
    f.close()


#sim='/net/virgo/data5/oppenheimer/Halo_x001/data_001_x001/'
sim= '.'
input_filename_base = sys.argv[1]
snapname = sys.argv[2]
snip = sys.argv[3]
xcoord = float(sys.argv[4])
ycoord = float(sys.argv[5])
zcoord = float(sys.argv[6])
lgrid = float(sys.argv[7])
ngrid = int(sys.argv[8])
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
    h1 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
    print 'h1= ', h1
    c4 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
    print 'c4= ', c4
    o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
    o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
    o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
    mg2 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)
else:
    chem = eagle.readArray("SNAP", sim, input_filename_base, "/PartType0/ChemistryAbundances",physicalUnits=1,noH=0,numThreads=1)
    h1 = chem[:,1]
    c2 = chem[:,8]
    c3 = chem[:,9]
    c4 = chem[:,10]
    o1 = chem[:,23]
    o6 = chem[:,28]
#    o7 = chem[:,29]
#    o8 = chem[:,30]
    mg2 = chem[:,45] 
    si2 = chem[:,58]
    si3 = chem[:,59]
    si4 = chem[:,60]

result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H',npix=ngrid)
result_o_z = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O',npix=ngrid)
result_h1_z = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',npix=ngrid)
result_c2_z = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CII',npix=ngrid)
result_c3_z = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIII',npix=ngrid)
result_c4_z = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV',npix=ngrid)
result_o1_z = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o1.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OI',npix=ngrid)
result_o6_z = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',npix=ngrid)
result_mg2_z = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',npix=ngrid)
result_si2_z = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiII',npix=ngrid)
result_si3_z = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIII',npix=ngrid)
result_si4_z = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIV',npix=ngrid)

f = file('coldens_map.%s.z.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f\n'%(i,j,result_h_z[i,j],result_o_z[i,j],result_h1_z[i,j],result_c2_z[i,j],result_c3_z[i,j],result_c4_z[i,j],result_o1_z[i,j],result_o6_z[i,j],result_mg2_z[i,j],result_si2_z[i,j],result_si3_z[i,j],result_si4_z[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03))
f.close()

write_radial_column('columnave.%s.z.%3.1f.dat'%(snapname,lgrid),result_h_z,result_o_z,result_h1_z,result_c2_z,result_c3_z,result_c4_z,result_o1_z,result_o6_z,result_mg2_z,result_si2_z,result_si3_z,result_si4_z,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.z.%3.1f.dat'%(snapname,lgrid),10**result_h_z,10**result_o_z,10**result_h1_z,10**result_c2_z,10**result_c3_z,10**result_c4_z,10**result_o1_z,10**result_o6_z,10**result_mg2_z,10**result_si2_z,10**result_si3_z,10**result_si4_z,ngrid,lgrid,1)



result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H',theta=90,npix=ngrid)
result_o_y = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O',theta=90,npix=ngrid)
result_h1_y = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',theta=90,npix=ngrid)
result_c2_y = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CII',theta=90,npix=ngrid)
result_c3_y = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIII',theta=90,npix=ngrid)
result_c4_y = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV',theta=90,npix=ngrid)
result_o1_y = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o1.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OI',theta=90,npix=ngrid)
result_o6_y = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',theta=90,npix=ngrid)
result_mg2_y = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',theta=90,npix=ngrid)
result_si2_y = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiII',theta=90,npix=ngrid)
result_si3_y = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIII',theta=90,npix=ngrid)
result_si4_y = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIV',theta=90,npix=ngrid)

f = file('coldens_map.%s.y.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f\n'%(i,j,result_h_y[i,j],result_o_y[i,j],result_h1_y[i,j],result_c2_y[i,j],result_c3_y[i,j],result_c4_y[i,j],result_o1_y[i,j],result_o6_y[i,j],result_mg2_y[i,j],result_si2_y[i,j],result_si3_y[i,j],result_si4_y[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03))
f.close()

write_radial_column('columnave.%s.y.%3.1f.dat'%(snapname,lgrid),result_h_y,result_o_y,result_h1_y,result_c2_y,result_c3_y,result_c4_y,result_o1_y,result_o6_y,result_mg2_y,result_si2_y,result_si3_y,result_si4_y,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.y.%3.1f.dat'%(snapname,lgrid),10**result_h_y,10**result_o_y,10**result_h1_y,10**result_c2_y,10**result_c3_y,10**result_c4_y,10**result_o1_y,10**result_o6_y,10**result_mg2_y,10**result_si2_y,10**result_si3_y,10**result_si4_y,ngrid,lgrid,1)


result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=16, Vmax=21,ion='H',phi=90,npix=ngrid)
result_o_x = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.oxygen.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=18,ion='O',phi=90,npix=ngrid)
result_h1_x = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.h1.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='HI',phi=90,npix=ngrid)
result_c2_x = coldens.main(coords, hsmooth, mass*hydrogen*c2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CII',phi=90,npix=ngrid)
result_c3_x = coldens.main(coords, hsmooth, mass*hydrogen*c3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIII',phi=90,npix=ngrid)
result_c4_x = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.c4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='CIV',phi=90,npix=ngrid)
result_o1_x = coldens.main(coords, hsmooth, mass*hydrogen*o1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o1.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OI',phi=90,npix=ngrid)
result_o6_x = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.o6.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='OVI',phi=90,npix=ngrid)
result_mg2_x = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.mg2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgII',phi=90,npix=ngrid)
result_si2_x = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiII',phi=90,npix=ngrid)
result_si3_x = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIII',phi=90,npix=ngrid)
result_si4_x = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, 25,fig_name='coldens.%s.si4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiIV',phi=90,npix=ngrid)

f = file('coldens_map.%s.x.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f\n'%(i,j,result_h_x[i,j],result_o_x[i,j],result_h1_x[i,j],result_c2_x[i,j],result_c3_x[i,j],result_c4_x[i,j],result_o1_x[i,j],result_o6_x[i,j],result_mg2_x[i,j],result_si2_x[i,j],result_si3_x[i,j],result_si4_x[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03))
f.close()

write_radial_column('columnave.%s.x.%3.1f.dat'%(snapname,lgrid),result_h_x,result_o_x,result_h1_x,result_c2_x,result_c3_x,result_c4_x,result_o1_x,result_o6_x,result_mg2_x,result_si2_x,result_si3_x,result_si4_x,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.x.%3.1f.dat'%(snapname,lgrid),10**result_h_x,10**result_o_x,10**result_h1_x,10**result_c2_x,10**result_c3_x,10**result_c4_x,10**result_o1_x,10**result_o6_x,10**result_mg2_x,10**result_si2_x,10**result_si3_x,10**result_si4_x,ngrid,lgrid,1)

