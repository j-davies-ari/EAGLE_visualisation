#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py
import turtle

def write_radial_column(filename,harr,oarr,h1arr,c4arr,o6arr,o7arr,o8arr,mg2arr,ne8arr,mg10arr,si12arr,n5arr,ngrid,lgrid,lin):
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
    ne8hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=ne8arr)[0]
    mg10hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=mg10arr)[0]
    si12hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=si12arr)[0]
    n5hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=n5arr)[0]

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
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],np.log10(hhist[i]/npart[i]),np.log10(ohist[i]/npart[i]),np.log10(h1hist[i]/npart[i]),np.log10(c4hist[i]/npart[i]),np.log10(o6hist[i]/npart[i]),np.log10(o7hist[i]/npart[i]),np.log10(o8hist[i]/npart[i]),np.log10(mg2hist[i]/npart[i]),np.log10(ne8hist[i]/npart[i]),np.log10(mg10hist[i]/npart[i]),np.log10(si12hist[i]/npart[i]),np.log10(n5hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5d \n'%(bins[i],bins[i+1],hhist[i]/npart[i],ohist[i]/npart[i],h1hist[i]/npart[i],c4hist[i]/npart[i],o6hist[i]/npart[i],o7hist[i]/npart[i],o8hist[i]/npart[i],mg2hist[i]/npart[i],ne8hist[i]/npart[i],mg10hist[i]/npart[i],si12hist[i]/npart[i],n5hist[i]/npart[i],npart[i]))
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
    mh = float(haloname.split("_")[2])
    ms = float(haloname.split("_")[3])
    sfr = float(haloname.split("_")[4])
    haloinfostr = 'lg M$_{200}=%4.1f$, lg M$_{*}=%3.1f$, SFR$=%4.2f$'%(mh,ms,sfr)
    mhalo = float(mh)
    #haloinfostr = 'lg M$_{200}=' + mh + '$, lg M$_{*}=' + ms + '$, SFR$=' + sfr + '$'
else:
    haloinfostr = runlabel
    haloinfostr = ''
    haloname = snapname.split("halo")[1]
    mhalo = float(haloname.split("_")[2])


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

omegaM = 0.307
omegaL = 1-omegaM
omegaratio = (omegaM+omegaL/(1+redshift)**3)

h = 0.6777
R200 = 1.63e-2*(10**float(mhalo)*h)**0.333/omegaratio**0.333/(1+redshift)/h
R200 *= 1e-03 # In Mpc.

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
    h1 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
    print 'h1= ', h1
    c4 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
    print 'c4= ', c4
    o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
    o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
    o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
    mg2 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)
    ne8 = mg10 = si12 = n5 = si2 = si3 = si4 = o8*0.0
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
    si2 = chem[:,58]
    si3 = chem[:,59]
    si4 = chem[:,60]


result_h_z = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.z.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',mhalo=mhalo,extralabel='H',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o_z = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.z.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',mhalo=mhalo,extralabel='O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_h1_z = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.z.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='H I',mhalo=mhalo,extralabel='H I',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_c4_z = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='C IV',mhalo=mhalo,extralabel='C IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o6_z = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.z.l%3.1f.png'%(snapname,lgrid),Vmin=12.7, Vmax=16.7,ion='O VI',mhalo=mhalo,extralabel='O VI',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o7_z = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.z.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VII',mhalo=mhalo,extralabel='O VII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o8_z = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.z.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VIII',mhalo=mhalo,extralabel='O VIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_mg2_z = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Mg II',mhalo=mhalo,extralabel='Mg II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si2_z = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si II',mhalo=mhalo,extralabel='Si II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si3_z = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si III',mhalo=mhalo,extralabel='Si III',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si4_z = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.z.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si IV',mhalo=mhalo,extralabel='Si IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
if(snip=="0"):
    result_ne8_z = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Ne VIII',mhalo=mhalo,extralabel='Ne VIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_mg10_z = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Mg X',mhalo=mhalo,extralabel='Mg X',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_si12_z = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Si XII',mhalo=mhalo,extralabel='Si XII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_n5_z = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.z.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='N V',mhalo=mhalo,extralabel='N V',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
else:
    result_ne8_z = result_mg10_z = result_si12_z = result_n5_z = 0*result_o6_z

coldens.make_colourmap('coldens.%s.fo6.z.l%3.1f.png'%(snapname,lgrid), result_o6_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'xOVI', ngrid, mhalo=mhalo, extralabel='OVI/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo7.z.l%3.1f.png'%(snapname,lgrid), result_o7_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVII', ngrid, mhalo=mhalo, extralabel='OVII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo8.z.l%3.1f.png'%(snapname,lgrid), result_o8_z-result_o_z, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVIII', ngrid, mhalo=mhalo, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)

f = file('colion_map.%s.z.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_z[i,j],result_o_z[i,j],result_h1_z[i,j],result_c4_z[i,j],result_o6_z[i,j],result_o7_z[i,j],result_o8_z[i,j],result_mg2_z[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_z[i,j],result_mg10_z[i,j],result_si12_z[i,j],result_n5_z[i,j]))
f.close()

write_radial_column('columnave.%s.z.%3.1f.dat'%(snapname,lgrid),result_h_z,result_o_z,result_h1_z,result_c4_z,result_o6_z,result_o7_z,result_o8_z,result_mg2_z,result_ne8_z,result_mg10_z,result_si12_z,result_n5_z,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.z.%3.1f.dat'%(snapname,lgrid),10**result_h_z,10**result_o_z,10**result_h1_z,10**result_c4_z,10**result_o6_z,10**result_o7_z,10**result_o8_z,10**result_mg2_z,10**result_ne8_z,10**result_mg10_z,10**result_si12_z,10**result_n5_z,ngrid,lgrid,1)

result_h_y = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.y.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',theta=90,mhalo=mhalo,extralabel='H',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o_y = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.y.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',theta=90,mhalo=mhalo,extralabel='O',haloinfostr=haloinfostr, docbar=docbar, R200=R200)
result_h1_y = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.y.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='H I',theta=90,mhalo=mhalo,extralabel='H I',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_c4_y = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='C IV',theta=90,mhalo=mhalo,extralabel='C IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o6_y = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.y.l%3.1f.png'%(snapname,lgrid),Vmin=12.7, Vmax=16.7,ion='O VI',theta=90,mhalo=mhalo,extralabel='O VI',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o7_y = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.y.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VII',theta=90,mhalo=mhalo,extralabel='O VII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o8_y = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.y.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VIII',theta=90,mhalo=mhalo,extralabel='O VIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_mg2_y = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Mg II',theta=90,mhalo=mhalo,extralabel='Mg II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si2_y = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si II',theta=90,mhalo=mhalo,extralabel='Si II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si3_y = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si III',theta=90,mhalo=mhalo,extralabel='Si III',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si4_y = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.y.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si IV',theta=90,mhalo=mhalo,extralabel='Si IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
if(snip=="0"):
    result_ne8_y = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Ne VIII',theta=90,mhalo=mhalo,extralabel='Ne VIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_mg10_y = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Mg X',theta=90,mhalo=mhalo,extralabel='Mg X',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_si12_y = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Si XII',theta=90,mhalo=mhalo,extralabel='Si XII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_n5_y = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.y.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='N V',theta=90,mhalo=mhalo,extralabel='N V',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
else:
    result_ne8_y = result_mg10_y = result_si12_y = result_n5_y = 0*result_o6_y

coldens.make_colourmap('coldens.%s.fo6.y.l%3.1f.png'%(snapname,lgrid), result_o6_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'xOVI', ngrid, mhalo=mhalo, extralabel='OVI/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo7.y.l%3.1f.png'%(snapname,lgrid), result_o7_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVII', ngrid, mhalo=mhalo, extralabel='OVII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo8.y.l%3.1f.png'%(snapname,lgrid), result_o8_y-result_o_y, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVIII', ngrid, mhalo=mhalo, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)

f = file('colion_map.%s.y.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_y[i,j],result_o_y[i,j],result_h1_y[i,j],result_c4_y[i,j],result_o6_y[i,j],result_o7_y[i,j],result_o8_y[i,j],result_mg2_y[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_y[i,j],result_mg10_y[i,j],result_si12_y[i,j],result_n5_y[i,j]))
f.close()

write_radial_column('columnave.%s.y.%3.1f.dat'%(snapname,lgrid),result_h_y,result_o_y,result_h1_y,result_c4_y,result_o6_y,result_o7_y,result_o8_y,result_mg2_y,result_ne8_y,result_mg10_y,result_si12_y,result_n5_y,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.y.%3.1f.dat'%(snapname,lgrid),10**result_h_y,10**result_o_y,10**result_h1_y,10**result_c4_y,10**result_o6_y,10**result_o7_y,10**result_o8_y,10**result_mg2_y,10**result_ne8_y,10**result_mg10_y,10**result_si12_y,10**result_n5_y,ngrid,lgrid,1)

result_h_x = coldens.main(coords, hsmooth, mass*hydrogen, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.hydrogen.x.l%3.1f.png'%(snapname,lgrid),Vmin=17.5, Vmax=21.5,ion='H',phi=90,mhalo=mhalo,extralabel='H',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o_x = coldens.main(coords, hsmooth, mass*oxygen/16.0, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.oxygen.x.l%3.1f.png'%(snapname,lgrid),Vmin=14, Vmax=18,ion='O',phi=90,mhalo=mhalo,extralabel='O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_h1_x = coldens.main(coords, hsmooth, mass*hydrogen*h1, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.h1.x.l%3.1f.png'%(snapname,lgrid),Vmin=13, Vmax=21,ion='H I',phi=90,mhalo=mhalo,extralabel='H I',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_c4_x = coldens.main(coords, hsmooth, mass*hydrogen*c4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.c4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='C IV',phi=90,mhalo=mhalo,extralabel='C IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o6_x = coldens.main(coords, hsmooth, mass*hydrogen*o6, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o6.x.l%3.1f.png'%(snapname,lgrid),Vmin=12.7, Vmax=16.7,ion='O VI',phi=90,mhalo=mhalo,extralabel='O VI',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o7_x = coldens.main(coords, hsmooth, mass*hydrogen*o7, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o7.x.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VII',phi=90,mhalo=mhalo,extralabel='O VII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_o8_x = coldens.main(coords, hsmooth, mass*hydrogen*o8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.o8.x.l%3.1f.png'%(snapname,lgrid),Vmin=13.0, Vmax=17.0,ion='O VIII',phi=90,mhalo=mhalo,extralabel='O VIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_mg2_x = coldens.main(coords, hsmooth, mass*hydrogen*mg2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='Mg II',phi=90,mhalo=mhalo,extralabel='Mg II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si2_y = coldens.main(coords, hsmooth, mass*hydrogen*si2, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si2.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si II',phi=90,mhalo=mhalo,extralabel='Si II',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si3_y = coldens.main(coords, hsmooth, mass*hydrogen*si3, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si3.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si III',phi=90,mhalo=mhalo,extralabel='Si III',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
result_si4_y = coldens.main(coords, hsmooth, mass*hydrogen*si4, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si4.x.l%3.1f.png'%(snapname,lgrid),Vmin=11.5, Vmax=15.5,ion='Si IV',phi=90,mhalo=mhalo,extralabel='Si IV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
if(snip=="0"):
    result_ne8_x = coldens.main(coords, hsmooth, mass*hydrogen*ne8, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.ne8.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NeVIII',phi=90,mhalo=mhalo,extralabel='NeVIII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_mg10_x = coldens.main(coords, hsmooth, mass*hydrogen*mg10, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.mg10.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='MgX',phi=90,mhalo=mhalo,extralabel='MgX',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_si12_x = coldens.main(coords, hsmooth, mass*hydrogen*si12, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.si12.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='SiXII',phi=90,mhalo=mhalo,extralabel='SiXII',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
    result_n5_x = coldens.main(coords, hsmooth, mass*hydrogen*n5, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.n5.x.l%3.1f.png'%(snapname,lgrid),Vmin=11, Vmax=15,ion='NV',phi=90,mhalo=mhalo,extralabel='NV',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
else:
    result_ne8_x = result_mg10_x = result_si12_x = result_n5_x = 0*result_o6_x

coldens.make_colourmap('coldens.%s.fo6.x.l%3.1f.png'%(snapname,lgrid), result_o6_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -4.0, -1.0, 'xOVI', ngrid, mhalo=mhalo, extralabel='OVI/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo7.x.l%3.1f.png'%(snapname,lgrid), result_o7_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVII', ngrid, mhalo=mhalo, extralabel='OVII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)
coldens.make_colourmap('coldens.%s.fo8.x.l%3.1f.png'%(snapname,lgrid), result_o8_x-result_o_x, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -3.0, 0.0, 'xOVIII', ngrid, mhalo=mhalo, extralabel='OVIII/O',haloinfostr=haloinfostr,docbar=docbar, R200=R200)

f = file('colion_map.%s.x.l%3.1f.dat'%(snapname,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f % 5.2f %5.1f % 5.2f % 5.2f % 5.2f % 5.2f\n'%(i,j,result_h_x[i,j],result_o_x[i,j],result_h1_x[i,j],result_c4_x[i,j],result_o6_x[i,j],result_o7_x[i,j],result_o8_x[i,j],result_mg2_x[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03,result_ne8_x[i,j],result_mg10_x[i,j],result_si12_x[i,j],result_n5_x[i,j]))
f.close()

write_radial_column('columnave.%s.x.%3.1f.dat'%(snapname,lgrid),result_h_x,result_o_x,result_h1_x,result_c4_x,result_o6_x,result_o7_x,result_o8_x,result_mg2_x,result_ne8_x,result_mg10_x,result_si12_x,result_n5_x,ngrid,lgrid,0)
write_radial_column('columnave_linsum.%s.x.%3.1f.dat'%(snapname,lgrid),10**result_h_x,10**result_o_x,10**result_h1_x,10**result_c4_x,10**result_o6_x,10**result_o7_x,10**result_o8_x,10**result_mg2_x,10**result_ne8_x,10**result_mg10_x,10**result_si12_x,10**result_n5_x,ngrid,lgrid,1)

