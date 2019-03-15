#import numpy
import eagle
import sys
#import os
import coldens
import numpy as np
import h5py

def write_radial_column(filename,arr,ngrid,lgrid,lin):
#def write_radial_column(filename,harr,ngrid,lgrid):

    nbins= 20
    dist = arr*0.0

    for i in range(0,ngrid-1):
        for j in range(0,ngrid-1):
            xcoord = abs(ngrid/2.-0.5-i)/ngrid*lgrid*1e+03
            ycoord = abs(ngrid/2.-0.5-j)/ngrid*lgrid*1e+03
            dist[i,j] = np.sqrt(xcoord**2+ycoord**2)
            #if(dist[i,j]<15.0): print dist[i,j],i,j,xcoord,ycoord

    npart,bins = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.])
    print npart
    print bins

    hist = np.histogram(dist,nbins,range=[1e-05,lgrid*1e+03/2.],weights=arr)[0]

    f = file(filename, 'w')
    f.write('#kpclo kpchi col n\n')
    for i in xrange(nbins):
        if lin > 0:
            f.write('%5.1f %5.1f %5.2f %5d\n'%(bins[i],bins[i+1],np.log10(hist[i]/npart[i]),npart[i]))
        else:
            f.write('%5.1f %5.1f %5.2f %5d\n'%(bins[i],bins[i+1],hist[i]/npart[i],npart[i]))
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
ion = sys.argv[9]
direction = sys.argv[10]
modelname = sys.argv[11]
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
lgridz = lgrid*2 #Now 2 times as of 12/30/14 # Four times longer than Lgrid in zdirection.

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

redshift = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/Redshift")
aex = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/ExpansionFactor")
hubble_param = eagle.readAttribute(sniptag, sim, input_filename_base, "/Header/HubbleParam")
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
    hsmooth = (density/mass)**(-0.3333)*2.39  #2.39 conversion factor for 58 nei ghbors?  
    print "hsmooth= ",hsmooth
else:
    hsmooth = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/SmoothingLength",numThreads=1)

mass *= 1.e+10
print "mass= ", mass

hydrogen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Hydrogen",numThreads=1)

#h1 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
#c4 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
#o6 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
#o7 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
#o8 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
#mg2 = eagle.readArray("SNIP", sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)

f = file('colion_map.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid), 'w')
f.close()

if(ion=="hydrogen"):
    indexes = np.where(mass>0.0)
    denominator = mass[indexes]*hydrogen[indexes]
    ionname = 'H'
    colmin = 18
    colmax = 21
if(ion=="oxygen"):
    oxygen = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ElementAbundance/Oxygen",numThreads=1)
    indexes = np.where(oxygen>1e-08)
    denominator = mass[indexes]*oxygen[indexes]
    ionname = 'O'
    colmin = 14
    colmax = 18
if(ion=="h1"):
    h1 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/HydrogenI",numThreads=1)
    indexes = np.where(h1>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*h1[indexes]
    ionname = 'HI'
    colmin = 13
    colmax = 20
if(ion=="c2"):
    c2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonII",numThreads=1)
    indexes = np.where(c2>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*c2[indexes]
    ionname = 'CII'
    colmin = 11
    colmax = 15
if(ion=="c4"):
    c4 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/CarbonIV",numThreads=1)
    indexes = np.where(c4>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*c4[indexes]
    ionname = 'CIV'
    colmin = 12
    colmax = 15.5
if(ion=="n5"):
    n5 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/NitrogenV",numThreads=1)
    indexes = np.where(n5>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*n5[indexes]
    ionname = 'NV'
    colmin = 11
    colmax = 15
if(ion=="o6" or ion=="To6"):
    o6 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVI",numThreads=1)
    indexes = np.where(o6>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*o6[indexes]
    ionname = 'OVI'
    colmin = 11.5
    colmax = 15.5
    if(ion=="To6"):
        temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1)
        numerator = temperature[indexes]*mass[indexes]*hydrogen[indexes]*o6[indexes]
        colmin = 4.0
        colmax = 7.0
        ionname = 'T(OVI)'
if(ion=="o7"):
    o7 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVII",numThreads=1)
    indexes = np.where(o7>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*o7[indexes]
    ionname = 'OVII'
    colmin = 12.5
    colmax = 16.5
if(ion=="o8"):
    o8 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/OxygenVIII",numThreads=1)
    indexes = np.where(o8>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*o8[indexes]
    ionname = 'OVIII'
    colmin = 12.5
    colmax = 16.5
if(ion=="mg2"):
    mg2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/MagnesiumII",numThreads=1)
    indexes = np.where(mg2>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*mg2[indexes]
    ionname = 'MgII'
    colmin = 11
    colmax = 14.5
if(ion=="T"):
    indexes = np.where(mass>0.0)
    temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1)
    numerator = mass[indexes]*temperature[indexes]
    denominator = mass[indexes]
    ionname = 'T'
    colmin = 4.0
    colmax = 7.0
if(ion=="si2"):
    si2 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconII",numThreads=1)
    indexes = np.where(si2>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*si2[indexes]
    ionname = 'SiII'
    colmin = 11
    colmax = 15
if(ion=="si3"):
    si3 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIII",numThreads=1)
    indexes = np.where(si3>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*si3[indexes]
    ionname = 'SiIII'
    colmin = 11
    colmax = 15
if(ion=="si4"):
    si4 = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/ChemicalAbundances/SiliconIV",numThreads=1)
    indexes = np.where(si4>1e-08)
    denominator = mass[indexes]*hydrogen[indexes]*si4[indexes]
    ionname = 'SiIV'
    colmin = 11
    colmax = 15
if(ion=="temperature"):
    indexes = np.where(mass>0.0)
    temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1)
    entropy = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Entropy",numThreads=1)
    internalenergy = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/InternalEnergy",numThreads=1)
    nh = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Density",numThreads=1)
    numerator = mass[indexes]*temperature[indexes]/nh[indexes]**(2/3.)*internalenergy[indexes]/(entropy[indexes]*3/2.)
    denominator = mass[indexes]
    ionname = 'T'
    colmin = 4.0
    colmax = 7.0
if(ion=="nH"):
    indexes = np.where(mass>0.0)
    unit_mass_in_cgs = 1.989e33 * 1.0e10 
    unit_length_in_cgs = 3.0857e24 
    proton_mass_cgs = 1.67e-24
    unitDensity_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3
    density = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Density",numThreads=1)*unitDensity_in_cgs/proton_mass_cgs
    numerator = mass[indexes]*density[indexes]*hydrogen[indexes]
    denominator = mass[indexes]
    ionname = 'n$_{\mathrm{H}}$'
    colmin = -6.0
    colmax = -1.0
if(ion=="entropy"):
    indexes = np.where(mass>0.0)
    unit_mass_in_cgs = 1.989e33 * 1.0e10 
    unit_length_in_cgs = 3.0857e24 
    proton_mass_cgs = 1.67e-24
    unitDensity_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3
    density = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Density",numThreads=1)*unitDensity_in_cgs/proton_mass_cgs
    temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1)/1.16e+07 # keV 
    numerator = mass[indexes]*temperature[indexes]/(density[indexes]*hydrogen[indexes]*1.2)**(2/3.)
    denominator = mass[indexes]
    ionname = 'S'
    colmin = 1.0
    colmax = 3.0
if(ion=="Z"):
    metallicity = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Metallicity",numThreads=1)
    indexes = np.where(metallicity>2e-08)
    numerator = mass[indexes]*metallicity[indexes]
    denominator = mass[indexes]
    ionname = 'Z'
    colmin = -4.0
    colmax = -1.5
if(ion=="sfr"):
    unit_mass_in_cgs = 1.989e33 * 1.0e10 
    unit_length_in_cgs = 3.0857e24 
    proton_mass_cgs = 1.67e-24
    unitVelocity_in_cm_per_s = 1e+05

    unitTime_in_s =  unit_length_in_cgs/ unitVelocity_in_cm_per_s

    unitDensity_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3 #/proton_mass_cgs #*hubble_param**2*aex**-3

    unitPressure_in_cgs = unit_mass_in_cgs / unit_length_in_cgs / unitTime_in_s**2

    SF_SchmidtLawCoeff_MSUNpYRpKPC2       =  1.515e-4   # % Schaye & Della Vecchia 2008 (For Chabrier IMF)
    SF_SchmidtLawExponent                 =  1.4        
    SF_SchmidtLawHighDensExponent         =  2.0
    SF_SchmidtLawHighDensThresh_HpCM3     =  1000.0
    CM_PER_MPC = 3.086e+24
    SOLAR_MASS = 1.989e+33
    SEC_PER_YEAR = 3.155e+07
    GAMMA = 5.0/3.0
    GRAVITY = 6.672e-8
    PROTONMASS = 1.6726e-24
        
    lgridz = lgrid # we are making this same size in z direction
    temperature = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Temperature",numThreads=1)
    entropy = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Entropy",numThreads=1)
    internalenergy = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/InternalEnergy",numThreads=1)
    density = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Density",numThreads=1)
    nh = density*unitDensity_in_cgs/proton_mass_cgs
    metallicity = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/Metallicity",numThreads=1)

    KS_norm_cgs   = SF_SchmidtLawCoeff_MSUNpYRpKPC2 * ((CM_PER_MPC / 1.e3)**2 / SOLAR_MASS)**(SF_SchmidtLawExponent - 1)/ (1.e6 * SEC_PER_YEAR); #    ks_norm_cgs = KS_norm_cgs;
    KS_exponent = SF_SchmidtLawExponent;
    KSHighDens_thres_cgs = SF_SchmidtLawHighDensThresh_HpCM3 * PROTONMASS / 0.752 #check this- should be XH

    print 'unitPressure_in_cgs, unitDensity_cgs = ',unitPressure_in_cgs, unitDensity_in_cgs
    #press = entropy * pow(SphP[i].cky.WeightedDensity, GAMMA);
    density = density*unitDensity_in_cgs
    entropy = entropy*unitPressure_in_cgs*unitDensity_in_cgs**-GAMMA
    pressure_cgs = entropy * (density)**GAMMA
    #pressure_cgs = temperature * nh

    #numerator = mass[indexes]*temperature[indexes]/nh[indexes]**(2/3.)*internalenergy[indexes]/(entropy[indexes]*3/2.)

    #if(density  > KSHighDens_thres_cgs): #Need to fix this to generalize it.
    #    KS_norm_HighDens_cgs = KS_norm_cgs /(GAMMA * pressure_cgs / GRAVITY)**((SF_SchmidtLawHighDensExponent-SF_SchmidtLawExponent)/2)
    #    KS_norm_cgs = KS_norm_HighDens_cgs + mass*0.0;
    #    KS_exponent = SF_SchmidtLawHighDensExponent + mass*0.0;

    inv_tsfr_cgs = mass*0.0
    inv_tsfr_cgs = KS_norm_cgs * (pressure_cgs * GAMMA / GRAVITY)**((KS_exponent - 1) / 2.) #/* inverse star formation time scale in seconds */

    inv_tsfr_cgs = 0.606*5.99e-10*(nh*temperature/1e+03)**0.2
    sfr = mass*inv_tsfr_cgs #sm_dot       = P[i].Mass * inv_tsfr_cgs;

    #indexes = np.where((nh>1e-02) & (temperature<1e+06))
    #indexes = np.where((nh>1e-01) & (temperature<1e+06))

    indexes = np.where((nh>1e-01*(metallicity/0.002)**-0.64) & (temperature<1e+06))

    print 'nh_cgs= ', nh[indexes]
    #print 'density= ', density[indexes]
    print 'temperature= ', temperature[indexes]
    print 'entropy= ', entropy[indexes]
    print 'mass= ', mass[indexes]
    print 'metallicity= ', metallicity[indexes]

    print 'pressure_cgs= ', pressure_cgs[indexes]
    print 'inv_tsfr_cgs= ', inv_tsfr_cgs[indexes]
    print 'sfr= ', sfr[indexes]
    print 'maxsfr= ', max(sfr)
    print 'maxinvtimescale= ', max(inv_tsfr_cgs)

    #sfr = eagle.readArray(sniptag, sim, input_filename_base, "/PartType0/StarFormationRate",numThreads=1)



    denominator = sfr[indexes]
    print 'SFR= ',np.sum(sfr[indexes])
    #print 'inv_tsfr_cgs=' inv_tsfr_cgs[indexes]

    ionname = 'SFR'
    colmin = 10
    colmax = 15
    f = file('colion_map.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid), 'w')
    f.write('#TOTSFR %7.5f %7.3f\n'%(redshift, np.sum(sfr[indexes])))
    f.close()


if(direction=="z"):
    thetaangle = 0
    phiangle = 0
if(direction=="y"):
    thetaangle = 90
    phiangle = 0
if(direction=="x"):
    thetaangle = 0
    phiangle = 90


result_den = coldens.main(coords[indexes], hsmooth[indexes], denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid),Vmin=colmin, Vmax=colmax,ion=ion,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)

if(ion=="temperature" or ion=="T" or ion=="nH" or ion=="Z" or ion=="To6" or ion=="entropy"):
    result_num = coldens.main(coords[indexes], hsmooth[indexes], numerator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='num.png',Vmin=colmin, Vmax=colmax,ion=ion,theta=thetaangle,phi=phiangle,redshift=redshift,extralabel=ionname)
#make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix, redshift)
    ###coldens.make_colourmap('coldens.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, colmin, colmax, ion, ngrid, redshift=redshift, extralabel=ionname)
    coldens.make_colourmap('coldens.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, colmin, colmax, ion, ngrid, topbottomlabel=modelname, extralabel=ionname)
    print "result_num= ", result_num
    print "result_den= ", result_den

N10lgrid = 0.0
nN10lgrid = 0
N05lgrid = 0.0
nN05lgrid = 0
N025lgrid = 0.0
nN025lgrid = 0
N01lgrid = 0.0
nN01lgrid = 0
N033lgrid = 0.0
nN033lgrid = 0
N016lgrid = 0.0
nN016lgrid = 0
Nlog10lgrid = 0.0
Nlog05lgrid = 0.0
Nlog025lgrid = 0.0
Nlog01lgrid = 0.0
Nlog033lgrid = 0.0
Nlog016lgrid = 0.0

f = file('colion_map.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid), 'w')
for i in xrange(ngrid):
    for j in xrange(ngrid):
        f.write('%3d %3d % 5.2f %5.1f\n'%(i,j,result_den[i,j],np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03))
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/2.):
            N10lgrid += 10**(result_den[i,j]+1e-02)
            Nlog10lgrid += result_den[i,j]+1e-02
            nN10lgrid += 1
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/4.):
            N05lgrid += 10**(result_den[i,j]+1e-02)
            Nlog05lgrid += result_den[i,j]+1e-02
            nN05lgrid += 1
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/8.):
            N025lgrid += 10**(result_den[i,j]+1e-02)
            Nlog025lgrid += result_den[i,j]+1e-02
            nN025lgrid += 1
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/20.):
            N01lgrid += 10**(result_den[i,j]+1e-02)
            Nlog01lgrid += result_den[i,j]+1e-02
            nN01lgrid += 1
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/12.):
            N016lgrid += 10**(result_den[i,j]+1e-02)
            Nlog016lgrid += result_den[i,j]+1e-02
            nN016lgrid += 1
        if(np.sqrt((ngrid/2.-0.5-j)**2+(ngrid/2.-0.5-i)**2)/ngrid*lgrid*1e+03<lgrid*1e+03/6.):
            N033lgrid += 10**(result_den[i,j]+1e-02)
            Nlog033lgrid += result_den[i,j]+1e-02   
            nN033lgrid += 1
            
f.write('#TOTCOL_%s_%s %7.5f 1.0 %5.3f %6d 0.5 %5.3f %6d 0.25 %5.3f %6d 0.10 %5.3f %6d 0.33 %5.3f %6d 0.16 %5.3f %6d\n'%(ion,direction, redshift, np.log10(N10lgrid/nN10lgrid), nN10lgrid, np.log10(N05lgrid/nN05lgrid), nN05lgrid, np.log10(N025lgrid/nN025lgrid), nN025lgrid, np.log10(N01lgrid/nN01lgrid), nN01lgrid, np.log10(N033lgrid/nN033lgrid), nN033lgrid, np.log10(N016lgrid/nN016lgrid), nN016lgrid))
f.write('#TOTCLOG_%s_%s %7.5f 1.0 %5.3f %6d 0.5 %5.3f %6d 0.25 %5.3f %6d 0.10 %5.3f %6d 0.33 %5.3f %6d 0.16 %5.3f %6d\n'%(ion,direction, redshift, (Nlog10lgrid/nN10lgrid), nN10lgrid, (Nlog05lgrid/nN05lgrid), nN05lgrid, (Nlog025lgrid/nN025lgrid), nN025lgrid, (Nlog01lgrid/nN01lgrid), nN01lgrid, (Nlog033lgrid/nN033lgrid), nN033lgrid, (Nlog016lgrid/nN016lgrid), nN016lgrid))

f.close()

write_radial_column('columnave.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid),result_den,ngrid,lgrid,0)

write_radial_column('columnave_lin.%s.%s.%s.l%3.1f.dat'%(snapname,ion,direction,lgrid),10**result_den,ngrid,lgrid,1)


#    if(ion=="temperature")
#    if(ion=="Z"):
#        coldens.make_colourmap('coldens.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, colmin, colmax, "Z", ngrid, redshift)
#        print "result_num= ", result_num
#        print "result_den= ", result_den
#    if(ion=="To6"):
#        coldens.make_colourmap('coldens.%s.%s.%s.l%3.1f.png'%(snapname,ion,direction,lgrid), result_num-result_den, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, colmin, colmax, "T$_{OVI}$", ngrid, redshift)


#if(direction=="z"):
#    result_den = coldens.main(coords[indexes], hsmooth[indexes], denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.z.l%3.1f.png'%(snapname,ion,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=0,phi=0,redshift=redshift)
#if(direction=="y"):
#    result_den = coldens.main(coords[indexes], hsmooth[indexes], denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.y.l%3.1f.png'%(snapname,ion,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=90,phi=0,redshift=redshift)
#if(direction=="x"):
#    coldens.main(coords[indexes], hsmooth[indexes], denominator, center, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize,fig_name='coldens.%s.%s.x.l%3.1f.png'%(snapname,ion,lgrid),Vmin=colmin, Vmax=colmax,ion=ionname,npix=ngrid,theta=0,phi=90,redshift=redshift)


