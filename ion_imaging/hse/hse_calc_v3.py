import sys
import tables
import h5py
import glob as glob
import os
import numpy as np
import scipy.stats as scist
sys.path.append('/cosma/home/dp004/dc-oppe1/pythonCode/galmanip/')
import readGalaxy as rg
import rotateGalaxy as rotg
import writeGalaxy as wg
import binRadial as br
import binSpherical as bs
#import galmanip.readGalaxy as rg
#import galmanip.rotateGalaxy as rotg
#import galmanip.writeGalaxy as wg
#import galmanip.binRadial as br
#import galmanip.binSpherical as bs
#import coldens_ben.coldens as coldens
#from astropy import constants as const
import matplotlib.pyplot as plt
import hse_spherical 

G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24

def vels_calc(coords,vels,coords_center,vels_center):
    coords = coords-coords_center
    vels = vels-vels_center
    r = np.sqrt(coords[:,0]**2+coords[:,1]**2+coords[:,2]**2)

    theta = np.arccos(coords[:,2]/r)
    phi = np.arctan2((coords[:,1]),(coords[:,0]))

    vrad = np.sin(theta)*np.cos(phi)*vels[:,0] + np.sin(theta)*np.sin(phi)*vels[:,1] + np.cos(theta)*vels[:,2]
    vrad_mag = np.sqrt(vrad**2)

    vtan = np.cross(vels,coords)
    vtan[:,0] = vtan[:,0]/r
    vtan[:,1] = vtan[:,1]/r
    vtan[:,2] = vtan[:,2]/r
    
    vtan_mag = np.sqrt(vtan[:,0]**2+vtan[:,1]**2+vtan[:,2]**2)
    #vtangantial = np.sqrt(np.mean(vtan[:,0])**2+np.mean(vtan[:,1])**2+np.mean(vtan[:,2])**2)

    return(vrad_mag,vtan_mag,vrad,vtan)

def massfrac_radial_plot(rbins,DM_mass_bin,gas_mass_bin,stars_mass_bin,cool_mass_bin,hot_mass_bin,cool_noISM_mass_bin,galidstr,M200,R200):

    fig_mfrac = plt.figure(figsize=(5.0,4.0))
    ax_mfrac = fig_mfrac.add_subplot(111)

    if(rbins[0] < 0):
        rbins = 10**rbins
    rbins_cgs = rbins*cmpermpc
    rbins_M = np.zeros(len(DM_mass_bin))
    
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]

    fmass = open("mass_%s.dat"%galidstr,"w")

    fmass.write("#rbins DMmass gasmass starmass coolmass hotmass coolnoISMmass\n")

    for i in xrange(len(rbins_M)):    
        fmass.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e %5.3e  %5.3e\n"%(np.log10(rbins_M[i]/cmpermpc),DM_mass_bin[i],gas_mass_bin[i],stars_mass_bin[i],cool_mass_bin[i],hot_mass_bin[i],cool_noISM_mass_bin[i]))

    fmass.close()



    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(DM_mass_bin+gas_mass_bin+stars_mass_bin), color = 'black', lw=4,label="$M_{\mathrm{tot}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(stars_mass_bin), color = 'blue', lw=4,label="$M_{\mathrm{stars}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(gas_mass_bin), color = 'red', lw=4,label="$M_{\mathrm{gas}}$")

    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(gas_mass_bin*0.307/0.04825), color = 'red', lw=2, ls=":",label="$M_{\mathrm{gas}} \Omega_{\mathrm{M}}/\Omega_{\mathrm{b}}$")

    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(hot_mass_bin), color = 'orange', lw=2,label="$M_{\mathrm{hot}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(cool_mass_bin), color = 'cyan', lw=2,label="$M_{\mathrm{cool}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(cool_noISM_mass_bin), color = 'teal', ls=":", lw=2,label="$M_{\mathrm{cool-noISM}}$")

    M200_DM = np.log10( 10**M200 * (0.307-0.04825)/0.307)
    ax_mfrac.plot([-10,20], [M200_DM,M200_DM], color = 'black', lw=2,ls=':')
    ax_mfrac.plot([np.log10(R200),np.log10(R200)], [-10,20], color = 'black', lw=2, ls=':')

    ax_mfrac.set_xlim(-3.0,0.6)
    ax_mfrac.set_ylim(8.0,16.0)
    ax_mfrac.set_xlabel('log r [Mpc]', fontsize=16)
    ax_mfrac.set_ylabel('log M [$M_{\odot}$]', fontsize=16)
    ax_mfrac.legend(loc="upper left", ncol=2)
    
    fig_mfrac.subplots_adjust(bottom=0.15,left=0.15,top=0.98,right=0.98)
    fig_mfrac.savefig("mfrac_%s.png"%(galidstr))
            

def vels_radial_plot(rbins, vrad_mag_bin, vrad_bin, vtan_mag_bin, vtan_bin, galidstr, prefix, M200, R200, v200):

    fig_vel = plt.figure(figsize=(5.0,4.0))
    ax_vel = fig_vel.add_subplot(111)

    if(rbins[0] < 0):
        rbins = 10**rbins
    rbins_cgs = rbins*cmpermpc
    rbins_M = np.zeros(len(vrad_bin))
    
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]

    fvel = open("%s_%s.dat"%(prefix,galidstr),"w")

    fvel.write("#rbins sigrad vrad vtan sigtan\n")

    for i in xrange(len(rbins_M)):    
        fvel.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e\n"%(np.log10(rbins_M[i]/cmpermpc),vrad_mag_bin[i],vrad_bin[i],vtan_bin[i], vtan_mag_bin[i]))

    fvel.close()


    ax_vel.plot(np.log10(rbins_M/cmpermpc), np.log10(vrad_mag_bin), color = 'cyan', lw=2,label="$\sigma_{\mathrm{rad}}$")

    ax_vel.plot(np.log10(rbins_M/cmpermpc), np.log10(-vrad_bin), color = 'blue', lw=2,label="$v_{\mathrm{rad}}$")
    ax_vel.plot(np.log10(rbins_M/cmpermpc), np.log10(vrad_bin), color = 'blue', lw=2, ls=':')

    ax_vel.plot(np.log10(rbins_M/cmpermpc), np.log10(vtan_mag_bin), color = 'gold', lw=2,label="$\sigma_{\mathrm{tan}}$")

    ax_vel.plot(np.log10(rbins_M/cmpermpc), np.log10(vtan_bin), color = 'lime', lw=2,label="$v_{\mathrm{tan}}$")


    ax_vel.plot([-10,10], [np.log10(v200),np.log10(v200)], color = 'black', lw=2, ls=':')
    ax_vel.plot([np.log10(R200),np.log10(R200)], [-10,10], color = 'black', lw=2, ls=':')


    ax_vel.set_xlim(-3.0,0.6)
    ax_vel.set_ylim(0.5,3.5)
    ax_vel.set_xlabel('log r [Mpc]', fontsize=16)
    ax_vel.set_ylabel('log v [km s$^{-1}]$', fontsize=16)
    ax_vel.legend(loc="upper left", ncol=2)
    
    fig_vel.subplots_adjust(bottom=0.15,left=0.15,top=0.98,right=0.98)
    fig_vel.savefig("%s_%s.png"%(prefix,galidstr))
            

def radial_profile_calc(nH, T, Z, logR, Rloglow, Rloghi, nbins):


    XH = 0.752
    mu = 0.59
    nenh = 1.16
    keVtoK = 8.6e-08
    
    P = nH/XH/mu*T
    P_bin, bins = br.radialphysical(P,logR,Rloglow,Rloghi,nbins)

    S = T*keVtoK*(nH*nenh)**(-2/3.)
    S_bin = br.radialphysical(S,logR,Rloglow,Rloghi,nbins)[0]
    T_bin = br.radialphysical(T,logR,Rloglow,Rloghi,nbins)[0]
    nH_bin = br.radialphysical(nH,logR,Rloglow,Rloghi,nbins)[0]
    Z_bin = br.radialphysical(Z,logR,Rloglow,Rloghi,nbins)[0] 
    
    return(T_bin, nH_bin, Z_bin, P_bin, S_bin, bins)

#def energy_radial_plot(rbins,vradial_mag_bin,vtangential_mag_bin,T_bin,totmass_bin,M200,R200):
def energy_radial_plot(rbins,Ekinetic,Ethermal,totmass_bin,gasmass_bin,galidstr,M200,R200):


    mu = 0.59

    if(rbins[0] < 0):
        rbins = 10**rbins
    rbins_cgs = rbins*cmpermpc
    rbins_M = np.zeros(len(Ekinetic))
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]


    print "LENGTHS: rbins_M= ", len(rbins_M), "totmass, gasmass ", totmass_bin, gasmass_bin
    #Ekinetic = 0.5*(vradial_mag_bin**2 + vtangential_mag_bin**2)*1e+10
    #Ethermal = 3/2.*K_Boltz*T_bin/M_P/mu # ADD MU
    Egravity = 3/5.*G_Grav*(totmass_bin*gasmass_bin)/rbins_M ### BDO- this is wrong, but maybe works at R200.  
    kin_ther_ratio = Ekinetic/Ethermal
    Egravcarry = 0
    for i in xrange(len(rbins_M)-1):
        if(i > 0):
            Egravcarry += G_Grav*totmass_bin[i]*(gasmass_bin[i]-gasmass_bin[i-1])/((rbins_cgs[i]+rbins_cgs[i-1])/2) 
            Egravity[i] = Egravcarry 
    
    print "Ekinetic = ", np.log10(Ekinetic)
    print "Ethermal = ", np.log10(Ethermal)
    print "Egravity = ", np.log10(Egravity)

    Vir_ratio = (Ekinetic+Ethermal)/Egravity

    fegy = open("E_%s.dat"%galidstr,"w")

    fegy.write("#rbins Ekin Etherm Egrav\n")

    for i in xrange(len(rbins_M)):    
        fegy.write("% 5.3f % 5.3e % 5.3e % 5.3e\n"%(np.log10(rbins_M[i]/cmpermpc),Ekinetic[i],Ethermal[i],Egravity[i])) 

    fegy.close()
                       
    fig_E = plt.figure(figsize=(5.0,4.0))
    ax_E = fig_E.add_subplot(111)

    ax_E.plot(np.log10(rbins_M/cmpermpc), np.log10(Egravity), color = 'black', lw=4,label="$E_{\mathrm{gravity}}$")
    ax_E.plot(np.log10(rbins_M/cmpermpc), np.log10(Ethermal), color = 'red', lw=4,label="$E_{\mathrm{thermal}}$")
    ax_E.plot(np.log10(rbins_M/cmpermpc), np.log10(Ekinetic), color = 'blue', lw=4,label="$E_{\mathrm{kinetic}}$")

    ax_E.plot([np.log10(R200),np.log10(R200)], [-10,20], color = 'black', lw=2, ls=':')

    ax_E.set_xlim(-3.0,0.6)
    ax_E.set_ylim(54,64)
    ax_E.set_xlabel('log R [Mpc]', fontsize=16)
    ax_E.set_ylabel('log E [erg]', fontsize=16)
    ax_E.legend(loc="upper left", ncol=2)

    Ebindgashalo = np.log10(3/5.*G_Grav*10**M200*10**M200*M_Solar*M_Solar*0.04825/0.307/(R200*cmpermpc))
    print "Ebindgashalo= ", Ebindgashalo 
    ax_E.plot([-10,10], [Ebindgashalo,Ebindgashalo], color = 'black', lw=2, ls=':')
    ax_E.plot([np.log10(R200),np.log10(R200)], [-10,100], color = 'black', lw=2, ls=':')

    fig_E.subplots_adjust(bottom=0.15,left=0.15,top=0.98,right=0.98)
    fig_E.savefig("E_%s.png"%(galidstr))

    fig_Erat = plt.figure(figsize=(5.0,4.0))
    ax_Erat = fig_Erat.add_subplot(111)

    ax_Erat.plot(np.log10(rbins_M/cmpermpc), Vir_ratio, color = 'black', lw=4,label="($E_{\mathrm{kinetic}}$+$E_{\mathrm{thermal}}$)/$E_{\mathrm{gravity}}$")
    ax_Erat.plot(np.log10(rbins_M/cmpermpc), kin_ther_ratio, color = 'blueviolet', lw=4,label="$E_{\mathrm{kinetic}}$/$E_{\mathrm{thermal}}$")

    ax_Erat.plot([np.log10(R200),np.log10(R200)], [-10,20], color = 'black', lw=2, ls=':')
    ax_Erat.plot([-10,10], [1,1], color = 'black', lw=2, ls=':')

    ax_Erat.set_xlim(-3.0,0.6)
    ax_Erat.set_ylim(0,2.0)
    ax_Erat.set_xlabel('log R [Mpc]', fontsize=16)
    ax_Erat.set_ylabel('Ratio', fontsize=16)
    ax_Erat.legend(loc="upper left", ncol=2)

    fig_Erat.subplots_adjust(bottom=0.15,left=0.15,top=0.98,right=0.98)
    fig_Erat.savefig("Erat_%s.png"%(galidstr))


def calc_hydrostatic_masses(rbins, totmass, P, nH, vrad, vrad_mag, vtan, vtan_mag, accel, galidstr, M200, R200):

    XH = 0.752
    if(rbins[0] < 0):
        rbins = 10**rbins
    rbins_cgs = rbins*cmpermpc
    Mtot = np.zeros(len(P)+1)
    Mtherm = np.zeros(len(P)+1)
    Mrot = np.zeros(len(P)+1)
    Mraddisp = np.zeros(len(P)+1)
    Mcurl = np.zeros(len(P)+1)
    Mstream = np.zeros(len(P)+1)
    Msum = np.zeros(len(P)+1)
    Macc = np.zeros(len(P)+1)
    Mstream_pred = np.zeros(len(P)+1)

    rbins_M = np.zeros(len(P)+1)
    
    fhse = open("hse_%s.dat"%galidstr,"w")
    fhse.write("#rbins Mtot Mtherm Mrot Mstream Macc Msum Mcurl Mraddisp Mstream-pred\n")

    fig_hse = plt.figure(figsize=(5.0,4.0))
    ax_hse = fig_hse.add_subplot(111)

    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]
    
    for i in xrange(len(totmass)):
        if ((i > 0) & (i < len(totmass)-1)):
            Mtot[i] = totmass[i]*M_Solar
            Mtherm[i] = -1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)/(nH[i]/XH*M_P)*(P[i-1]-P[i+1])*K_Boltz/(rbins_M[i-1]-rbins_M[i+1])
            Mcurl[i] = 1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)*(vtan[i]*1e+05)**2/rbins_M[i]
            Mraddisp[i] = 1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)*(vrad_mag[i]*1e+05)**2/rbins_M[i]
            Mrot[i] = 1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)*(vtan_mag[i]*1e+05)**2/rbins_M[i]
            Mstream[i] = -1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)*(vrad[i]*1e+05)*((vrad[i-1]-vrad[i+1])*1e+05)/(rbins_M[i-1]-rbins_M[i+1])
            Macc[i] = -1/(4*np.pi*G_Grav)*(4*np.pi*rbins_M[i]**2)*accel[i]
            Msum[i] = Mtherm[i]+Macc[i]+Mrot[i]  #+Mstream[i]
            Mstream_pred[i] = Mtot[i]-(Mtherm[i]+Macc[i]+Mrot[i])
            #Mtherm[i] = float(-1/(4*np.pi*const.G.cgs)*(4*np.pi*rbins_M[i]**2)/(nH[i]/XH*const.m_p.cgs)*(P[i-1]-P[i+1])*const.k_B.cgs/(rbins_M[i-1]-rbins_M[i+1]))
            #print "Mtherm= ", rbins[i+1], rbins[i], Mtherm[i], np.log10(Mtot[i]/M_Solar), np.log10(Mtherm[i]/M_Solar), np.log10(Mrot[i]/M_Solar), np.log10(Mstream[i]/M_Solar)
            fhse.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e   % 4.2e\n"%(np.log10(rbins_M[i]/cmpermpc), Mtot[i]/M_Solar, Mtherm[i]/M_Solar, Mrot[i]/M_Solar, Mstream[i]/M_Solar, Macc[i]/M_Solar, Msum[i]/M_Solar, Mcurl[i]/M_Solar, Mraddisp[i]/M_Solar, Mstream_pred[i]/M_Solar,accel[i]))

    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mtot/M_Solar), color = 'black', lw=4,label="$M_{\mathrm{tot}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Msum/M_Solar), color = 'orange', lw=2, label="$M_{\mathrm{sum}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mtherm/M_Solar), color = 'red', lw=2,label="$M_{\mathrm{therm}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(-Mtherm/M_Solar), color = 'red', lw=2, ls=':')
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mrot/M_Solar), color = 'lime', lw=2,label="$M_{\mathrm{rot}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Macc/M_Solar), color = 'magenta', lw=2,label="$M_{\mathrm{acc}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(-Macc/M_Solar), color = 'magenta', lw=2,ls=":")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mstream/M_Solar), color = 'blue', lw=2,label="$M_{\mathrm{stream}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(-Mstream/M_Solar), color = 'blue', lw=2,ls=":")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mstream_pred/M_Solar), color = 'blueviolet', lw=2,label="$M_{\mathrm{stream-pred}}$")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(-Mstream_pred/M_Solar), color = 'blueviolet', lw=2,ls=":")
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mraddisp/M_Solar), color = 'cyan', lw=2,label="$M_{\mathrm{rad-disp}}$",alpha=0.3)
    ax_hse.plot(np.log10(rbins_M/cmpermpc), np.log10(Mcurl/M_Solar), color = 'gold', lw=2,label="$M_{\mathrm{curl}}$",alpha=0.3)

    ax_hse.plot([-10,10], [M200,M200], color = 'black', lw=2,ls=':')
    ax_hse.plot([np.log10(R200),np.log10(R200)], [-10,20], color = 'black', lw=2, ls=':')
    
    ax_hse.set_xlim(-3.0,0.6)
    ax_hse.set_ylim(8.0,16.0)
    ax_hse.set_xlabel('log r [Mpc]', fontsize=16)
    ax_hse.set_ylabel('log M [$M_{\odot}$]', fontsize=16)
    ax_hse.legend(loc="upper left", ncol=3)

    fig_hse.subplots_adjust(bottom=0.15,left=0.15,top=0.98,right=0.98)
    fig_hse.savefig("hse_%s.png"%galidstr)
                
    fhse.close()
    return(Mtot,Mtherm)


####  BEGINNING OF MAIN PROGRAM ###


snap_intag = sys.argv[1] # snapshot identifier
galidstr = sys.argv[2] # Galaxy ID string-  has a specific format that is parsed below for halo mass, stellar mass, SFR, etc.  
    #snap_outtag = sys.argv[2]
    #sim_name = sys.argv[2]  # e.g. m1e11f30/data
x_gal = float(sys.argv[3]) # x-position of galaxy in cosmological box in comoving Mpc/h
y_gal = float(sys.argv[4]) 
z_gal = float(sys.argv[5])
v_x_gal = float(sys.argv[6]) # x velocity of galaxy in cosmological box in km/s  
v_y_gal = float(sys.argv[7])
v_z_gal = float(sys.argv[8])

write_hdf5 = 0
metals = 1
chem = 0 # tag to read in chemical ions, turn off to 0  
dofigs = 0 # tag to make column density map figures.  
###box = 0 # replaced if hash table
cosmoowls = 0 # this was to read a certain type of simulation.   
###acceleration = 1
do_angmom_calcs = 1 # Do angular momentum calculation.  
do_hse_calcs = 0


runlabel = 'halo' # read specific gal id string to get halo information.  

unit_mass_in_cgs = 1.989e33 * 1.0e10 
unit_length_in_cgs = 3.0857e24 
proton_mass_cgs = 1.67e-24

unit_Density_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3/proton_mass_cgs 


if (runlabel=='halo'): # Parsing of galidstr  
    haloname = galidstr.split("halo")[1]
    halo_id = haloname.split("_")[1]
    M200 = float(haloname.split("_")[2])
    ms = float(haloname.split("_")[3])
    sfr = float(haloname.split("_")[4])
    haloinfostr = 'lg M$_{200}=%4.1f$, lg M$_{*}=%3.1f$, SFR$=%4.2f$'%(M200,ms,sfr)
    z_hold = galidstr.split("_z")[1]
    zname = z_hold.split("_")[0]
    zint = zname.split("p")[0]
    zfrac = zname.split("p")[1]
    redshift = float(zint)*1.0 + float(zfrac)/1000.        
else:
    haloinfostr = ''


data_dir = "" #% sim_name
sniptag = "SNAP"  # These are snapshots.  
sim = "." # /projects/beop5934/CosmoSim/[simname]

hdf5filename0 = sim + "/snapshot_" + snap_intag + "/snap_" + snap_intag + ".0.hdf5"
print "hdf5filename0 = ", hdf5filename0

hashtable =  "/HashTable/PartType0/FirstKeyInFile" in  h5py.File(hdf5filename0) # This alerts the code if we are selecting from a periodic box instead of a zoom.  

#h5file_inbase = "%ssnapshot_%s/snap_%s" % (data_dir, snap_intag, snap_intag)

hubbleparam,aex,redshift,boxsize,mass_table = rg.header(sniptag,sim,snap_intag) # read header

omegaM = 0.307
omegaL = 1-omegaM
omegaratio = (omegaM+omegaL/(1+redshift)**3) # proportional to the Hubble expansion at the given redshift, important for you Ezra.  :) 
R200 = 1.63e-5*(10**float(M200)*hubbleparam)**0.333/omegaratio**0.333/(1+redshift)/hubbleparam # R200 is the virial radius in Mpc, 200x the critcal overdensity
v200 = np.sqrt(G_Grav*10**M200*M_Solar/(R200*cmpermpc))/1e+05 # v200 is the virial velocity in km/s.  

boxsize_pMpc = boxsize*aex/hubbleparam # Even zooms live in periodic boxes, so this gives the boxsize of the simulation, in Mpc.  
print "boxsize_pMpc = ", boxsize_pMpc

# aex = 1/(1+redshift)- the scale of the Universe relative today- smaller at high-z, hubbleparm is H_0/100(km/s/Mpc)- and is 0.6777 since H_0=67.77 for our simulations
gal_coords_pMpc = ((float(x_gal*aex/hubbleparam),float(y_gal*aex/hubbleparam),float(z_gal*aex/hubbleparam))) # Headaches can be had converting between comoving Mpc, physical Mpc, and physical Mpc/h.  We'll discuss this in detail.  
gal_vels_kms = ((float(v_x_gal),float(v_y_gal),float(v_z_gal)))

#We are going to shift the galaxy to the center of the box.  Our zooms are selected from a 25 or 100 Mpc box.  We shift to the center, but we use Mpc/h, which is 25/2*0.6777 = 8.47125 Mpc.  
boxcenter_pMpc = ((float(boxsize_pMpc/2.),float(boxsize_pMpc/2.),float(boxsize_pMpc/2.))) 

grabregion = 3.0*R200/aex*hubbleparam # We are going to grab only gas, DM, and stars within 3 times the virial radius.  
xlow = np.where(x_gal - grabregion > 0, x_gal - grabregion, 0.0)
xhi = np.where(x_gal + grabregion < boxsize, x_gal + grabregion, boxsize)
ylow = np.where(y_gal - grabregion > 0, y_gal - grabregion, 0.0)
yhi = np.where(y_gal + grabregion < boxsize, y_gal + grabregion, boxsize)
zlow = np.where(z_gal - grabregion > 0, z_gal - grabregion, 0.0)
zhi = np.where(z_gal + grabregion < boxsize, z_gal + grabregion, boxsize)

#xlow = x_gal - grabregion
#xhi = x_gal + grabregion
#ylow = y_gal - grabregion
#yhi = y_gal + grabregion
#zlow = z_gal - grabregion
#zhi = z_gal + grabregion
#hashtable = 0

#Reading in particles starting with dark matter (DM) using my readGalaxy (rg) library
if(hashtable):
    print "Using hastables" 
    coords_DM,vels_DM = rg.coords(sniptag,sim,snap_intag,1,region=[xlow,xhi,ylow,yhi,zlow,zhi])
else:
    print "Not using hastables" 
    coords_DM,vels_DM = rg.coords(sniptag,sim,snap_intag,1)

print "length = ",len(coords_DM[:,0]),len(vels_DM[:,0]) 
#mass_DM = rg.mass(sniptag,sim,snap_intag,1)

if(hashtable):
    coords_stars,vels_stars = rg.coords(sniptag,sim,snap_intag,4,region=[xlow,xhi,ylow,yhi,zlow,zhi])
else:
    coords_stars,vels_stars = rg.coords(sniptag,sim,snap_intag,4)
    
if(hashtable):
    mass_stars = rg.mass(sniptag,sim,snap_intag,4,region=[xlow,xhi,ylow,yhi,zlow,zhi])
    Z_stars,stellarformationtime_stars = rg.starprops(sniptag,sim,snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi]) # 8/24/18 this is broken.  
    ##print "len(stars)= ", len(mass_stars), len(Z_stars)
else:
    mass_stars = rg.mass(sniptag,sim,snap_intag,4)
    Z_stars,stellarformationtime_stars = rg.starprops(sniptag,sim,snap_intag)
                                                      
coords_BH, vels_BH, mass_BH = rg.BH_zoom_read(sniptag,sim,snap_intag)

# Reading in gas coords, vels, mass, and physical parameters- nH is hydrogen numebr density in terms of cm^-3, T is temperature in Kelvin, Z is metallicity
if(hashtable):    
    coords_gas,vels_gas = rg.coords(sniptag,sim,snap_intag,0,region=[xlow,xhi,ylow,yhi,zlow,zhi])
else:
    coords_gas,vels_gas = rg.coords(sniptag,sim,snap_intag,0)

if(hashtable):    
    mass_gas = rg.mass(sniptag,sim,snap_intag,0,region=[xlow,xhi,ylow,yhi,zlow,zhi])
else:
    mass_gas = rg.mass(sniptag,sim,snap_intag,0)
    
if(cosmoowls):
    if(hashtable):
        density, T = rg.gasprops_cosmoowls(sniptag,sim,snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi])
    else:
        density, T = rg.gasprops_cosmoowls(sniptag,sim,snap_intag)

    SFR = T * 0
    nH = density*unit_Density_in_cgs*0.752
else:
    if(hashtable):
        density, T, Z, SFR, hsmooth = rg.gasprops(sniptag,sim,snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi])
        hydrogen = rg.hydrogen(sniptag,sim,snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi])
        if(metals):
            hydrogen,carbon,oxygen,magnesium = rg.metalsprimary(sniptag, sim, snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi])
            helium,nitrogen,neon,silicon,iron = rg.metalssecondary(sniptag, sim, snap_intag,region=[xlow,xhi,ylow,yhi,zlow,zhi])
    else:        
        density, T, Z, SFR, hsmooth = rg.gasprops(sniptag,sim,snap_intag)
        hydrogen = rg.hydrogen(sniptag,sim,snap_intag)
        if(metals):
            hydrogen,carbon,oxygen,magnesium = rg.metalsprimary(sniptag, sim, snap_intag)
            helium,nitrogen,neon,silicon,iron = rg.metalssecondary(sniptag, sim, snap_intag)

    nH = density*unit_Density_in_cgs*hydrogen
                
                
### Finding the indexes of stars within 0.03 Mpc (30 kpc) and centering them on the galaxy coordinates.  
angmom_indexes = rg.distanceindexes_limitedshift(coords_stars,gal_coords_pMpc,0.03,boxsize_pMpc) # Now Angular momentum axis is 30 kpc (not 30/h kpc) 5-18-18
print "angmom_array = ", coords_stars[angmom_indexes], "length = ", len(coords_stars[angmom_indexes])
# My rotatge galaxy (rotg) library calculates the angular momentum vector and then calculates the Euler angle by which the galaxy and velocity coordinate need to be rotated so the the angular momentum vector is along z direction.  
phi, theta, psi, normed_axis = rotg.calc_angmom_axisangles(gal_coords_pMpc, gal_vels_kms, coords_stars[angmom_indexes], vels_stars[angmom_indexes], boxsize_pMpc)



# Older code when I used cool gas at within 100 kpc of the galaxy center to calculate angular momentum.  
###cool_gas = np.where((T<1e+05), mass_gas, 0.0)
#### With Gas
###cool_indexes = np.where((T<1e+05))
###coords_cool = coords_gas[cool_indexes]
###vels_cool = vels_gas[cool_indexes]

###angmom_indexes = rg.distanceindexes(coords_cool,gal_coords_pMpc,0.10/hubbleparam,boxsize_pMpc)
###print "angmom_array = ", coords_cool[angmom_indexes], "length = ", len(coords_cool[angmom_indexes])
###phi, theta, psi = rotg.calc_angmom_axisangles(gal_coords_pMpc, gal_vels_kms, coords_cool[angmom_indexes], vels_cool[angmom_indexes], boxsize_pMpc)

print "Finished ang mom calculation "

Rcut = 3.0*R200 # New as of 11/1/17- BDO 1/23/18 this is sort of redundant with grabregion code from above.  
###Rcut = 4.0 # Used to us 4 Mpc.  

# Grab indexes for gas/DM/stars/BHs within Rcut  
gas_indexes = rg.distanceindexes_limitedshift(coords_gas,gal_coords_pMpc,Rcut,boxsize_pMpc)
print "Finished gas indexing, gas_indexes length = ", len(coords_gas[gas_indexes]) 
DM_indexes = rg.distanceindexes_limitedshift(coords_DM,gal_coords_pMpc,Rcut,boxsize_pMpc)
print "Finished DM indexing, DM_indexes length = ", len(coords_DM[DM_indexes]) 
#quit()
stars_indexes = rg.distanceindexes_limitedshift(coords_stars,gal_coords_pMpc,Rcut,boxsize_pMpc)
print "Finished stellar indexing, star_indexes length = ", len(coords_stars[stars_indexes]) 
BH_indexes = rg.distanceindexes_limitedshift(coords_BH,gal_coords_pMpc,Rcut,boxsize_pMpc)
print "Finished BH indexing, BH_indexes length = ", len(coords_BH[BH_indexes]) 

# rotated coordinates and velocities using angles calculated from angular momentum calculation.  
coords_gas_rot, vels_gas_rot = rotg.rotategalaxy(coords_gas[gas_indexes],vels_gas[gas_indexes],gal_coords_pMpc,gal_vels_kms,boxsize_pMpc,phi,theta,psi)
coords_stars_rot, vels_stars_rot = rotg.rotategalaxy(coords_stars[stars_indexes],vels_stars[stars_indexes],gal_coords_pMpc,gal_vels_kms,boxsize_pMpc,phi,theta,psi)
coords_DM_rot, vels_DM_rot = rotg.rotategalaxy(coords_DM[DM_indexes],vels_DM[DM_indexes],gal_coords_pMpc,gal_vels_kms,boxsize_pMpc,phi,theta,psi)
coords_BH_rot, vels_BH_rot = rotg.rotategalaxy(coords_BH[BH_indexes],vels_BH[BH_indexes],gal_coords_pMpc,gal_vels_kms,boxsize_pMpc,phi,theta,psi)


#Gas temepratures
Tgas = T[gas_indexes]
#Cool gas is anything that is <10^5 K.  
cool_indexes = np.where((Tgas<1e+05))
coords_cool_rot = coords_gas_rot[cool_indexes]
vels_cool_rot = vels_gas_rot[cool_indexes]
mass_gas_cool = mass_gas[cool_indexes]

hot_indexes = np.where((Tgas>=1e+05))
coords_hot_rot = coords_gas_rot[hot_indexes]
vels_hot_rot = vels_gas_rot[hot_indexes]
mass_gas_hot = mass_gas[hot_indexes]


print "Finished rotating galaxy"

# Now shifting the galaxy to the center of the box and 0 velocity using rotg library  
coords_gas_shift, vels_gas_shift = rotg.shiftgalaxy(coords_gas_rot,vels_gas_rot,gal_coords_pMpc,gal_vels_kms,boxcenter_pMpc,[0,0,0],boxsize_pMpc)
coords_DM_shift, vels_DM_shift = rotg.shiftgalaxy(coords_DM_rot,vels_DM_rot,gal_coords_pMpc,gal_vels_kms,boxcenter_pMpc,[0,0,0],boxsize_pMpc)
coords_stars_shift, vels_stars_shift = rotg.shiftgalaxy(coords_stars_rot,vels_stars_rot,gal_coords_pMpc,gal_vels_kms,boxcenter_pMpc,[0,0,0],boxsize_pMpc)
coords_BH_shift, vels_BH_shift = rotg.shiftgalaxy(coords_BH_rot,vels_BH_rot,gal_coords_pMpc,gal_vels_kms,boxcenter_pMpc,[0,0,0],boxsize_pMpc)

print "Finished shifting galaxy"

# This is to test again angular momentum indices and if the galaxy is rotated properly.  Just print some stuff out.  
print "coords_stars_rot= ", coords_stars_rot
angmom_indexes = rg.distanceindexes_limitedshift(coords_stars_rot,gal_coords_pMpc,0.03/hubbleparam,boxsize_pMpc)
print "angmom_indexes = ", coords_stars_rot[angmom_indexes], "length = ", len(coords_stars_rot[angmom_indexes])
print "mean x,y,z = ", np.mean(coords_stars_rot[angmom_indexes,0]),np.mean(coords_stars_rot[angmom_indexes,1]),np.mean(coords_stars_rot[angmom_indexes,2])
        
if (cosmoowls): # hopefully code that you shouldn't need to worry about.  
    if(write_hdf5):
        wg.gas("snap_%s.0.hdf5"%galidstr, coords_gas_shift,vels_gas_shift,mass_gas[gas_indexes],nH[gas_indexes]/0.752/unit_Density_in_cgs,T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],aex,redshift,boxsize,hubbleparam,mass_table,cosmoowls)
    #print "hsmooth_read = ", hsmooth, len(hsmooth)
    hsmooth = (nH[gas_indexes]/0.752/unit_Density_in_cgs/mass_gas[gas_indexes])**(-0.3333)*0.00111066
    print "hsmooth_calc = ", hsmooth, len(hsmooth)
    hydrogen = nH[gas_indexes]*0.0 + 0.752
    #SFR = nH[gas_indexes]*0.0 
else: # this is code to read in the metal abundances.          
    # I write out a snapshot of rotated and centered gas particles for other uses.  
    if(write_hdf5):
        wg.gas("snap_%s.0.hdf5"%galidstr, coords_gas_shift,vels_gas_shift,mass_gas[gas_indexes],nH[gas_indexes]/hydrogen[gas_indexes]/unit_Density_in_cgs,T[gas_indexes],Z[gas_indexes],SFR[gas_indexes],hsmooth[gas_indexes],hydrogen[gas_indexes],helium[gas_indexes],carbon[gas_indexes],nitrogen[gas_indexes],oxygen[gas_indexes],neon[gas_indexes],magnesium[gas_indexes],silicon[gas_indexes],iron[gas_indexes],aex,redshift,boxsize,hubbleparam,mass_table,cosmoowls)
    ###wg.gas("snap_%s.0.hdf5"%galidstr, coords_gas_shift,vels_gas_shift,mass_gas[gas_indexes],nH[gas_indexes]/hydrogen[gas_indexes]/unit_Density_in_cgs,T[gas_indexes],Z[gas_indexes],SFR[gas_indexes],hsmooth[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],hydrogen[gas_indexes],aex,redshift,boxsize,hubbleparam,mass_table,cosmoowls)


    hsmooth = hsmooth[gas_indexes] # hsmooth is the smoothing length of SPH gas particles, i.e. the length of their smoothing kernel.  
    hydrogen = hydrogen[gas_indexes]

    if(write_hdf5):
        print "Stellar lengths= ", len(mass_stars[stars_indexes]) #, len(Z_stars[stars_indexes])
        ###wg.stars("snap_%s.0.hdf5"%galidstr,coords_stars_shift,vels_stars_shift,mass_stars[stars_indexes],Z_stars[stars_indexes],stellarformationtime_stars[stars_indexes],aex,hubbleparam)
        #wg.stars("snap_%s.0.hdf5"%galidstr,coords_stars_shift,vels_stars_shift,mass_stars,Z_stars,stellarformationtime_stars,aex,hubbleparam)

        ###wg.darkmatter("snap_%s.0.hdf5"%galidstr,coords_DM_shift,vels_DM_shift,aex,hubbleparam)

if (chem): # individual ions- don't use.  
    oxygen = oxygen[gas_indexes]

    h1,c2,c3,c4,n5,o6,o7,o8,ne8,mg2,si2,si3,si4 = rg.chemistryprimary(sniptag, sim, snap_intag)
    h1 = h1[gas_indexes]
    c2 = c2[gas_indexes]
    c3 = c3[gas_indexes]
    c4 = c4[gas_indexes]
    n5 = n5[gas_indexes]
    o6 = o6[gas_indexes]
    o7 = o7[gas_indexes]
    o8 = o8[gas_indexes]
    ne8 = ne8[gas_indexes]
    mg2 = mg2[gas_indexes]
    si2 = si2[gas_indexes]
    si3 = si3[gas_indexes]
    si4 = si4[gas_indexes]

    if(write_hdf5):
        wg.chemistryprimary("snap_%s.0.hdf5"%galidstr, h1, c2, c3, c4, n5, o6, o7, o8, ne8, mg2, si2, si3, si4)

if(write_hdf5):
    wg.blackholes("snap_%s.0.hdf5"%galidstr,coords_BH_shift,vels_BH_shift,mass_BH[BH_indexes],aex,hubbleparam)

mass_stars = mass_stars[stars_indexes]
###Z_stars = Z_stars[stars_indexes]
###stellarformationtime_stars = stellarformationtime_stars[stars_indexes]
mass_gas = mass_gas[gas_indexes]
mass_BH = mass_BH[BH_indexes]
nH = nH[gas_indexes]
T = T[gas_indexes]
if(cosmoowls):
    Z = T
    SFR = T
else:
    Z = Z[gas_indexes]
    SFR = SFR[gas_indexes] # Star formation rate
P = T*nH/(0.752*0.59) # Pressure- 0.59 in mu- mean molecular weight of ionized gas.  0.752 is mass fraction in hydrogran.  
#cool_gas = np.where((T<1e+05) & (SFR <= 0), mass_gas, 0.0)
#hot_gas = np.where((T>=1e+05) & (SFR <= 0), mass_gas, 0.0)
cool_gas = np.where((T<1e+05), mass_gas, 0.0) # cool gas mass mask (set mass to 0 if hot gas 
hot_gas = np.where((T>=1e+05), mass_gas, 0.0) # hot gas mass mask  

hot_indexes = np.where((T>=1e+05))
cool_indexes = np.where((T<1e+05))

fangmomaxis = open("angmomaxis_%s.dat"%galidstr,"w")
fangmomaxis.write("# 1: Dist(kpc) 2: StarMass 3: StarAngle 4: CoolMass 5: CoolAngle 6: HotMass 7: HotAngle 8: HotCoolAngle\n")

zvector = [0.,0.,1.]
if(do_angmom_calcs):
    for i in range(50):
        dist_Mpc_min = (i)*0.01
        dist_Mpc = (i+1)*0.01
        star_indexes = rg.distanceindexes_limitedshift(coords_stars_shift,boxcenter_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min)
        phi_star, theta_star, psi_star, normed_axis_star = rotg.calc_angmom_axisangles(boxcenter_pMpc, [0,0,0], coords_stars_shift[star_indexes], vels_stars_shift[star_indexes], boxsize_pMpc)
        print "angmom_stars= ", dist_Mpc, np.sum(mass_stars[star_indexes]),np.arccos(np.dot(normed_axis_star,zvector))*180./np.pi 
        print "lengths = ", len(coords_cool_rot[:,0]), len(vels_cool_rot[:,0]), len(mass_gas_cool)
        #coolgas_indexes = rg.distanceindexes(coords_gas_shift[cool_indexes],boxcenter_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min) 
        ##coolgas_indexes = rg.distanceindexes((coords_gas_shift,boxcenter_pMpc,dist_Mpc,boxsize_pMpc) & (T<1e+05))#,mindist=dist_Mpc_min) 
        coolgas_indexes = rg.distanceindexes_limitedshift(coords_cool_rot,gal_coords_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min) 
        phi_cool, theta_cool, psi_cool, normed_axis_cool = rotg.calc_angmom_axisangles(gal_coords_pMpc, gal_vels_kms, coords_cool_rot[coolgas_indexes], vels_cool_rot[coolgas_indexes], boxsize_pMpc)
        print "angmom_cool= ", dist_Mpc, np.sum(mass_gas_cool[coolgas_indexes]),np.arccos(np.dot(normed_axis_cool,zvector))*180./np.pi 
        #hotgas_indexes = rg.distanceindexes(coords_gas_shift[hot_indexes],boxcenter_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min)
        ##hotgas_indexes = rg.distanceindexes(coords_gas_shift,boxcenter_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min) 
        hotgas_indexes = rg.distanceindexes_limitedshift(coords_hot_rot,gal_coords_pMpc,dist_Mpc,boxsize_pMpc)#,mindist=dist_Mpc_min) 
        phi_hot, theta_hot, psi_hot, normed_axis_hot = rotg.calc_angmom_axisangles(gal_coords_pMpc, gal_vels_kms, coords_hot_rot[hotgas_indexes], vels_hot_rot[hotgas_indexes], boxsize_pMpc)
        print "angmom_hot= ", dist_Mpc, np.sum(mass_gas_hot[hotgas_indexes]),np.arccos(np.dot(normed_axis_hot,zvector))*180./np.pi 
        
        fangmomaxis.write("%5.1f  %5.3e  %6.1f  %5.3e  %6.1f  %5.3e  %6.1f  %6.1f\n"%(dist_Mpc*1e+03,np.sum(mass_stars[star_indexes]), np.arccos(np.dot(normed_axis_star,zvector))*180./np.pi, np.sum(mass_gas_cool[coolgas_indexes]),  np.arccos(np.dot(normed_axis_cool,zvector))*180./np.pi , np.sum(mass_gas_hot[hotgas_indexes]),  np.arccos(np.dot(normed_axis_hot,zvector))*180./np.pi, np.arccos(np.dot(normed_axis_hot,normed_axis_cool))*180./np.pi))

    fangmomaxis.close()

#    quit()



#To be honest, I used the masking method above, and I also made indexes of gas that is hot and cool.  Really should be better about consistency.  

lt104_indexes = np.where(T<1e+04) # <10^4 K
gt106_indexes = np.where(T>=1e+06) # >10^6 K.  

# calculate radial distance of  particle from center of galaxy.  
R_gas = np.sqrt((coords_gas_shift[:,0]-boxcenter_pMpc[0])**2+(coords_gas_shift[:,1]-boxcenter_pMpc[1])**2+(coords_gas_shift[:,2]-boxcenter_pMpc[2])**2)
R_DM = np.sqrt((coords_DM_shift[:,0]-boxcenter_pMpc[0])**2+(coords_DM_shift[:,1]-boxcenter_pMpc[1])**2+(coords_DM_shift[:,2]-boxcenter_pMpc[2])**2)
logR_DM = np.log10(R_DM)   
print "logR_DM = ", logR_DM
R_stars = np.sqrt((coords_stars_shift[:,0]-boxcenter_pMpc[0])**2+(coords_stars_shift[:,1]-boxcenter_pMpc[1])**2+(coords_stars_shift[:,2]-boxcenter_pMpc[2])**2)
logR_stars = np.log10(R_stars)   
print "logR_stars = ", logR_stars
R_BH = np.sqrt((coords_BH_shift[:,0]-boxcenter_pMpc[0])**2+(coords_BH_shift[:,1]-boxcenter_pMpc[1])**2+(coords_BH_shift[:,2]-boxcenter_pMpc[2])**2)

# indexes where I select different gas- noISM is no star-forming (i.e. only CGM) gas, and nosatISM is more obscure-  this isn't important for you Ezra, but I wanted to calculate gas masses using these different categories below.  
cool_nosatISM_indexes = np.where(((T<1e+05) & (R_gas<=0.03)) | ((T<1e+05) & (SFR <= 0) & (R_gas>0.03)))
cool_noISM_indexes = np.where((T<1e+05) & (SFR <= 0))

# calculating stellar mass inside R200, and inside 30 kpc, which is defined as the central galaxy.  
M_stars = np.sum(mass_stars[np.where(R_stars<R200)])
M_stars_30kpc = np.sum(mass_stars[np.where(R_stars<0.03)])

# M_DM- dark matter mass, M_CGM- gas mass that is not in ISM, in the CGM.  M_ISM is ISM, M_ISM_30kpc is ISM in central galaxy, then CGM gas masses in various temperature bins.  
M_DM =  len(R_DM[np.where(R_DM<R200)])*mass_table[1]*1e+10/hubbleparam
M_CGM = np.sum(mass_gas[np.where((R_gas<R200) & (SFR<=0))])
M_ISM = np.sum(mass_gas[np.where((R_gas<R200) & (SFR>0))])
M_ISM_30kpc = np.sum(mass_gas[np.where((R_gas<0.03) & (SFR>0))])
print "MASSES= ", np.log10(M_DM), np.log10(M_CGM), " ISMs= ", np.log10(M_ISM), np.log10(M_ISM_30kpc), " Stars= ", np.log10(M_stars),ms,np.log10(M_stars_30kpc)
M_lt104 = np.sum(mass_gas[np.where((R_gas<R200) & (T<1e+04))])
M_104_105 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+04) & (T<1e+05))])
M_105_106 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+05) & (T<1e+06))])
M_106_107 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+06) & (T<1e+07))])
M_gt107 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+07))])
M_BH = np.sum(mass_BH[np.where(R_BH<R200)])
max_M_BH = np.max(mass_BH[np.where(R_BH<R200)])

print "MASSES_gas= ", np.log10(M_lt104), np.log10(M_104_105), np.log10(M_105_106), np.log10(M_106_107), np.log10(M_gt107)

# Writing one line file of masses in the halo- useful for your work Ezra.  
fmasses = open("masses_%s.dat"%galidstr,"w")
fmasses.write("galidstr= %5.2f %5.2f %5.3f MDM= %5.3e Mstars=  %5.3e %5.3e MISMs= %5.3e %5.3e MCGM= %5.3e MCGM_sub= %5.3e %5.3e %5.3e %5.3e %5.3e M_BH= %5.3e %5.3e\n"%(M200,ms,sfr,M_DM,M_stars,M_stars_30kpc,M_ISM,M_ISM_30kpc,M_CGM,M_lt104,M_104_105,M_105_106,M_106_107,M_gt107,M_BH,max_M_BH))
fmasses.close()

#BDO- needs to write Z (metal) masses too (at some point).  

# This defines my spherical coordinate binning scheme.  I use logarithmic bins from 0.003 R_200 to 3 R_200 in an adjustable bin size, 120 log bins between these radii- have used 60 and 240 as well.   
rlow = np.log10(R200/316.)
rhi = np.log10(R200*3.16)
nradial= 120 #60
ntheta = 10 #36
nphi = 18 #20

if (M200>15.85): # if a higher mass halo- use fewer radial bins (for my 2018 paper)
    nradial= 60 

# I experimented with linear radial bins (dolin) where I would use 5 kpc or 10 kpc bins- this was too noisy.
dolin = 0
rlow_lin = 0
rhi_lin = 3
nradial_lin = 600
if(hashtable):
    nradial_lin = 300

if(dolin):
    rlow = rlow_lin
    rhi = rhi_lin
    nradial = nradial_lin
else:
    R_gas = np.log10(R_gas)
    R_DM = np.log10(R_DM)
    R_stars = np.log10(R_stars)

    
R_gas_cool = R_gas[cool_indexes]
R_gas_hot = R_gas[hot_indexes]

# I calculate radial and tangential velocities using a function at the top of this file.  I calculate mean v and rms (root mean square) v.  
vradial_mag, vtangential_mag, vradial, vtangential = vels_calc(coords_gas_shift,vels_gas_shift,boxcenter_pMpc,[0,0,0])
print "v_means (<R>,<tan>,R,tan) ", np.mean(vradial_mag), np.mean(vtangential_mag), np.mean(vradial), np.mean(vtangential)

vradial_mag_hot, vtangential_mag_hot, vradial_hot, vtangential_hot = vels_calc(coords_gas_shift[hot_indexes],vels_gas_shift[hot_indexes],boxcenter_pMpc,[0,0,0])
vradial_mag_cool, vtangential_mag_cool, vradial_cool, vtangential_cool = vels_calc(coords_gas_shift[cool_indexes],vels_gas_shift[cool_indexes],boxcenter_pMpc,[0,0,0])

# This is redundant code where I binned radially (using my binRadial library) radial and tangential velocities.   A bunch of overkill on velocities.  Really needs to be cleaned up or just ignored.  
vradial_all_bin = br.radialphysical(vradial,R_gas,rlow,rhi,nradial)[0]
vradial_mag_all_bin = br.radialphysical(vradial_mag,R_gas,rlow,rhi,nradial)[0]
vtangential_all_bin_x = br.radialphysical(vtangential[:,0],R_gas,rlow,rhi,nradial)[0]
vtangential_all_bin_y = br.radialphysical(vtangential[:,1],R_gas,rlow,rhi,nradial)[0]
vtangential_all_bin_z = br.radialphysical(vtangential[:,2],R_gas,rlow,rhi,nradial)[0]
vtangential_mag_all_bin = br.radialphysical(vtangential_mag,R_gas,rlow,rhi,nradial)[0]
vtangential_all_bin = np.sqrt(vtangential_all_bin_x**2+vtangential_all_bin_y**2+vtangential_all_bin_z**2)

vradial_hot_bin = br.radialphysical(vradial_hot,R_gas_hot,rlow,rhi,nradial)[0]
vradial_mag_hot_bin = br.radialphysical(vradial_mag_hot,R_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_x = br.radialphysical(vtangential_hot[:,0],R_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_y = br.radialphysical(vtangential_hot[:,1],R_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_z = br.radialphysical(vtangential_hot[:,2],R_gas_hot,rlow,rhi,nradial)[0]
vtangential_mag_hot_bin = br.radialphysical(vtangential_mag_hot,R_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin = np.sqrt(vtangential_hot_bin_x**2+vtangential_hot_bin_y**2+vtangential_hot_bin_z**2)

vradial_cool_bin = br.radialphysical(vradial_cool,R_gas_cool,rlow,rhi,nradial)[0]
vradial_mag_cool_bin = br.radialphysical(vradial_mag_cool,R_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_x = br.radialphysical(vtangential_cool[:,0],R_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_y = br.radialphysical(vtangential_cool[:,1],R_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_z = br.radialphysical(vtangential_cool[:,2],R_gas_cool,rlow,rhi,nradial)[0]
vtangential_mag_cool_bin = br.radialphysical(vtangential_mag_cool,R_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin = np.sqrt(vtangential_cool_bin_x**2+vtangential_cool_bin_y**2+vtangential_cool_bin_z**2)

# I calculate velocities in spherical coordinates as well.   
vradial_hot, vtheta_hot, vphi_hot = bs.sphericalvels(coords_gas_shift[hot_indexes]-boxcenter_pMpc,vels_gas_shift[hot_indexes])
vradial_cool, vtheta_cool, vphi_cool = bs.sphericalvels(coords_gas_shift[cool_indexes]-boxcenter_pMpc,vels_gas_shift[cool_indexes])
vradial_all, vtheta_all, vphi_all = bs.sphericalvels(coords_gas_shift-boxcenter_pMpc,vels_gas_shift)

mu =0.59 # mean molecular weight
XH = 0.752 # hydrogen mass fraction of primordial (or most) gas.  

# Almost nobody outputs acceleration of gas particles in simulations, but I reran them to do so, this checks if the acceleration vector exists.  
acceleration =  "/PartType0/Acceleration" in  h5py.File(hdf5filename0)

print "ACCELERATION = ", acceleration

if(acceleration):
    acc_CGSConv = 3.24078e-15
    if(hashtable):
        accel_gas = rg.accel(sniptag,sim,snap_intag,0,region=[xlow,xhi,ylow,yhi,zlow,zhi])*acc_CGSConv  
    else:
        accel_gas = rg.accel(sniptag,sim,snap_intag,0)*acc_CGSConv  

    print "L677: phi, theta, psi = ", phi, theta, psi
    accel_gas_rot = rotg.rotatevector(accel_gas[gas_indexes],phi,theta,psi)
    accel_gas_shift = accel_gas_rot
    print "ACCELERATION, accel_gas_shift = ", accel_gas_shift
    accradial_hot_mag, acctangential_hot_mag, accradial_hot, acctangential_hot = vels_calc(coords_gas_shift[hot_indexes],accel_gas_shift[hot_indexes],boxcenter_pMpc,[0,0,0])
    accradial_hot_bin = br.radialphysical(accradial_hot,R_gas_hot,rlow,rhi,nradial)[0]

    accradial_cool_mag, acctangential_cool_mag, accradial_cool, acctangential_cool = vels_calc(coords_gas_shift[cool_indexes],accel_gas_shift[cool_indexes],boxcenter_pMpc,[0,0,0])
    accradial_cool_bin = br.radialphysical(accradial_cool,R_gas_cool,rlow,rhi,nradial)[0]

    accradial_all_mag, acctangential_all_mag, accradial_all, acctangential_all = vels_calc(coords_gas_shift,accel_gas_shift,boxcenter_pMpc,[0,0,0])
    accradial_all_bin = br.radialphysical(accradial_all,R_gas,rlow,rhi,nradial)[0]
    print "ACCELERATION_bin, radial = ", accradial_all_bin

else: # I use placeholders for the acceleration vectors if they are not available.  
    accradial_hot_bin = vradial_hot_bin*0.0
    accradial_cool_bin = vradial_cool_bin*0.0
    accradial_all_bin = vradial_all_bin*0.0
    accradial_hot = np.zeros(len(P[hot_indexes]))
    accradial_cool = np.zeros(len(P[cool_indexes]))
    accradial_all = np.zeros(len(P))

#My binRadial library has a radial cumulative summing in the radial bins.  
gas_mass_radial = br.radialcumsum(mass_gas,R_gas,rlow,rhi,nradial)[0]
DM_mass_radial = br.radialcumsum(coords_DM[:,0]*0.0+mass_table[1]*1e+10/hubbleparam,R_DM,rlow,rhi,nradial)[0]
stars_mass_radial = br.radialcumsum(mass_stars,R_stars,rlow,rhi,nradial)[0]
Totmass_radial = (gas_mass_radial+DM_mass_radial+stars_mass_radial)*M_Solar

# These function calls to radial_profile_calc (function above) calculate some physical parameters as a function of radius for various gas cuts- P is pressure and S is astronomer's entropy.  
T_hot_bin, nH_hot_bin, Z_hot_bin, P_hot_bin, S_hot_bin, rbins = radial_profile_calc(nH[hot_indexes], T[hot_indexes], Z[hot_indexes], R_gas[hot_indexes],rlow,rhi,nradial)
T_cool_bin, nH_cool_bin, Z_cool_bin, P_cool_bin, S_cool_bin, rbins = radial_profile_calc(nH[cool_indexes], T[cool_indexes], Z[cool_indexes], R_gas[cool_indexes],rlow,rhi,nradial)
T_bin, nH_bin, Z_bin, P_bin, S_bin, rbins = radial_profile_calc(nH, T, Z, R_gas,rlow,rhi,nradial)

if(dolin): # if we were doing linear bins, now is the time to switch back to logarithmic.  
    rbins = np.log10(rbins)

# I put a bunch of physical paramters into spherical coordinates.  
P_hot_3D, bin_edges_3D = bs.sphericalmean(P[hot_indexes], coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)
density_hot_3D = bs.sphericalmean(nH[hot_indexes]/XH*M_P, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
vrad_hot_3D = bs.sphericalmean(vradial_hot, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
vtheta_hot_3D = bs.sphericalmean(vtheta_hot, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
vphi_hot_3D = bs.sphericalmean(vphi_hot, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
acc_hot_3D = bs.sphericalmean(accradial_hot, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
nhot_3D = bs.sphericalcount(P[hot_indexes], coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
vrad_mag_hot_3D = bs.sphericalmean(vradial_mag_hot, coords_gas_shift[hot_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
sigrad_hot_3D = np.sqrt(vrad_mag_hot_3D**2 - vrad_hot_3D**2)

print 'bin_edges_3D = ', bin_edges_3D
if(dolin==0):
    bin_edges_R = 10**bin_edges_3D[0]
else:
    bin_edges_R = bin_edges_3D[0]

# This label below uses different gas cuts to do some HSE (hydrostatic equilibrium) calculations, I also have a label to put in the filename    

if(do_hse_calcs):
    print 'bin_edges_R = ', bin_edges_R
    label = "Hot"
    label = "Hot." + str(nradial) + "_" + str(ntheta) + "_" + str(nphi)

# Ezra this calculates and makes hydrostatic equilibrium plots which went into figure 1 of my 2018 paper.  This is rather involved code in hse_spherical.py that I would skip for the time being.  
    hse_spherical.calc_Euler_spherical(bin_edges_R*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_hot_3D, density_hot_3D, vrad_hot_3D*1e+05, vtheta_hot_3D*1e+05, vphi_hot_3D*1e+05, sigrad_hot_3D*1e+05, acc_hot_3D, galidstr, label, M200, R200)
    hse_spherical.calc_Euler_radial(bin_edges_R*cmpermpc, Totmass_radial, P_hot_bin, nH_hot_bin/XH*M_P, vradial_mag_hot_bin*1e+05, vtangential_mag_hot_bin*1e+05, np.sqrt(vradial_mag_hot_bin**2-vradial_hot_bin**2)*1e+05, accradial_hot_bin, galidstr, "Hot." + str(nradial), M200, R200)

# Doing the same as above for cool gas.  
    P_cool_3D, bin_edges_3D = bs.sphericalmean(P[cool_indexes], coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)
    density_cool_3D = bs.sphericalmean(nH[cool_indexes]/XH*M_P, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vrad_cool_3D = bs.sphericalmean(vradial_cool, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vtheta_cool_3D = bs.sphericalmean(vtheta_cool, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vphi_cool_3D = bs.sphericalmean(vphi_cool, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    acc_cool_3D = bs.sphericalmean(accradial_cool, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    ncool_3D = bs.sphericalcount(P[cool_indexes], coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vrad_mag_cool_3D = bs.sphericalmean(vradial_mag_cool, coords_gas_shift[cool_indexes]-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    sigrad_cool_3D = np.sqrt(vrad_mag_cool_3D**2 - vrad_cool_3D**2)

    label = "Cool"
    label = "Cool." + str(nradial) + "_" + str(ntheta) + "_" + str(nphi)
    
    hse_spherical.calc_Euler_spherical(bin_edges_R*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_cool_3D, density_cool_3D, vrad_cool_3D*1e+05, vtheta_cool_3D*1e+05, vphi_cool_3D*1e+05, sigrad_cool_3D*1e+05, acc_cool_3D, galidstr, label, M200, R200)
    fcool_3D = ncool_3D/(nhot_3D+ncool_3D)

# Now doing the same as above for all gas (hot and cool).  
    P_all_3D, bin_edges_3D = bs.sphericalmean(P, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)
    density_all_3D = bs.sphericalmean(nH/XH*M_P, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vrad_all_3D = bs.sphericalmean(vradial, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vtheta_all_3D = bs.sphericalmean(vtheta_all, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    vphi_all_3D = bs.sphericalmean(vphi_all, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    acc_all_3D = bs.sphericalmean(accradial_all, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0] 
    vrad_mag_all_3D = bs.sphericalmean(vradial_mag, coords_gas_shift-boxcenter_pMpc, rlow,rhi,nradial,ntheta,nphi)[0]
    sigrad_all_3D = np.sqrt(vrad_mag_all_3D**2 - vrad_all_3D**2)
    
    label = "All"
    label = "All." + str(nradial) + "_" + str(ntheta) + "_" + str(nphi)
    
    hse_spherical.calc_Euler_spherical(bin_edges_R*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_all_3D, density_all_3D, vrad_all_3D*1e+05, vtheta_all_3D*1e+05, vphi_all_3D*1e+05, sigrad_all_3D*1e+05, acc_all_3D, galidstr, label, M200, R200)
    hse_spherical.calc_Euler_radial(bin_edges_R*cmpermpc, Totmass_radial, P_bin, nH_bin/XH*M_P, vradial_mag_all_bin*1e+05, vtangential_mag_all_bin*1e+05, np.sqrt(vradial_mag_all_bin**2-vradial_all_bin**2)*1e+05, accradial_all_bin, galidstr, "All." + str(nradial), M200, R200)
    
    
# I also did something where I applied a filter to pick out spherical coordianate cells with >90% hot gas.  I'd skip this.  
    label = "Hotfilter"
    label = "Hotfilter." + str(nradial) + "_" + str(ntheta) + "_" + str(nphi)
    
    P_hotfilter_3D = np.where(fcool_3D<0.1, P_hot_3D, -1.0)

    hse_spherical.calc_Euler_spherical(bin_edges_R*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_hotfilter_3D, density_all_3D, vrad_all_3D*1e+05, vtheta_all_3D*1e+05, vphi_all_3D*1e+05, sigrad_all_3D*1e+05, acc_all_3D, galidstr, label, M200, R200)

# Now making cumulative sums.  
gas_mass_bin,rbins = br.radialcumsum(mass_gas,R_gas,rlow,rhi,nradial)
DM_mass_bin = br.radialcumsum(coords_DM[:,0]*0.0+mass_table[1]*1e+10/hubbleparam,R_DM,rlow,rhi,nradial)[0]
stars_mass_bin = br.radialcumsum(mass_stars,R_stars,rlow,rhi,nradial)[0]
hot_mass_bin = br.radialcumsum(mass_gas[hot_indexes],R_gas[hot_indexes],rlow,rhi,nradial)[0]
cool_mass_bin = br.radialcumsum(mass_gas[cool_indexes],R_gas[cool_indexes],rlow,rhi,nradial)[0]
cool_nosatISM_mass_bin = br.radialcumsum(mass_gas[cool_nosatISM_indexes],R_gas[cool_nosatISM_indexes],rlow,rhi,nradial)[0]
cool_noISM_mass_bin = br.radialcumsum(mass_gas[cool_noISM_indexes],R_gas[cool_noISM_indexes],rlow,rhi,nradial)[0]

print "Finished cumulative sums"

acc_hot_bin = P_hot_bin*0.0 

# I write files of radial profiles for all the physcal parameters above (for hot, cool, and all gas).  
fradial = open("rbins_%s.dat"%galidstr,"w")
fradial_hot = open("rbins_hot_%s.dat"%galidstr,"w")
fradial_cool = open("rbins_cool_%s.dat"%galidstr,"w")

for i in xrange(len(gas_mass_bin)):
    print rbins[i], rbins[i+1],gas_mass_bin[i],T_hot_bin[i],nH_hot_bin[i],P_hot_bin[i],S_hot_bin[i]
    fradial.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f %5.3e\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],gas_mass_bin[i],T_bin[i],nH_bin[i],P_bin[i],S_bin[i],vradial_all_bin[i],vradial_mag_all_bin[i],vtangential_all_bin[i],vtangential_mag_all_bin[i],Z_bin[i]))
    fradial_hot.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f %5.3e\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],hot_mass_bin[i],T_hot_bin[i],nH_hot_bin[i],P_hot_bin[i],S_hot_bin[i],vradial_hot_bin[i],vradial_mag_hot_bin[i],vtangential_hot_bin[i],vtangential_mag_hot_bin[i],Z_hot_bin[i]))
    fradial_cool.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f %5.3e  %5.3e\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],cool_mass_bin[i],T_cool_bin[i],nH_cool_bin[i],P_cool_bin[i],S_cool_bin[i],vradial_cool_bin[i],vradial_mag_cool_bin[i],vtangential_cool_bin[i],vtangential_mag_cool_bin[i], Z_cool_bin[i], cool_noISM_mass_bin[i]))

fradial.close()
fradial_hot.close()
fradial_cool.close()

print "Finished basic radial profiles"

# Function calls to above vels_radial_plot to output velocity files and also to make plots.  
vels_radial_plot(rbins, vradial_mag_hot_bin, vradial_hot_bin, vtangential_mag_hot_bin, vtangential_hot_bin, galidstr, "vel_hot", M200, R200, v200)
vels_radial_plot(rbins, vradial_mag_cool_bin, vradial_cool_bin, vtangential_mag_cool_bin, vtangential_cool_bin, galidstr, "vel_cool", M200, R200, v200)

### DM velocities- dark matter has velocity info as well.  I also calculate DM angular momentum too.  
vradial_DM_mag, vtangential_DM_mag, vradial_DM, vtangential_DM = vels_calc(coords_DM_shift,vels_DM_shift,boxcenter_pMpc,[0,0,0])
print "v_DM_means (<R>,<tan>,R,tan) ", np.mean(vradial_DM_mag), np.mean(vtangential_DM_mag), np.mean(vradial_DM), np.mean(vtangential_DM)

vradial_DM_bin = br.radialphysical(vradial_DM,R_DM,rlow,rhi,nradial)[0]
vradial_DM_mag_bin = br.radialphysical(vradial_DM_mag,R_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_x = br.radialphysical(vtangential_DM[:,0],R_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_y = br.radialphysical(vtangential_DM[:,1],R_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_z = br.radialphysical(vtangential_DM[:,2],R_DM,rlow,rhi,nradial)[0]
vtangential_DM_mag_bin = br.radialphysical(vtangential_DM_mag,R_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin = np.sqrt(vtangential_DM_bin_x**2+vtangential_DM_bin_y**2+vtangential_DM_bin_z**2)

vels_radial_plot(rbins, vradial_DM_mag_bin, vradial_DM_bin, vtangential_DM_mag_bin, vtangential_DM_bin, galidstr, "vel_DM", M200, R200, v200)

vradial_stars_mag, vtangential_stars_mag, vradial_stars, vtangential_stars = vels_calc(coords_stars_shift,vels_stars_shift,boxcenter_pMpc,[0,0,0])
print "LENGTHS VECTORS: ", len(vtangential_stars_mag), len(R_stars), len(mass_stars), len(R_stars)

angmom_DM = np.cross(vels_DM_shift,(coords_DM_shift-boxcenter_pMpc)*1e+03)
angmom_DM_bin_x = br.radialcumsum(angmom_DM[:,0]*(R_DM*0.0+mass_table[1]*1e+10/hubbleparam),R_DM,rlow,rhi,nradial)[0]
angmom_DM_bin_y = br.radialcumsum(angmom_DM[:,1]*(R_DM*0.0+mass_table[1]*1e+10/hubbleparam),R_DM,rlow,rhi,nradial)[0]
angmom_DM_bin_z = br.radialcumsum(angmom_DM[:,2]*(R_DM*0.0+mass_table[1]*1e+10/hubbleparam),R_DM,rlow,rhi,nradial)[0]
angmom_DM_bin = np.sqrt(angmom_DM_bin_x**2+angmom_DM_bin_y**2+angmom_DM_bin_z**2)
angmom_DM_bin /= DM_mass_bin
print "ANGMOM_DM = ", angmom_DM_bin

angmom_stars = np.cross(vels_stars_shift,(coords_stars_shift-boxcenter_pMpc)*1e+03)
angmom_stars_bin_x = br.radialcumsum(angmom_stars[:,0]*mass_stars,R_stars,rlow,rhi,nradial)[0]
angmom_stars_bin_y = br.radialcumsum(angmom_stars[:,1]*mass_stars,R_stars,rlow,rhi,nradial)[0]
angmom_stars_bin_z = br.radialcumsum(angmom_stars[:,2]*mass_stars,R_stars,rlow,rhi,nradial)[0]
angmom_stars_bin = np.sqrt(angmom_stars_bin_x**2+angmom_stars_bin_y**2+angmom_stars_bin_z**2)
angmom_stars_bin /= stars_mass_bin

angmom_gas = np.cross(vels_gas_shift,(coords_gas_shift-boxcenter_pMpc)*1e+03)
angmom_gas_bin_x = br.radialcumsum(angmom_gas[:,0]*mass_gas,R_gas,rlow,rhi,nradial)[0]
angmom_gas_bin_y = br.radialcumsum(angmom_gas[:,1]*mass_gas,R_gas,rlow,rhi,nradial)[0]
angmom_gas_bin_z = br.radialcumsum(angmom_gas[:,2]*mass_gas,R_gas,rlow,rhi,nradial)[0]
angmom_gas_bin = np.sqrt(angmom_gas_bin_x**2+angmom_gas_bin_y**2+angmom_gas_bin_z**2)
angmom_gas_bin /= gas_mass_bin

angmom_hot = np.cross(vels_gas_shift[hot_indexes],(coords_gas_shift[hot_indexes]-boxcenter_pMpc)*1e+03)
angmom_hot_bin_x = br.radialcumsum(angmom_hot[:,0]*mass_gas[hot_indexes],R_gas[hot_indexes],rlow,rhi,nradial)[0]
angmom_hot_bin_y = br.radialcumsum(angmom_hot[:,1]*mass_gas[hot_indexes],R_gas[hot_indexes],rlow,rhi,nradial)[0]
angmom_hot_bin_z = br.radialcumsum(angmom_hot[:,2]*mass_gas[hot_indexes],R_gas[hot_indexes],rlow,rhi,nradial)[0]
angmom_hot_bin = np.sqrt(angmom_hot_bin_x**2+angmom_hot_bin_y**2+angmom_hot_bin_z**2)
angmom_hot_bin /= hot_mass_bin

angmom_cool = np.cross(vels_gas_shift[cool_indexes],(coords_gas_shift[cool_indexes]-boxcenter_pMpc)*1e+03)
angmom_cool_bin_x = br.radialcumsum(angmom_cool[:,0]*mass_gas[cool_indexes],R_gas[cool_indexes],rlow,rhi,nradial)[0]
angmom_cool_bin_y = br.radialcumsum(angmom_cool[:,1]*mass_gas[cool_indexes],R_gas[cool_indexes],rlow,rhi,nradial)[0]
angmom_cool_bin_z = br.radialcumsum(angmom_cool[:,2]*mass_gas[cool_indexes],R_gas[cool_indexes],rlow,rhi,nradial)[0]
angmom_cool_bin = np.sqrt(angmom_cool_bin_x**2+angmom_cool_bin_y**2+angmom_cool_bin_z**2)
angmom_cool_bin /= cool_mass_bin

angmom_cool_nosatISM = np.cross(vels_gas_shift[cool_nosatISM_indexes],(coords_gas_shift[cool_nosatISM_indexes]-boxcenter_pMpc)*1e+03)
angmom_cool_nosatISM_bin_x = br.radialcumsum(angmom_cool_nosatISM[:,0]*mass_gas[cool_nosatISM_indexes],R_gas[cool_nosatISM_indexes],rlow,rhi,nradial)[0]
angmom_cool_nosatISM_bin_y = br.radialcumsum(angmom_cool_nosatISM[:,1]*mass_gas[cool_nosatISM_indexes],R_gas[cool_nosatISM_indexes],rlow,rhi,nradial)[0]
angmom_cool_nosatISM_bin_z = br.radialcumsum(angmom_cool_nosatISM[:,2]*mass_gas[cool_nosatISM_indexes],R_gas[cool_nosatISM_indexes],rlow,rhi,nradial)[0]
angmom_cool_nosatISM_bin = np.sqrt(angmom_cool_nosatISM_bin_x**2+angmom_cool_nosatISM_bin_y**2+angmom_cool_nosatISM_bin_z**2)
angmom_cool_nosatISM_bin /= cool_nosatISM_mass_bin



##angmom_DM_bin = br.radialcumsum(vtangential_DM*(R_DM*1e+03)*(R_DM*0.0+mass_table[1]*1e+10/hubbleparam),R_DM,rlow,rhi,nradial)[0]
#angmom_gas_bin = br.radialcumsum(vtangential_mag*(R_gas*1e+03)*(mass_gas),R_gas,rlow,rhi,nradial)[0]
#angmom_hot_bin = br.radialcumsum(vtangential_mag[hot_indexes]*(R_gas[hot_indexes]*1e+03)*(mass_gas[hot_indexes]),R_gas[hot_indexes],rlow,rhi,nradial)[0]
#angmom_cool_bin = br.radialcumsum(vtangential_mag[cool_indexes]*(R_gas[cool_indexes]*1e+03)*(mass_gas[cool_indexes]),R_gas[cool_indexes],rlow,rhi,nradial)[0]

#if(dolin):
#    rbins_cgs = rbins*cmpermpc
#else:
rbins_cgs = 10**rbins*cmpermpc

# rbins_M is sihifted by one index because I believe the rbins_m listed encompasses everything inside.  
rbins_M = np.zeros(len(angmom_DM_bin))
for i in xrange(len(rbins_M)-1):
    rbins_M[i] = rbins_cgs[i+1]

#This ia a messy file that I write out angular momentum in different quantities.  I wouldn't worry about this for now.  In fact I should write a separate function for this.  
fangmom = open("angmom_%s.dat"%galidstr,"w")
fangmom.write("#1: rbins 2: DM_lambda 3: gas_lambda 4: star_lambda 5: cool_lambda 6:hot_lambda 7: DM_sp_angmom 8: gas_sp_angmom 9: star_sp_angmom 10: cool_sp_angmom 11: hot_sp_angmom 12: DM_angmom 13: gas_angmom 14: star_angmom 15: cool_angmom 16: hot_angmom 17: coolnosatISM_lambda 18: coonosatISM__sp__angmom\n")
for i in xrange(len(angmom_gas_bin)):    
    fangmom.write("% 5.3f %6.4f %6.4f %6.4f %6.4f %6.4f   %5.3e %5.3e %5.3e %5.3e %5.3e   %5.3e %5.3e %5.3e %5.3e %5.3e  %6.4f %5.3e\n"%(np.log10(rbins_M[i]/cmpermpc),angmom_DM_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_gas_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_stars_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_cool_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_hot_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_DM_bin[i],angmom_gas_bin[i],angmom_stars_bin[i],angmom_cool_bin[i],angmom_hot_bin[i],angmom_DM_bin[i]*DM_mass_bin[i],angmom_gas_bin[i]*gas_mass_bin[i],angmom_stars_bin[i]*stars_mass_bin[i],angmom_cool_bin[i]*cool_mass_bin[i],angmom_hot_bin[i]*hot_mass_bin[i],angmom_cool_nosatISM_bin[i]/(np.sqrt(2)*(R200*1e+03)*v200),angmom_cool_nosatISM_bin[i]))

fangmom.close()

print "Finished velocity radial profiles"

### Plot mass as a function of radius and write it out.  Probably redundant with above.  
massfrac_radial_plot(rbins,DM_mass_bin,gas_mass_bin,stars_mass_bin,cool_mass_bin,hot_mass_bin,cool_noISM_mass_bin,galidstr,M200,R200)

### Energy densities (of kinetic motions and thermal) cumulatively binned.  
mu = 0.59
Ekinetic_cum_bin = br.radialcumsum(0.5*(vradial_mag[hot_indexes]**2+vtangential_mag[hot_indexes]**2)*1e+10*mass_gas[hot_indexes]*M_Solar, R_gas[hot_indexes],rlow,rhi,nradial)[0]
Ethermal_cum_bin = br.radialcumsum(3/2.*K_Boltz*T[hot_indexes]/M_P/mu*(mass_gas[hot_indexes]*M_Solar), R_gas[hot_indexes],rlow,rhi,nradial)[0]
Ethermal_cum_bin = br.radialcumsum(T[hot_indexes]*mass_gas[hot_indexes], R_gas[hot_indexes],rlow,rhi,nradial)[0]
Ethermal_cum_bin = Ethermal_cum_bin*3/2.*(K_Boltz)/M_P/mu*M_Solar

#Egravity_cum_bin = br.radialcumsum(3/5.*G_Grav*tot R_gas[hot_indexes],rlow,rhi,nradial)[0]
#Egravity_cum_bin = 3/5.*G_Grav*(DM_mass_bin+stars_mass_bin+gas_mass_bin)*M_Solar/rbins_M

# energy plot as listed above.  Gravitational energy is calculated from sum of DM stars and gas that is binned.  
energy_radial_plot(rbins,Ekinetic_cum_bin,Ethermal_cum_bin,(DM_mass_bin+stars_mass_bin+gas_mass_bin)*M_Solar,gas_mass_bin*M_Solar,galidstr,M200,R200)

print "Finished energy radial profiles"

# Obselete function call here, so I'm going to comment out 
###Mtot, Mtherm = calc_hydrostatic_masses(rbins, DM_mass_bin+stars_mass_bin+gas_mass_bin, P_hot_bin, nH_hot_bin, vradial_hot_bin, vradial_mag_hot_bin, vtangential_hot_bin, vtangential_mag_hot_bin,accradial_hot_bin,galidstr,M200,R200)

print "Finished hydrostatic radial profiles"


#exit()

# This is to do figures of column density maps and other maps on a 2-D projection using coldens library- but don't bother, dofigs=0.

lgrid = 1.5*R200 #3*R200
#lgrid = 600/1e+03
lgridz = lgrid*4.
ngrid = 600

acc200 = (v200*1e+05)**2./(R200*cmpermpc)

if dofigs:
    for k in xrange(3):
        if(k == 0):
            theta = 90.
            phi = 0.
            psi = 90.
            dir = 'x'
        if(k == 1):
            theta = 90.
            phi = 0.
            psi = 0.
            dir = 'y'
        if(k == 2):
            theta = 0.
            phi = 0.
            psi = 0.
            dir = 'z'

            
        vels_gas_dir = vels_gas_rot[:,k] - gal_vels_kms[k]
        sigsq_gas_dir = (vels_gas_shift[:,k])**2 #(vels_gas_rot[:,k] - gal_vels_kms[k])**2.
        print "vels_gas_dir = ", vels_gas_dir
        print "sigsq_gas_dir = ", sigsq_gas_dir
        
        result_h = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.hydrogen.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=19.0, Vmax=21.5,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)
        
        result_hvel = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*vels_gas_dir, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.vel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)
        print "result_h = ", result_h
        ###print "result_hvel = ", result_hvel
        ###print "result_hvel/10**result_h", result_hvel/10**result_h
    
        coldens.make_colourmap('coldens.%s.vel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hvel/10**result_h, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

        if(acceleration):
            accel_gas_dir = accel_gas_rot[:,k]
            result_hacc = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*accel_gas_dir, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.acc.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-acc200*1e+19, Vmax=acc200*1e+19,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='Hacc',haloinfostr=haloinfostr,docbar=1)

            coldens.make_colourmap('coldens.%s.acc.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hacc/10**result_h, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -acc200, acc200, 'Acceleration', ngrid, redshift=redshift, extralabel='Acceleration',haloinfostr=haloinfostr,docbar=1)

            
        result_hsig = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*sigsq_gas_dir, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.sig.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=0, Vmax=v200,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='Hsig',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.sig.%s.l%3.1f.png'%(galidstr,dir,lgrid), np.sqrt(10**(result_hsig-result_h)), -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, 0, v200, 'Dispersion', ngrid, redshift=redshift, extralabel='Vel Disp.',haloinfostr=haloinfostr,docbar=1)
    
        if(chem):
            result_h1 = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*h1, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.h1.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=13.0, Vmax=21.0,ion='HI',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='HI',haloinfostr=haloinfostr,docbar=1)

        result_hcool = coldens.main(coords_gas_rot, hsmooth, cool_gas*hydrogen, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.hcool.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=17.5, Vmax=21.5,ion='H-cool',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)

        result_hhot = coldens.main(coords_gas_rot, hsmooth, hot_gas*hydrogen, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.hhot.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=17.5, Vmax=21.5,ion='H-hot',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)

        result_hcoolvel = coldens.main(coords_gas_rot, hsmooth, cool_gas*hydrogen*vels_gas_dir, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.coolvel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.coolvel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hcoolvel/10**result_hcool, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Cool Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

        result_hhotvel = coldens.main(coords_gas_rot, hsmooth, hot_gas*hydrogen*vels_gas_dir, gal_coords_pMpc, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_pMpc,fig_name='coldens.%s.hotvel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,psi=psi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.hotvel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hhotvel/10**result_hhot, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Hot Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

    
