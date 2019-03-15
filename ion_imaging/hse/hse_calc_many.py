# We need to move the original output temperatures to a
# new array, create a new array with the corrected
# temperatures and remove any equilibrium chemistry
# abundance arrays.
import sys
import tables
import eagle
import h5py
import glob as glob
import os
import numpy as np
import scipy.stats as scist
import galmanip.readGalaxy as rg
import galmanip.rotateGalaxy as rotg
import galmanip.writeGalaxy as wg
import galmanip.binRadial as br
import galmanip.binSpherical as bs
import coldens_ben.coldens as coldens
from astropy import constants as const
import matplotlib.pyplot as plt
import hse_spherical 

G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24

def vels_calc(coords,vels,coords_center,vels_center):
    #r = np.sqrt((coords[:,0]-coords_center[0])**2+(coords[:,1]-coords_center[1])**2+(coords[:,2]-center_coords[2])**2)
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

def massfrac_radial_plot(rbins,DM_mass_bin,gas_mass_bin,stars_mass_bin,cool_mass_bin,hot_mass_bin,galidstr,M200,R200):

    fig_mfrac = plt.figure(figsize=(5.0,4.0))
    ax_mfrac = fig_mfrac.add_subplot(111)

    rbins_cgs = 10**rbins*cmpermpc
    rbins_M = np.zeros(len(DM_mass_bin))
    
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]

    fmass = open("mass_%s.dat"%galidstr,"w")

    fmass.write("#rbins DMmass gasmass starmass coolmass hotmass\n")

    for i in xrange(len(rbins_M)):    
        fmass.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e %5.3e\n"%(np.log10(rbins_M[i]/cmpermpc),DM_mass_bin[i],gas_mass_bin[i],stars_mass_bin[i],cool_mass_bin[i],hot_mass_bin[i]))

    fmass.close()



    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(DM_mass_bin+gas_mass_bin+stars_mass_bin), color = 'black', lw=4,label="$M_{\mathrm{tot}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(stars_mass_bin), color = 'blue', lw=4,label="$M_{\mathrm{stars}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(gas_mass_bin), color = 'red', lw=4,label="$M_{\mathrm{gas}}$")

    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(gas_mass_bin*0.307/0.04825), color = 'red', lw=2, ls=":",label="$M_{\mathrm{gas}} \Omega_{\mathrm{M}}/\Omega_{\mathrm{b}}$")

    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(hot_mass_bin), color = 'orange', lw=2,label="$M_{\mathrm{hot}}$")
    ax_mfrac.plot(np.log10(rbins_M/cmpermpc), np.log10(cool_mass_bin), color = 'cyan', lw=2,label="$M_{\mathrm{cool}}$")

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

    rbins_cgs = 10**rbins*cmpermpc
    rbins_M = np.zeros(len(vrad_bin))
    
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]

    fvel = open("%s_%s.dat"%(prefix,galidstr),"w")

    fvel.write("#rbins sigrad vrad sigtan vtan\n")

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
            

def radial_profile_calc(nH, T, logR, Rloglow, Rloghi, nbins):


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

    return(T_bin, nH_bin, P_bin, S_bin, bins)

#def energy_radial_plot(rbins,vradial_mag_bin,vtangential_mag_bin,T_bin,totmass_bin,M200,R200):
def energy_radial_plot(rbins,Ekinetic,Ethermal,totmass_bin,gasmass_bin,galidstr,M200,R200):


    mu = 0.59

    rbins_cgs = 10**rbins*cmpermpc
    rbins_M = np.zeros(len(Ekinetic))
    for i in xrange(len(rbins_M)-1):
        rbins_M[i] = rbins_cgs[i+1]


    print "LENGTHS: rbins_M= ", len(rbins_M), "totmass, gasmass ", totmass_bin, gasmass_bin
    #Ekinetic = 0.5*(vradial_mag_bin**2 + vtangential_mag_bin**2)*1e+10
    #Ethermal = 3/2.*K_Boltz*T_bin/M_P/mu # ADD MU
    Egravity = 3/5.*G_Grav*(totmass_bin*gasmass_bin)/rbins_M
    Vir_ratio = (Ekinetic+Ethermal)/Egravity
    kin_ther_ratio = Ekinetic/Ethermal
    print "Ekinetic = ", np.log10(Ekinetic)
    print "Ethermal = ", np.log10(Ethermal)
    print "Egravity = ", np.log10(Egravity)

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
    rbins_cgs = 10**rbins*cmpermpc
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

##def whole_sim_read()



chem = 0
dofigs = 0
cosmoowls = 0
###acceleration = 1

runlabel = 'halo'

unit_mass_in_cgs = 1.989e33 * 1.0e10 
unit_length_in_cgs = 3.0857e24 
proton_mass_cgs = 1.67e-24

unit_Density_in_cgs = unit_mass_in_cgs/unit_length_in_cgs**3/proton_mass_cgs 




if (runlabel=='halo'):
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
sniptag = "SNAP"
sim = "."

#h5file_inbase = "%ssnapshot_%s/snap_%s" % (data_dir, snap_intag, snap_intag)

hubbleparam,aex,redshift,boxsize,mass_table = rg.header(sniptag,sim,snap_intag)

omegaM = 0.307
omegaL = 1-omegaM
omegaratio = (omegaM+omegaL/(1+redshift)**3)
R200 = 1.63e-5*(10**float(M200)*hubbleparam)**0.333/omegaratio**0.333/(1+redshift)/hubbleparam
v200 = np.sqrt(G_Grav*10**M200*M_Solar/(R200*cmpermpc))/1e+05

boxsize_proper = boxsize*aex/hubbleparam
print "boxsize_proper = ", boxsize_proper

gal_coords = ((float(x_gal*aex/hubbleparam),float(y_gal*aex/hubbleparam),float(z_gal*aex/hubbleparam)))
gal_vels = ((float(v_x_gal),float(v_y_gal),float(v_z_gal)))

box_center = ((float(boxsize_proper/2.),float(boxsize_proper/2.),float(boxsize_proper/2.)))

coords_DM,vels_DM = rg.coords(sniptag,sim,snap_intag,1)
#mass_DM = rg.mass(sniptag,sim,snap_intag,1)


coords_stars,vels_stars = rg.coords(sniptag,sim,snap_intag,4)
mass_stars = rg.mass(sniptag,sim,snap_intag,4)*1e+10

coords_gas,vels_gas = rg.coords(sniptag,sim,snap_intag,0)
mass_gas = rg.mass(sniptag,sim,snap_intag,0)*1e+10


if(cosmoowls):
    nH, T = rg.gasprops_cosmoowls(sniptag,sim,snap_intag)
else:
    nH, T, Z, SFR, hsmooth = rg.gasprops(sniptag,sim,snap_intag)



datadir, snapshot, x_gal, y_cal, z_gal, galidstr, halodir, v_x_gal, v_y_gal, v_z_gal = np.loadtxt(sys.argv[1], usecols=(1,2,3,4,5,6,7,11,12,13))

sim = "/net/virgo/data5/oppenheimer/" + halodir[0] + "/" + datadir[0]



snap_intag = sys.argv[1]
galidstr = sys.argv[2]
x_gal = float(sys.argv[3])
y_gal = float(sys.argv[4])
z_gal = float(sys.argv[5])
v_x_gal = float(sys.argv[6])
v_y_gal = float(sys.argv[7])
v_z_gal = float(sys.argv[8])


    
#print "coords_stars = ", coords_stars
angmom_indexes = rg.distanceindexes(coords_stars,gal_coords,0.03/hubbleparam,boxsize_proper)
print "angmom_array = ", coords_stars[angmom_indexes], "length = ", len(coords_stars[angmom_indexes])
#print "mean x,y,z = ", np.mean(coords_stars[angmom_indexes,0]),np.mean(coords_stars[angmom_indexes,1]),np.mean(coords_stars[angmom_indexes,2])

phi, theta, psi = rotg.calc_angmom_axisangles(gal_coords, gal_vels, coords_stars[angmom_indexes], vels_stars[angmom_indexes], boxsize_proper)
print "Finished ang mom calculation "

gas_indexes = rg.distanceindexes(coords_gas,gal_coords,4.0,boxsize_proper)
print "Finished gas indexing, gas_indexes length = ", len(coords_gas[gas_indexes]) 
DM_indexes = rg.distanceindexes(coords_DM,gal_coords,4.0,boxsize_proper)
print "Finished DM indexing, DM_indexes length = ", len(coords_DM[DM_indexes]) 
stars_indexes = rg.distanceindexes(coords_stars,gal_coords,4.0,boxsize_proper)
print "Finished stellar indexing, star_indexes length = ", len(coords_stars[stars_indexes]) 

coords_gas_rot, vels_gas_rot = rotg.rotategalaxy(coords_gas[gas_indexes],vels_gas[gas_indexes],gal_coords,gal_vels,boxsize_proper,phi,theta,psi)
coords_stars_rot, vels_stars_rot = rotg.rotategalaxy(coords_stars[stars_indexes],vels_stars[stars_indexes],gal_coords,gal_vels,boxsize_proper,phi,theta,psi)
coords_DM_rot, vels_DM_rot = rotg.rotategalaxy(coords_DM[DM_indexes],vels_DM[DM_indexes],gal_coords,gal_vels,boxsize_proper,phi,theta,psi)

print "Finished rotating galaxy"

coords_gas_shift, vels_gas_shift = rotg.shiftgalaxy(coords_gas_rot,vels_gas_rot,gal_coords,gal_vels,box_center,[0,0,0],boxsize_proper)
coords_DM_shift, vels_DM_shift = rotg.shiftgalaxy(coords_DM_rot,vels_DM_rot,gal_coords,gal_vels,box_center,[0,0,0],boxsize_proper)
coords_stars_shift, vels_stars_shift = rotg.shiftgalaxy(coords_stars_rot,vels_stars_rot,gal_coords,gal_vels,box_center,[0,0,0],boxsize_proper)

print "Finished shifting galaxy"

print "coords_stars_rot= ", coords_stars_rot
angmom_indexes = rg.distanceindexes(coords_stars_rot,gal_coords,0.03/hubbleparam,boxsize_proper)
print "angmom_indexes = ", coords_stars_rot[angmom_indexes], "length = ", len(coords_stars_rot[angmom_indexes])
print "mean x,y,z = ", np.mean(coords_stars_rot[angmom_indexes,0]),np.mean(coords_stars_rot[angmom_indexes,1]),np.mean(coords_stars_rot[angmom_indexes,2])

### CHECK ON ANG MOM AXIS
phi, theta, psi = rotg.calc_angmom_axisangles(gal_coords, gal_vels, coords_stars_rot[angmom_indexes], vels_stars_rot[angmom_indexes], boxsize_proper)

if (cosmoowls):
    wg.write_gas("snap_%s.0.hdf5"%galidstr, coords_gas_shift,vels_gas_shift,mass_gas[gas_indexes],nH[gas_indexes]/0.752/unit_Density_in_cgs,T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],T[gas_indexes],aex,redshift,boxsize,hubbleparam,mass_table,cosmoowls)
    #print "hsmooth_read = ", hsmooth, len(hsmooth)
    hsmooth = (nH[gas_indexes]/0.752/unit_Density_in_cgs/mass_gas[gas_indexes])**(-0.3333)*0.00111066
    print "hsmooth_calc = ", hsmooth, len(hsmooth)
    hydrogen = nH[gas_indexes]*0.0 + 0.752
    SFR = nH[gas_indexes]*0.0 
else:
    hydrogen,carbon,oxygen,magnesium = rg.metalsprimary(sniptag, sim, snap_intag)
    helium,nitrogen,neon,silicon,iron = rg.metalssecondary(sniptag, sim, snap_intag)

    wg.write_gas("snap_%s.0.hdf5"%galidstr, coords_gas_shift,vels_gas_shift,mass_gas[gas_indexes],nH[gas_indexes]/hydrogen[gas_indexes]/unit_Density_in_cgs,T[gas_indexes],Z[gas_indexes],SFR[gas_indexes],hsmooth[gas_indexes],hydrogen[gas_indexes],helium[gas_indexes],carbon[gas_indexes],nitrogen[gas_indexes],oxygen[gas_indexes],neon[gas_indexes],magnesium[gas_indexes],silicon[gas_indexes],iron[gas_indexes],aex,redshift,boxsize,hubbleparam,mass_table,cosmoowls)


    hsmooth = hsmooth[gas_indexes]
    hydrogen = hydrogen[gas_indexes]

    oxygen = oxygen[gas_indexes]

wg.write_stars("snap_%s.0.hdf5"%galidstr,coords_stars_shift,vels_stars_shift,mass_stars[stars_indexes],aex,redshift,boxsize,hubbleparam)
wg.write_DM("snap_%s.0.hdf5"%galidstr,coords_DM_shift,vels_DM_shift,aex,redshift,boxsize,hubbleparam)

if (chem):
    h1,c4,o6,o7,o8,mg2 = rg.chemistryprimary(sniptag, sim, snap_intag)
    h1 = h1[gas_indexes]
    o6 = o6[gas_indexes]

mass_gas = mass_gas[gas_indexes]
nH = nH[gas_indexes]
T = T[gas_indexes]
P = T*nH/(0.752*0.59)
cool_gas = np.where(T<1e+05, mass_gas, 0.0)
hot_gas = np.where(T>=1e+05, mass_gas, 0.0)

hot_indexes = np.where(T>=1e+05)
cool_indexes = np.where(T<1e+05)
lt104_indexes = np.where(T<1e+04)
gt106_indexes = np.where(T>=1e+06)

R_gas = np.sqrt((coords_gas_shift[:,0]-box_center[0])**2+(coords_gas_shift[:,1]-box_center[1])**2+(coords_gas_shift[:,2]-box_center[2])**2)
logR_gas = np.log10(R_gas)
logR_gas_hot = logR_gas[hot_indexes]
logR_gas_cool = logR_gas[cool_indexes]
print "logR_gas = ", logR_gas   
R_DM = np.sqrt((coords_DM_shift[:,0]-box_center[0])**2+(coords_DM_shift[:,1]-box_center[1])**2+(coords_DM_shift[:,2]-box_center[2])**2)
logR_DM = np.log10(R_DM)   
print "logR_DM = ", logR_DM
R_stars = np.sqrt((coords_stars_shift[:,0]-box_center[0])**2+(coords_stars_shift[:,1]-box_center[1])**2+(coords_stars_shift[:,2]-box_center[2])**2)
logR_stars = np.log10(R_stars)   
print "logR_stars = ", logR_stars

M_stars = np.sum(mass_stars[np.where(R_stars<R200)])
M_stars_30kpc = np.sum(mass_stars[np.where(R_stars<0.03)])

M_DM =  len(R_DM[np.where(R_DM<R200)])*mass_table[1]*1e+10/hubbleparam
M_CGM = np.sum(mass_gas[np.where((R_gas<R200) & (SFR[gas_indexes]<=0))])
M_ISM = np.sum(mass_gas[np.where((R_gas<R200) & (SFR[gas_indexes]>0))])
M_ISM_30kpc = np.sum(mass_gas[np.where((R_gas<0.03) & (SFR[gas_indexes]>0))])
print "MASSES= ", np.log10(M_DM), np.log10(M_CGM), " ISMs= ", np.log10(M_ISM), np.log10(M_ISM_30kpc), " Stars= ", np.log10(M_stars),ms,np.log10(M_stars_30kpc)
M_lt104 = np.sum(mass_gas[np.where((R_gas<R200) & (T<1e+04))])
M_104_105 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+04) & (T<1e+05))])
M_105_106 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+05) & (T<1e+06))])
M_gt106 = np.sum(mass_gas[np.where((R_gas<R200) & (T>=1e+06))])

print "MASSES_gas= ", np.log10(M_lt104), np.log10(M_104_105), np.log10(M_105_106), np.log10(M_gt106)

fmasses = open("masses_%s.dat"%galidstr,"w")
fmasses.write("galidstr= %5.2f %5.2f %5.3f MDM= %5.3e Mstars=  %5.3e %5.3e MISMs= %5.3e %5.3e MCGM= %5.3e MCGM_sub= %5.3e %5.3e %5.3e %5.3e\n"%(M200,ms,sfr,M_DM,M_stars,M_stars_30kpc,M_ISM,M_ISM_30kpc,M_CGM,M_lt104,M_104_105,M_105_106,M_gt106))
fmasses.close()



rlow = np.log10(R200/316.)
rhi = np.log10(R200*3.16)
nradial= 60 #120
ntheta = 18 #36
nphi = 10 #20

vradial_mag, vtangential_mag, vradial, vtangential = vels_calc(coords_gas_shift,vels_gas_shift,box_center,[0,0,0])
print "v_means (<R>,<tan>,R,tan) ", np.mean(vradial_mag), np.mean(vtangential_mag), np.mean(vradial), np.mean(vtangential)

vradial_mag_hot, vtangential_mag_hot, vradial_hot, vtangential_hot = vels_calc(coords_gas_shift[hot_indexes],vels_gas_shift[hot_indexes],box_center,[0,0,0])
vradial_mag_cool, vtangential_mag_cool, vradial_cool, vtangential_cool = vels_calc(coords_gas_shift[cool_indexes],vels_gas_shift[cool_indexes],box_center,[0,0,0])

#vradial_bin = br.radialphysical(vradial,logR_gas,rlow,rhi,nradial)[0]
#vradial_mag_bin = br.radialphysical(vradial_mag,logR_gas,rlow,rhi,nradial)[0]
#vtangential_bin_x = br.radialphysical(vtangential[:,0],logR_gas,rlow,rhi,nradial)[0]
#vtangential_bin_y = br.radialphysical(vtangential[:,1],logR_gas,rlow,rhi,nradial)[0]
#vtangential_bin_z = br.radialphysical(vtangential[:,2],logR_gas,rlow,rhi,nradial)[0]
#vtangential_mag_bin = br.radialphysical(vtangential_mag,logR_gas,rlow,rhi,nradial)[0]
#vtangential_bin = np.sqrt(vtangential_bin_x**2+vtangential_bin_y**2+vtangential_bin_z**2)

vradial_all_bin = br.radialphysical(vradial_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vradial_mag_all_bin = br.radialphysical(vradial_mag_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_all_bin_x = br.radialphysical(vtangential_hot[:,0],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_all_bin_y = br.radialphysical(vtangential_hot[:,1],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_all_bin_z = br.radialphysical(vtangential_hot[:,2],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_mag_all_bin = br.radialphysical(vtangential_mag_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_all_bin = np.sqrt(vtangential_all_bin_x**2+vtangential_all_bin_y**2+vtangential_all_bin_z**2)

vradial_hot_bin = br.radialphysical(vradial_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vradial_mag_hot_bin = br.radialphysical(vradial_mag_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_x = br.radialphysical(vtangential_hot[:,0],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_y = br.radialphysical(vtangential_hot[:,1],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin_z = br.radialphysical(vtangential_hot[:,2],logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_mag_hot_bin = br.radialphysical(vtangential_mag_hot,logR_gas_hot,rlow,rhi,nradial)[0]
vtangential_hot_bin = np.sqrt(vtangential_hot_bin_x**2+vtangential_hot_bin_y**2+vtangential_hot_bin_z**2)

vradial_cool_bin = br.radialphysical(vradial_cool,logR_gas_cool,rlow,rhi,nradial)[0]
vradial_mag_cool_bin = br.radialphysical(vradial_mag_cool,logR_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_x = br.radialphysical(vtangential_cool[:,0],logR_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_y = br.radialphysical(vtangential_cool[:,1],logR_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin_z = br.radialphysical(vtangential_cool[:,2],logR_gas_cool,rlow,rhi,nradial)[0]
vtangential_mag_cool_bin = br.radialphysical(vtangential_mag_cool,logR_gas_cool,rlow,rhi,nradial)[0]
vtangential_cool_bin = np.sqrt(vtangential_cool_bin_x**2+vtangential_cool_bin_y**2+vtangential_cool_bin_z**2)

#print "vradial_mag_hot_bin= ", vradial_mag_hot_bin
#print "vtangential_mag_hot_bin= ", vtangential_mag_hot_bin
#print "vradial_hot_bin= ", vradial_hot_bin
#print "vtangential_hot_bin= ", vtangential_hot_bin
#print "v/sigma_radial = ",vradial_hot_bin/vradial_mag_hot_bin
#print "v/sigma_tangential = ",vtangential_hot_bin/vtangential_mag_hot_bin
#print "v_radial/v_tangential = ",vradial_mag_hot_bin/vtangential_mag_hot_bin

hdf5filename0 = sim + "/snapshot_" + snap_intag + "/snap_" + snap_intag + ".0.hdf5"
print "hdf5filename0 = ", hdf5filename0


acceleration =  "/PartType0/Acceleration" in  h5py.File(hdf5filename0)

print "ACCELERATION = ", acceleration
                             

vradial_hot, vtheta_hot, vphi_hot = bs.sphericalvels(coords_gas_shift[hot_indexes]-box_center,vels_gas_shift[hot_indexes])
vradial_cool, vtheta_cool, vphi_cool = bs.sphericalvels(coords_gas_shift[cool_indexes]-box_center,vels_gas_shift[cool_indexes])
vradial_all, vtheta_all, vphi_all = bs.sphericalvels(coords_gas_shift-box_center,vels_gas_shift)

mu =0.59
XH = 0.752

if(acceleration):
    acc_CGSConv = 3.24078e-15
    accel_gas = rg.accel(sniptag,sim,snap_intag,0)*acc_CGSConv  
    accel_gas_rot = rotg.rotatevector(accel_gas[gas_indexes],phi,theta,psi)
    accel_gas_shift = accel_gas_rot
    print "ACCELERATION, accel_gas_shift = ", accel_gas_shift
    accradial_hot_mag, acctangential_hot_mag, accradial_hot, acctangential_hot = vels_calc(coords_gas_shift[hot_indexes],accel_gas_shift[hot_indexes],box_center,[0,0,0])
    accradial_hot_bin = br.radialphysical(accradial_hot,logR_gas_hot,rlow,rhi,nradial)[0]
    print "ACCELERATION_bin, radial = ", accradial_hot_bin

    accradial_cool_mag, acctangential_cool_mag, accradial_cool, acctangential_cool = vels_calc(coords_gas_shift[cool_indexes],accel_gas_shift[cool_indexes],box_center,[0,0,0])
    accradial_cool_bin = br.radialphysical(accradial_cool,logR_gas_cool,rlow,rhi,nradial)[0]

    accradial_all_mag, acctangential_all_mag, accradial_all, acctangential_all = vels_calc(coords_gas_shift,accel_gas_shift,box_center,[0,0,0])
    accradial_all_bin = br.radialphysical(accradial_all,logR_gas,rlow,rhi,nradial)[0]

else:
    accradial_hot_bin = vradial_hot_bin*0.0
    accradial_cool_bin = vradial_cool_bin*0.0
    accradial_all_bin = vradial_all_bin*0.0
    accradial_hot = np.zeros(len(P[hot_indexes]))
    accradial_cool = np.zeros(len(P[cool_indexes]))
    accradial_all = np.zeros(len(P))


gas_mass_radial = br.radialcumsum(mass_gas,logR_gas,rlow,rhi,nradial)[0]
DM_mass_radial = br.radialcumsum(coords_DM[:,0]*0.0+mass_table[1]*1e+10/hubbleparam,logR_DM,rlow,rhi,nradial)[0]
stars_mass_radial = br.radialcumsum(mass_stars,logR_stars,rlow,rhi,nradial)[0]
Totmass_radial = (gas_mass_radial+DM_mass_radial+stars_mass_radial)*M_Solar

T_hot_bin, nH_hot_bin, P_hot_bin, S_hot_bin, rbins = radial_profile_calc(nH[hot_indexes], T[hot_indexes], logR_gas[hot_indexes],rlow,rhi,nradial)
T_cool_bin, nH_cool_bin, P_cool_bin, S_cool_bin, rbins = radial_profile_calc(nH[cool_indexes], T[cool_indexes], logR_gas[cool_indexes],rlow,rhi,nradial)
T_bin, nH_bin, P_bin, S_bin, rbins = radial_profile_calc(nH, T, logR_gas,rlow,rhi,nradial)

P_hot_3D, bin_edges_3D = bs.sphericalmean(P[hot_indexes], coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)
density_hot_3D = bs.sphericalmean(nH[hot_indexes]/XH/mu*M_P, coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vrad_hot_3D = bs.sphericalmean(vradial_hot, coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vtheta_hot_3D = bs.sphericalmean(vtheta_hot, coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vphi_hot_3D = bs.sphericalmean(vphi_hot, coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
acc_hot_3D = bs.sphericalmean(accradial_hot, coords_gas_shift[hot_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
label = "Hot"

hse_spherical.calc_Euler_spherical(10**bin_edges_3D[0]*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_hot_3D, density_hot_3D, vrad_hot_3D*1e+05, vtheta_hot_3D*1e+05, vphi_hot_3D*1e+05, acc_hot_3D, galidstr, label, M200, R200)
hse_spherical.calc_Euler_radial(10**bin_edges_3D[0]*cmpermpc, Totmass_radial, P_hot_bin, nH_hot_bin, vradial_hot_bin, vradial_mag_hot_bin, vtangential_hot_bin, vtangential_mag_hot_bin,accradial_hot_bin,galidstr,label,M200,R200)


P_cool_3D, bin_edges_3D = bs.sphericalmean(P[cool_indexes], coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)
density_cool_3D = bs.sphericalmean(nH[cool_indexes]/XH/mu*M_P, coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vrad_cool_3D = bs.sphericalmean(vradial_cool, coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vtheta_cool_3D = bs.sphericalmean(vtheta_cool, coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vphi_cool_3D = bs.sphericalmean(vphi_cool, coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
acc_cool_3D = bs.sphericalmean(accradial_cool, coords_gas_shift[cool_indexes]-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
label = "Cool"

hse_spherical.calc_Euler_spherical(10**bin_edges_3D[0]*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_cool_3D, density_cool_3D, vrad_cool_3D*1e+05, vtheta_cool_3D*1e+05, vphi_cool_3D*1e+05, acc_cool_3D, galidstr, label, M200, R200)

P_all_3D, bin_edges_3D = bs.sphericalmean(P, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)
density_all_3D = bs.sphericalmean(nH/XH/mu*M_P, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vrad_all_3D = bs.sphericalmean(vradial_all, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vtheta_all_3D = bs.sphericalmean(vtheta_all, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
vphi_all_3D = bs.sphericalmean(vphi_all, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
acc_all_3D = bs.sphericalmean(accradial_all, coords_gas_shift-box_center, rlow,rhi,nradial,ntheta,nphi)[0]
label = "All"

hse_spherical.calc_Euler_spherical(10**bin_edges_3D[0]*cmpermpc, bin_edges_3D[1], bin_edges_3D[2], Totmass_radial, P_all_3D, density_all_3D, vrad_all_3D*1e+05, vtheta_all_3D*1e+05, vphi_all_3D*1e+05, acc_all_3D, galidstr, label, M200, R200)



#print "bin_edges_3d = ", bin_edges_3D

#print "P_hot_3D = ", P_hot_3D
#for i in range(nradial):
#    for j in range(ntheta):
#        for k in range(nphi):
#            if(np.isnan(P_hot_3D[i][j][k])):
#                print "Missing P ", P_hot_3D[i][j][k], i, j, k
#            if((i==20) | (i==30)):
#                print "T,vrad,vtheta,vphi,acc= ", P_hot_3D[i][j][k]/density_hot_3D[i][j][k], vrad_hot_3D[i][j][k], vtheta_hot_3D[i][j][k], vphi_hot_3D[i][j][k], acc_hot_3D[i][j][k], i, j, k



gas_mass_bin,rbins = br.radialcumsum(mass_gas,logR_gas,rlow,rhi,nradial)
DM_mass_bin = br.radialcumsum(coords_DM[:,0]*0.0+mass_table[1]*1e+10/hubbleparam,logR_DM,rlow,rhi,nradial)[0]
stars_mass_bin = br.radialcumsum(mass_stars,logR_stars,rlow,rhi,nradial)[0]
hot_mass_bin = br.radialcumsum(mass_gas[hot_indexes],logR_gas[hot_indexes],rlow,rhi,nradial)[0]
cool_mass_bin = br.radialcumsum(mass_gas[cool_indexes],logR_gas[cool_indexes],rlow,rhi,nradial)[0]

print "Finished cumulative sums"

acc_hot_bin = P_hot_bin*0.0 

fradial = open("rbins_%s.dat"%galidstr,"w")
fradial_hot = open("rbins_hot_%s.dat"%galidstr,"w")
fradial_cool = open("rbins_cool_%s.dat"%galidstr,"w")

for i in xrange(len(gas_mass_bin)):
    print rbins[i], rbins[i+1],gas_mass_bin[i],T_hot_bin[i],nH_hot_bin[i],P_hot_bin[i],S_hot_bin[i]
    fradial_hot.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],hot_mass_bin[i],T_hot_bin[i],nH_hot_bin[i],P_hot_bin[i],S_hot_bin[i],vradial_hot_bin[i],vradial_mag_hot_bin[i],vtangential_hot_bin[i],vtangential_mag_hot_bin[i]))
    fradial_cool.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],cool_mass_bin[i],T_cool_bin[i],nH_cool_bin[i],P_cool_bin[i],S_cool_bin[i],vradial_cool_bin[i],vradial_mag_cool_bin[i],vtangential_cool_bin[i],vtangential_mag_cool_bin[i]))
    fradial.write("% 5.2f % 5.2f  %5.3e %5.3e %5.3e  %5.3e %5.3e %5.3e %5.3e %5.3e  % 6.1f %6.1f %6.1f %6.1f\n"%(rbins[i], rbins[i+1],DM_mass_bin[i],stars_mass_bin[i],gas_mass_bin[i],gas_mass_bin[i],T_bin[i],nH_bin[i],P_bin[i],S_bin[i],vradial_all_bin[i],vradial_mag_all_bin[i],vtangential_all_bin[i],vtangential_mag_all_bin[i]))

fradial.close()
fradial_hot.close()
fradial_cool.close()

print "Finished basic radial profiles"

vels_radial_plot(rbins, vradial_mag_hot_bin, vradial_hot_bin, vtangential_mag_hot_bin, vtangential_hot_bin, galidstr, "vel_hot", M200, R200, v200)

### DM velocities
vradial_DM_mag, vtangential_DM_mag, vradial_DM, vtangential_DM = vels_calc(coords_DM_shift,vels_DM_shift,box_center,[0,0,0])
print "v_DM_means (<R>,<tan>,R,tan) ", np.mean(vradial_DM_mag), np.mean(vtangential_DM_mag), np.mean(vradial_DM), np.mean(vtangential_DM)

vradial_DM_bin = br.radialphysical(vradial_DM,logR_DM,rlow,rhi,nradial)[0]
vradial_DM_mag_bin = br.radialphysical(vradial_DM_mag,logR_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_x = br.radialphysical(vtangential_DM[:,0],logR_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_y = br.radialphysical(vtangential_DM[:,1],logR_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin_z = br.radialphysical(vtangential_DM[:,2],logR_DM,rlow,rhi,nradial)[0]
vtangential_DM_mag_bin = br.radialphysical(vtangential_DM_mag,logR_DM,rlow,rhi,nradial)[0]
vtangential_DM_bin = np.sqrt(vtangential_DM_bin_x**2+vtangential_DM_bin_y**2+vtangential_DM_bin_z**2)

vels_radial_plot(rbins, vradial_DM_mag_bin, vradial_DM_bin, vtangential_DM_mag_bin, vtangential_DM_bin, galidstr, "vel_DM", M200, R200, v200)

print "Finished velocity radial profiles"


###
massfrac_radial_plot(rbins,DM_mass_bin,gas_mass_bin,stars_mass_bin,cool_mass_bin,hot_mass_bin,galidstr,M200,R200)

mu = 0.59
Ekinetic_cum_bin = br.radialcumsum(0.5*(vradial_mag[hot_indexes]**2+vtangential_mag[hot_indexes]**2)*1e+10*mass_gas[hot_indexes]*M_Solar, logR_gas[hot_indexes],rlow,rhi,nradial)[0]
Ethermal_cum_bin = br.radialcumsum(3/2.*K_Boltz*T[hot_indexes]/M_P/mu*(mass_gas[hot_indexes]*M_Solar), logR_gas[hot_indexes],rlow,rhi,nradial)[0]
Ethermal_cum_bin = br.radialcumsum(T[hot_indexes]*mass_gas[hot_indexes], logR_gas[hot_indexes],rlow,rhi,nradial)[0]
#print "Ethermal_cum_bin = ", Ethermal_cum_bin, 3/2.*K_Boltz/M_P/mu*M_Solar, K_Boltz, M_P, M_Solar
Ethermal_cum_bin = Ethermal_cum_bin*3/2.*(K_Boltz)/M_P/mu*M_Solar

#Egravity_cum_bin = br.radialcumsum(3/5.*G_Grav*tot logR_gas[hot_indexes],rlow,rhi,nradial)[0]
#Egravity_cum_bin = 3/5.*G_Grav*(DM_mass_bin+stars_mass_bin+gas_mass_bin)*M_Solar/rbins_M

energy_radial_plot(rbins,Ekinetic_cum_bin,Ethermal_cum_bin,(DM_mass_bin+stars_mass_bin+gas_mass_bin)*M_Solar,gas_mass_bin*M_Solar,galidstr,M200,R200)

print "Finished energy radial profiles"

Mtot, Mtherm = calc_hydrostatic_masses(rbins, DM_mass_bin+stars_mass_bin+gas_mass_bin, P_hot_bin, nH_hot_bin, vradial_hot_bin, vradial_mag_hot_bin, vtangential_hot_bin, vtangential_mag_hot_bin,accradial_hot_bin,galidstr,M200,R200)

print "Finished hydrostatic radial profiles"


#exit()

lgrid = 3*R200
#lgrid = 600/1e+03
lgridz = lgrid*4.
ngrid = 600

acc200 = (v200*1e+05)**2./(R200*cmpermpc)

if dofigs:
    for k in xrange(3):
        if(k == 0):
            theta = 0.
            phi = 90.
            dir = 'x'
        if(k == 1):
            theta = 90.
            phi = 0.
            dir = 'y'
        if(k == 2):
            theta = 0.
            phi = 0.
            dir = 'z'

        vels_gas_dir = vels_gas_rot[:,k] - gal_vels[k]
        sigsq_gas_dir = (vels_gas_rot[:,k] - gal_vels[k])**2.
        print "vels_gas_dir = ", vels_gas_dir
        print "sigsq_gas_dir = ", sigsq_gas_dir
        
        result_h = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.hydrogen.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=17.5, Vmax=21.5,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)
        
        result_hvel = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*vels_gas_dir, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.vel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)
        print "result_h = ", result_h
        print "result_hvel = ", result_hvel
        print "result_hvel/10**result_h", result_hvel/10**result_h
    
        coldens.make_colourmap('coldens.%s.vel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hvel/10**result_h, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

        if(acceleration):
            accel_gas_dir = accel_gas_rot[:,k]
            result_hacc = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*accel_gas_dir, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.acc.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-acc200*1e+19, Vmax=acc200*1e+19,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='Hacc',haloinfostr=haloinfostr,docbar=1)

            coldens.make_colourmap('coldens.%s.acc.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hacc/10**result_h, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -acc200, acc200, 'Acceleration', ngrid, redshift=redshift, extralabel='Acceleration',haloinfostr=haloinfostr,docbar=1)

        
        result_hsig = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*sigsq_gas_dir, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.sig.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=0, Vmax=v200,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='Hsig',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.sig.%s.l%3.1f.png'%(galidstr,dir,lgrid), np.sqrt(10**(result_hsig-result_h)), -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, 0, v200, 'Dispersion', ngrid, redshift=redshift, extralabel='Vel Disp.',haloinfostr=haloinfostr,docbar=1)
    
        if(chem):
            result_h1 = coldens.main(coords_gas_rot, hsmooth, mass_gas*hydrogen*h1, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.h1.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=13.0, Vmax=21.0,ion='HI',theta=theta,phi=phi,redshift=redshift,extralabel='HI',haloinfostr=haloinfostr,docbar=1)

        result_hcool = coldens.main(coords_gas_rot, hsmooth, cool_gas*hydrogen, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.hcool.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=17.5, Vmax=21.5,ion='H-cool',theta=theta,phi=phi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)

        result_hhot = coldens.main(coords_gas_rot, hsmooth, hot_gas*hydrogen, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.hhot.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=17.5, Vmax=21.5,ion='H-hot',theta=theta,phi=phi,redshift=redshift,extralabel='H',haloinfostr=haloinfostr,docbar=1)

        result_hcoolvel = coldens.main(coords_gas_rot, hsmooth, cool_gas*hydrogen*vels_gas_dir, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.coolvel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.coolvel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hcoolvel/10**result_hcool, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Cool Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

        result_hhotvel = coldens.main(coords_gas_rot, hsmooth, hot_gas*hydrogen*vels_gas_dir, gal_coords, lgrid, lgrid, lgridz, ngrid, ngrid, 58, boxsize_proper,fig_name='coldens.%s.hotvel.%s.l%3.1f.png'%(galidstr,dir,lgrid),Vmin=-v200*1e+19, Vmax=v200*1e+19,ion='H',theta=theta,phi=phi,redshift=redshift,extralabel='Hvel',haloinfostr=haloinfostr,docbar=1)

        coldens.make_colourmap('coldens.%s.hotvel.%s.l%3.1f.png'%(galidstr,dir,lgrid), result_hhotvel/10**result_hhot, -1.0*lgrid/2.0, lgrid/2.0, -1.0*lgrid/2.0, lgrid/2.0, -v200, v200, 'Hot Velocity', ngrid, redshift=redshift, extralabel='Velocity',haloinfostr=haloinfostr,docbar=1)

    
