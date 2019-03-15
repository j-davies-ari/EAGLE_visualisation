import sys
import tables
import h5py
import glob as glob
import os
import numpy as np
import scipy.stats as scist
import coldens_ben.coldens as coldens
import matplotlib.pyplot as plt
import matplotlib as mpl
import hse_spherical 

def slope_and_mean(R, val):
    val_slope = (val[-1]-val[0])/(R[-1]-R[0])
    val = np.mean(val)
    R = np.mean(R)

    return(val_slope,val,R)

def find_scale_and_norm(R,val,val_slope_inner,val_inner,R_inner,val_slope_outer,val_outer,R_outer):

    min_diff_slope = 1000.  
    for i in range(len(R)):
        val_slope_mean = (val_slope_inner + val_slope_outer)/2.
        if((R[i]>R_inner) & (R[i]<R_outer)):
            slope = (val[i+1]-val[i-1])/(R[i+1]-R[i-1])
            if(np.fabs(slope-val_slope_mean)<min_diff_slope):
                min_diff_slope = np.fabs(slope-val_slope_mean)
                middle_slope = slope
                i_min_diff_slope = i

    print "FIT: ", R[i_min_diff_slope],val[i_min_diff_slope],middle_slope,val_slope_inner,val_slope_outer, val_slope_mean

    val_0 = 10**val[i_min_diff_slope]*(1)**-val_slope_inner*(2)**-(val_slope_outer-val_slope_inner)
    val_0 = np.log10(val_0)
    print "RETURN: ", val_0, R[i_min_diff_slope]
    return(R[i_min_diff_slope],val_0)

def return_high_oxygen_ion_fractions(nH, T):

    #nH_table, T_table, fOVI_table, fOVII_table, fOVIII_table = np.loadtxt(oxygen_table_file)
    oxygen_table_file = "/net/virgo.strw.leidenuniv.nl/data2/oppenheimer/ionfiles/lt02HMQG.CIE.Oxygen.short"
    T_table, fOVI_table, fOVII_table, fOVIII_table = np.loadtxt(oxygen_table_file, usecols=(0,6,7,8),unpack=True)
    
    fOVI = np.interp(T,T_table, fOVI_table)
    fOVII = np.interp(T,T_table, fOVII_table)
    fOVIII = np.interp(T,T_table, fOVIII_table)

    return(10**fOVI,10**fOVII,10**fOVIII)

def McDonald_2017_nH_MH150(lM_200,lR,R_200):

    R200_to_R500_ratio = 1.54
    M_500 = 10**lM_200*(500/200)*(R200_to_R500_ratio)**-3.

    
    R_200_kpc = R_200*1e+03
    R_500_kpc = R_200_kpc/R200_to_R500_ratio
    R_kpc = 10**lR * R_200_kpc
    
    R_scale_500 = 10**lR*R200_to_R500_ratio
    R_scale_500_data, lrho_crit_500_data, err_500 = np.loadtxt("/home/oppenheimer/hse/profiles/McDonald_2017.rhocrit_r500.dat",usecols=(0,1,2),unpack=True)
    lrho_crit_500 = np.interp(R_scale_500,R_scale_500_data,lrho_crit_500_data)

    rhocrit = rhocrit_cgs/ M_P
    
    nH_cgs = 10**lrho_crit_500*rhocrit*XH

    print "McD= ", rhocrit, lrho_crit_500, nH_cgs
    
    M_500_hot = calc_spherical_mass(nH_cgs,R_kpc,R_500_kpc)
    print "M_500_hot= ", M_500, np.log10(M_500_hot), M_500_hot/M_500

    R_kpc_return = np.linspace(10.,R_200_kpc,num=50)
    nH_cgs_return = np.interp(R_kpc_return,R_kpc,nH_cgs)
    print "McDonald = ", R_kpc_return, nH_cgs_return
    return(R_kpc_return,nH_cgs_return)

def Miller_2015_nH_fit_MW(lM_200,lR,R_200):

    R_200_kpc = R_200*1e+03
    r_kpc = 10**lR * R_200_kpc
    ne_over_nh = 1.16
    
    n0_rc_normed = 1.35e-02
    correction = 1./0.3
    #correction = 1.1e-04/1.54e-05 # Corrected to Faerman LMC e^-2 value. 
    beta= 0.50
    
    n_e = n0_rc_normed*correction/(r_kpc)**(3*beta)

    print "r_kpc = ", r_kpc
    print "n_e = ", n_e

    nH = n_e/ne_over_nh

    M_200_hot = calc_spherical_mass(nH,r_kpc,R_200_kpc)
    print "M_200_hot_Miller = ", M_200_hot
    return(nH)

    

def Bregman_2018_nH_fit_MW(lM_200,lR,R_200):
    #n(r) = n_0 r_c^(3*beta)/r^(3*beta)
    #n_0 r_c^(3*beta) = 1.20e-02 (+2.13,-0.82)
    #beta = 0.56 (+0.10,-0.12)
    # Is this electron density?  XXX

    R_200_kpc = R_200*1e+03
    r_kpc = 10**lR * R_200_kpc
    ne_over_nh = 1.16
    
    n0_rc_normed = 1.2e-02
    correction = 1.1e-04/1.54e-05 # Corrected to Faerman LMC e^-2 value. 
    beta= 0.56
    
    n_e = n0_rc_normed*correction/(r_kpc)**(3*beta)

    print "r_kpc = ", r_kpc
    print "n_e = ", n_e

    nH = n_e/ne_over_nh

    M_200_hot = calc_spherical_mass(nH,r_kpc,R_200_kpc)
    print "M_200_Bregman2018 = ", M_200_hot
    return(nH)

    
def calc_spherical_mass(nH,r_kpc,R_outer_kpc):

    M_sum = 0.0
    r_cgs = r_kpc*cmperkpc
    
    for i in range(len(r_kpc)):
        if((i>0) & (r_kpc[i]<R_outer_kpc)):
            M_sum += 4*np.pi/3.*(r_cgs[i]**3-r_cgs[i-1]**3)*(nH[i]+nH[i-1])/2.*M_P/XH
        if((r_kpc[i]>R_outer_kpc) & (r_kpc[i-1]<R_outer_kpc)):
            nH_200 = 10**(np.log10(nH[i-1]) + (np.log10(nH[i])-np.log10(nH[i-1]))*(R_outer_kpc-r_kpc[i-1])/(r_kpc[i]-r_kpc[i-1]))
            M_sum_old = M_sum
            M_sum += 4*np.pi/3.*((R_outer_kpc*cmperkpc)**3-r_cgs[i-1]**3)*(nH_200+nH[i-1])/2.*M_P/XH
            print "Last integration= ", nH[i-1],nH[i],nH_200,r_kpc[i-1],r_kpc[i],R_outer_kpc,M_sum,M_sum/M_sum_old
           
    M_sum /= 2.0e+33
    return(M_sum)

def calc_density_vector(cum_mass,r_cgs):
    #r_cgs = r_kpc*cmperkpc
    cum_mass_cgs = cum_mass*M_Solar

    cum_density = r_cgs*0.0
    
    for i in range(len(r_cgs)):
        if((i>0) & i<len(r_cgs)):
            cum_density[i] = (cum_mass_cgs[i]-cum_mass_cgs[i-1])/(4*np.pi/3.*(r_cgs[i]**3-r_cgs[i-1]**3))
            #print cum_mass_cgs[i]-cum_mass_cgs[i-1], (4*np.pi/3.*(r_cgs[i]**3-r_cgs[i-1]**3)),r_cgs[i],r_cgs[i]
            
    return(cum_density)
            
def Arnaud_2010_P_fit_z0(lM_200,lR):

    R200_to_R500_ratio = 1.54
    M_500 = 10**lM_200*(500/200)*(R200_to_R500_ratio)**-3.
    
    R_scale_500 = 10**lR*R200_to_R500_ratio

    x = R_scale_500
    print "R_scale_500 = ", R_scale_500, lM_200, np.log10(M_500)
    keVtoK = 8.6e-08
    
    h_70 = 0.6777/0.7
    Mass_exponent = (2/3.+(0.12+0.10)-(0.12+0.10)*(x/0.5)**3/(1+(x/0.5)**3))

    print "Mass_exponent = ", Mass_exponent

    P_0 = 8.403*h_70**(-3/2.)
    c_500 = 1.177
    gamma = 0.3081
    alpha = 1.0510
    beta = 5.4905
    
    P_e = 1.65e-03*(M_500/(3e+14/h_70))**Mass_exponent * P_0*h_70**2/((c_500*x)**gamma*(1+(c_500*x)**alpha)**((beta-gamma)/alpha))/keVtoK  

    n_over_ne = 1/(1.16*0.752*0.59)

    P = P_e*n_over_ne
    print "Arnaud = ", lM_200, M_500, P, x, lR
    
    return(P)

def GNFW_fit(lR, alpha, beta, gamma, val_norm, c_200):

    R200_to_R500_ratio = 1.54
    #M_500 = 10**lM_200*(500/200)*(R200_to_R500_ratio)**-3.
    
    R_scale_500 = 10**lR*R200_to_R500_ratio
    c_500 = c_200*R200_to_R500_ratio 

    
    print "R= ", R_scale_500
    
    x = R_scale_500
    #print "R_scale_500 = ", R_scale_500, lM_200, np.log10(M_500)
    #keVtoK = 8.6e-08

    alpha = -alpha
    beta = -beta
    gamma = -gamma
    print "a,b,g= ", alpha, beta, gamma

    P = val_norm/((c_500*x)**gamma*(1+(c_500*x)**alpha)**((beta-gamma)/alpha))

    print "P= ", P

    return(P)
    
G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24
cmperkpc = 3.086e+21
XH = 0.752
omega_b = 0.04825
omega_M = 0.307
mu = 0.59
OtoH_Solar = 10**-3.31

hubbleparam = 0.6777
redshift = 0.0
omegaM = 0.307
omegaL = 1-omegaM
omegaratio = (omegaM+omegaL/(1+redshift)**3) 

rhocrit_cgs = 1.88e-29*hubbleparam**2*(1+redshift)**3

lRlow = -1.5

#fradial_file = sys.argv[1]
lsfilename = sys.argv[1]
fit = float(sys.argv[2])

lsfile = open(lsfilename,"r").readlines()

nhalos = len(lsfile)

fig_Phot_norm = plt.figure(figsize=(4.5,3.375))
ax_Phot_norm = fig_Phot_norm.add_subplot(111)

fig_nHhot_norm = plt.figure(figsize=(4.5,3.375))
ax_nHhot_norm = fig_nHhot_norm.add_subplot(111)

fig_Phot_all = plt.figure(figsize=(4.5,6))
ax_Phot_all = fig_Phot_all.add_subplot(111)

fig_nHhot_all = plt.figure(figsize=(4.5,6))
ax_nHhot_all = fig_nHhot_all.add_subplot(111)

fig_NOVI_all = plt.figure(figsize=(4.5,3.375))
ax_NOVI_all = fig_NOVI_all.add_subplot(111)

fig_NOVII_all = plt.figure(figsize=(4.5,3.375))
ax_NOVII_all = fig_NOVII_all.add_subplot(111)

fig_NOVIII_all = plt.figure(figsize=(4.5,3.375))
ax_NOVIII_all = fig_NOVIII_all.add_subplot(111)

fig_fhot_halo = plt.figure(figsize=(4.5,3.375))
ax_fhot_halo = fig_fhot_halo.add_subplot(111)

for h in range(nhalos):
    fradial_file = lsfile[h].split()[0]
    colour = lsfile[h].split()[1]
    label = lsfile[h].split()[2]
    colour_2 = lsfile[h].split()[3]
    put_label = lsfile[h].split()[4]
    alpha_plot = float(lsfile[h].split()[5])
    
    lR, Mtot, lnH_norm, lT_norm, lZ, lP_norm, lS_norm, cum_DM, fbar, fhot = np.loadtxt(fradial_file,usecols=(0,1,2,3,4,5,6,7,8,9),unpack=True)

    halolabel1 = fradial_file.split(".")[2]
    halolabel2 = fradial_file.split(".")[3]
    halolabel = halolabel1 + "." + halolabel2
    
    P_Arnaud_MH150_cgs = Arnaud_2010_P_fit_z0(15.0,lR)

    R_200_MH150_Mpc = 1.63e-5*(10**15.0*hubbleparam)**0.333/hubbleparam
    R_McDonald_MH150_kpc,nH_McDonald_MH150_cgs = McDonald_2017_nH_MH150(15.0,lR,R_200_MH150_Mpc)
    
    # Calculate slopes
    indexes_inner = np.where((lR>-1.35) & (lR<-1.05))
    indexes_outer = np.where((lR>0.15) & (lR<0.45))
    #indexes_inner = np.where((lR>-1.45) & (lR<-0.85))
    #indexes_outer = np.where((lR>-0.15) & (lR<0.45))

    lP_slope_inner, lP_inner, lR_inner = slope_and_mean(lR[indexes_inner],lP_norm[indexes_inner])
    lP_slope_outer, lP_outer, lR_outer = slope_and_mean(lR[indexes_outer],lP_norm[indexes_outer])

    lnH_slope_inner, lnH_inner, lR_inner = slope_and_mean(lR[indexes_inner],lnH_norm[indexes_inner])
    lnH_slope_outer, lnH_outer, lR_outer = slope_and_mean(lR[indexes_outer],lnH_norm[indexes_outer])

    lR_P_scale, lP_0 =  find_scale_and_norm(lR,lP_norm,lP_slope_inner,lP_inner,lR_inner,lP_slope_outer,lP_outer,lR_outer)
    
    lR_nH_scale, lnH_0 =  find_scale_and_norm(lR,lnH_norm,lnH_slope_inner,lnH_inner,lR_inner,lnH_slope_outer,lnH_outer,lR_outer)

    lP_fit = np.zeros(len(lR))
    lnH_fit = np.zeros(len(lR))

    print "lP_0= ", lP_0, " lnH_0= ", lnH_0
    
    for i in range(len(lR)):
        #print lP_0, 10**lR[i], 10**lR_P_scale, P_slope_inner, (lP_slope_outer-lP_slope_inner), (10**lR[i]/10**lR_P_scale), (1+10**lR[i]/10**lR_P_scale)
        lP_fit[i] = np.log10(10**lP_0/((10**lR[i]/10**lR_P_scale)**-lP_slope_inner*(1+10**lR[i]/10**lR_P_scale)**-(lP_slope_outer-lP_slope_inner)))
        lnH_fit[i] = np.log10(10**lnH_0/((10**lR[i]/10**lR_nH_scale)**-lnH_slope_inner*(1+10**lR[i]/10**lR_nH_scale)**-(lnH_slope_outer-lnH_slope_inner)))
        #print "Pressure_nH_fit = ", lR[i], lP_norm[i], lP_fit[i], ;lnH_norm[i], lnH_fit[i]
    
    #P_slope_inner = (lP_norm[indexes_inner][-1]-lP_norm[indexes_inner][0])/(lR[indexes_inner][-1]-lR[indexes_inner][0])
    #P_inner = np.mean(lP_norm[indexes_inner])
    #lR_inner = np.mean(lR[indexes_inner])

    print "Pinner= ", lP_slope_inner, lP_inner, lR_inner, lR[indexes_inner], lP_norm[indexes_inner], len(lR[indexes_inner])
    print "Pouter= ", lP_slope_outer, lP_outer, lR_outer, lR[indexes_outer], lP_norm[indexes_outer], len(lR[indexes_outer])

    #alpha at r_s, beta at r>>r_s, gamma at r<<r_s
    #P_GNFW = np.log10(GNFW_fit(lR, lP_slope_inner, lP_slope_outer-0.5, lP_slope_inner+0.5, lP_0, 1.0))
    #nH_GNFW = np.log10(GNFW_fit(lR, lnH_slope_inner, lnH_slope_outer-0.5, lnH_slope_inner+0.5, lnH_0, 1.0))
    lP_GNFW = np.log10(GNFW_fit(lR, (lP_slope_inner+lP_slope_outer)/2., lP_slope_outer, lP_slope_inner, lP_0, 10**lR_P_scale))
    lnH_GNFW = np.log10(GNFW_fit(lR, (lnH_slope_inner+lnH_slope_outer)/2., lnH_slope_outer, lnH_slope_inner, lnH_0, 10**lR_nH_scale))
    
    print "lP_GNFW= ", lP_GNFW
    
    print "nHinner= ", lnH_slope_inner, lnH_inner, lR_inner, lR[indexes_inner], lnH_norm[indexes_inner], len(lR[indexes_inner])
    print "nHouter= ", lnH_slope_outer, lnH_outer, lR_outer, lR[indexes_outer], lnH_norm[indexes_outer], len(lR[indexes_outer])

    indexes_mass_inner  =  np.where((lR>-1.65) & (lR<-1.35))
    Mtot_slope_inner, Mtot_inner, lR_inner = slope_and_mean(lR[indexes_inner],Mtot[indexes_inner])
    Mtot_slope_outer, Mtot_outer, lR_outer = slope_and_mean(lR[indexes_outer],Mtot[indexes_outer])

    #print "Mtotinner= ", Mtot_slope_inner, Mtot_inner, lR_inner, lR[indexes_inner], Mtot[indexes_inner], len(lR[indexes_inner])
    #print "Mtotouter= ", Mtot_slope_outer, Mtot_outer, lR_outer, lR[indexes_outer], Mtot[indexes_outer], len(lR[indexes_outer])


    indexes_R200 = np.where((lR>-0.05) & (lR<0.05))
#    lM_200 = Mtot[indexes_R200]

    lM_200 = np.interp(0.0,lR,Mtot)
    
    R_200_Mpc = 1.63e-5*(10**lM_200*hubbleparam)**0.333/omegaratio**0.333/(1+redshift)/hubbleparam
    R_200_kpc = R_200_Mpc*1e+03
    lTvir = 5.69 + 2/3.*(lM_200-12.0)+ 1/3.*np.log10(omegaratio) + np.log10(1+redshift)
    lrhocrit200 = np.log10(rhocrit_cgs*omega_b/omega_M*200*XH/M_P)

    labelwrite = "%s log[$M_{\mathrm{200}}$]$=%s$"%(label,halolabel)
    
    print "lTvir, lrhocrit200= ", lTvir, lrhocrit200

    R_cgs = R_200_Mpc*10**lR*cmpermpc

    if(fit == 0):
        lnH_norm = lnH_norm
        lP_nrom = lP_norm
    else:
        lnH_norm = lnH_GNFW #  nH_fit
        lP_norm = lP_GNFW   #  P_fit

    ax_Phot_norm.plot(10**lR, lP_norm, color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)

    ax_nHhot_norm.plot(10**lR, lnH_norm, color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)

    nH_cgs = 10**lrhocrit200*10**lnH_norm
    P_cgs = 10**lrhocrit200*10**lTvir*10**lP_norm
    
    nDM_cgs = calc_density_vector(10**cum_DM,R_cgs)

    print "nDM_cgs= ", nDM_cgs
    nDM_H = nDM_cgs*omega_b/(omega_M-omega_b)*XH/M_P
    print "nDM_H= ", nDM_H
                 
    if((lM_200>11.7) & (lM_200<12.3)):
        lM_200_MW = 12.2
        R_200_MW_Mpc = 1.63e-5*(10**lM_200_MW*hubbleparam)**0.333/hubbleparam
        nH_Bregman_MW_cgs = Bregman_2018_nH_fit_MW(lM_200_MW,lR,R_200_Mpc)
        lnH_Bregman_MW =  np.log10(nH_Bregman_MW_cgs/(10**lrhocrit200))
        print "nH_Bregman_MW = ", lnH_Bregman_MW, nH_Bregman_MW_cgs
        #ax_nHhot_norm.plot(10**lR, lnH_Bregman_MW, color="cyan",lw=4,alpha=1.0,ls=':',zorder=20+h)
        #ax_nHhot_all.plot(10**lR, lnH_Bregman_MW+lrhocrit200, color="cyan",lw=4,alpha=1.0,ls=':',zorder=20+h)
        #Miller_2015_nH_fit
        nH_Miller_MW_cgs = Miller_2015_nH_fit_MW(lM_200_MW,lR,R_200_Mpc)
        lnH_Miller_MW =  np.log10(nH_Miller_MW_cgs/(10**lrhocrit200))
        print "nH_Miller_MW = ", lnH_Miller_MW, nH_Miller_MW_cgs
        ax_nHhot_norm.plot(10**lR, lnH_Miller_MW, color="cyan",lw=2,alpha=1.0,ls='--',zorder=20+h, label="Miller 2015")
        ax_nHhot_all.plot(10**lR, lnH_Miller_MW+lrhocrit200, color="cyan",lw=2,alpha=1.0,ls='--',zorder=20+h, label="Miller 2015 MW")

        
    if(lM_200>14.8):
        P_0 = 3.0*1.65e-03/8.6e-08/(1.16*0.752*0.59)/(10**lrhocrit200*10**lTvir)
        c_500 = 1.177
        #gamma = 0.3081
        #alpha = 1.0510
        #beta = 5.4905
        print "P_0 = ", P_0
        
        lP_GNFW = np.log10(GNFW_fit(lR, -1.0510, -5.4905, -0.3081, P_0, c_500/1.54))
        ###ax_Phot_norm.plot(10**lR, lP_GNFW, color="green", lw=4,alpha=1.0,ls=':',zorder=20-h)
        lP_GNFW = np.log10(GNFW_fit(lR, -1.0510, -5.4905, 0.00, P_0, c_500/1.54))
        ###ax_Phot_norm.plot(10**lR, lP_GNFW, color="blue", lw=4,alpha=1.0,ls=':',zorder=20-h)
        lP_GNFW = np.log10(GNFW_fit(lR, -1.0510, -5.4905, 0.3081, P_0, c_500/1.54))
        ###ax_Phot_norm.plot(10**lR, lP_GNFW, color="purple", lw=4,alpha=1.0,ls=':',zorder=20-h)

        nH_0 = 3e-03/(10**lrhocrit200)
        c_500 = 1.54
        
        lnH_GNFW = np.log10(GNFW_fit(lR, -0.85, -3.5, 0.0, nH_0*2, c_500/1.54))
        ###ax_nHhot_norm.plot(10**lR, lnH_GNFW, color="green", lw=4,alpha=1.0,ls=':',zorder=20-h)
        lnH_GNFW = np.log10(GNFW_fit(lR, -0.85, -3.0, 0.0, nH_0*1, c_500/1.54))
        ###ax_nHhot_norm.plot(10**lR, lnH_GNFW, color="blue", lw=4,alpha=1.0,ls=':',zorder=20-h)
        lnH_GNFW = np.log10(GNFW_fit(lR, -0.85, -2.5, 0.0, nH_0*0.5, c_500/1.54))
        ###ax_nHhot_norm.plot(10**lR, lnH_GNFW, color="purple", lw=4,alpha=1.0,ls=':',zorder=20-h)


    #ax_Phot_all.plot(10**lR, np.log10(nDM_H)+lTvir, color=colour, lw=1,alpha=1.0,ls=':',zorder=10-h)
    
    lTvir_MH150 = 5.69 + 2/3.*(15.0-12.0)+ 1/3.*np.log10(omegaratio) + np.log10(1+redshift)
    #lP_Arnaud_MH150 =  np.log10(P_Arnaud_MH150_cgs/(10**lrhocrit200*10**lTvir_MH150))
    #ax_Phot_norm.plot(10**lR, lP_Arnaud_MH150, color="purple", lw=4,alpha=1.0,ls=':',zorder=20+h)
    #ax_Phot_all.plot(10**lR, lP_Arnaud_MH150+lrhocrit200+lTvir_MH150, color="purple", lw=4,alpha=1.0,ls=':',zorder=20+h)
    
    if(lM_200>14.8):
        ax_Phot_norm.plot(10**lR, np.log10(P_Arnaud_MH150_cgs)-lrhocrit200-lTvir_MH150, color="purple", lw=2,alpha=1.0,ls=':',zorder=20+h,label="Arnaud 2010 Clusters")
        ax_nHhot_norm.plot(R_McDonald_MH150_kpc/R_200_kpc*(R_200_Mpc/R_200_MH150_Mpc), np.log10(nH_McDonald_MH150_cgs)-lrhocrit200, color="purple",lw=2,alpha=1.0,ls=':',zorder=20+h,label="McDonald 2017 Clusters")


    ax_Phot_all.plot(10**lR, np.log10(P_cgs), color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)

    if(lM_200>14.8):    
        ax_Phot_all.plot(10**lR, np.log10(P_Arnaud_MH150_cgs)+lTvir-lTvir_MH150, color=colour, lw=2,alpha=1.0,ls=':',zorder=20+h,label="Arnaud 2010 Clusters")
    else:
        ax_Phot_all.plot(10**lR, np.log10(P_Arnaud_MH150_cgs)+lTvir-lTvir_MH150, color=colour, lw=2,alpha=1.0,ls=':',zorder=20+h)
        
    
    #ax_nHhot_all.plot(10**lR,np.log10(nDM_H)+(lM_200-12.0),color=colour,lw=1,ls=':',zorder=10-h)
    
    ax_nHhot_all.plot(10**lR, np.log10(nH_cgs)+(lM_200-12.0), color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)

    if(lM_200>14.8):
        ax_nHhot_all.plot(R_McDonald_MH150_kpc/R_200_kpc*(R_200_Mpc/R_200_MH150_Mpc),np.log10(nH_McDonald_MH150_cgs)+(lM_200-12.0), color=colour,lw=2,ls=':',zorder=20+h,label="McDonald 2017 Clusters")
    else:
        ax_nHhot_all.plot(R_McDonald_MH150_kpc/R_200_kpc*(R_200_Mpc/R_200_MH150_Mpc),np.log10(nH_McDonald_MH150_cgs)+(lM_200-12.0), color=colour,lw=2,ls=':',zorder=20+h)
    
    print "nH_cgs = ",nH_cgs
    print "R_cgs = ",R_cgs
    print "R_200_kpc = ",R_200_kpc
    
    M_200_hot = calc_spherical_mass(nH_cgs,R_cgs/cmperkpc,R_200_kpc)
    print "lM_200_hot= ", lM_200, np.log10(M_200_hot), M_200_hot/10**lM_200

    R200_to_R500_ratio = 1.54
    M_500 = 10**lM_200*(500/200)*(R200_to_R500_ratio)**-3.
    R_500_kpc = R_200_kpc/R200_to_R500_ratio
    M_500_hot = calc_spherical_mass(nH_cgs,R_cgs/cmperkpc,R_500_kpc)
    print "M_500_hot= ", np.log10(M_500), np.log10(M_500_hot), M_500_hot/M_500

    fhot_200 = np.interp(R_200_kpc,R_cgs/cmperkpc,fhot)
    fhot_500 = np.interp(R_500_kpc,R_cgs/cmperkpc,fhot)

    print "put_label= ", put_label, h
    if(h==0):
        fhot_500_label = ax_fhot_halo.plot(lM_200,fhot_500,color='k',marker="s",ls='')
        fhot_200_label = ax_fhot_halo.plot(lM_200,fhot_200,color='k',marker="o",ls='')

    if(float(put_label)):
        ax_fhot_halo.plot(lM_200,fhot_500,color=colour_2,marker="s",ls='',label=label,alpha=alpha_plot)
    else:
        ax_fhot_halo.plot(lM_200,fhot_500,color=colour_2,marker="s",ls='',alpha=alpha_plot)

    ax_fhot_halo.plot(lM_200,fhot_200,color=colour_2,marker="o",ls='',alpha=alpha_plot)
    
    T_cgs = P_cgs/(nH_cgs/XH)*mu

    nO_cgs = nH_cgs*10**(lZ+1.87)*OtoH_Solar
    
    f_OVI, f_OVII, f_OVIII = return_high_oxygen_ion_fractions(np.log10(nH_cgs),np.log10(T_cgs))

    #print "Oxygen ion fractions= ", T_cgs, f_OVI, f_OVII, f_OVIII
    
    #for i in range(len(lR)):
        #print "Oxygen= ", lR[i], nO_cgs[i], f_OVI[i]*nO_cgs[i], f_OVII[i]*nO_cgs[i], f_OVIII[i]*nO_cgs[i]
    
    b_kpc = np.linspace(10,1200,num=120)
    b_kpc_halfstep = np.linspace(5,1195,num=120)
    nOVI_kpc = np.interp(b_kpc,R_cgs/cmperkpc,f_OVI*nO_cgs)
    nOVII_kpc = np.interp(b_kpc,R_cgs/cmperkpc,f_OVII*nO_cgs)
    nOVIII_kpc = np.interp(b_kpc,R_cgs/cmperkpc,f_OVIII*nO_cgs)
    N_OVI = np.zeros(len(b_kpc))
    N_OVII = np.zeros(len(b_kpc))
    N_OVIII = np.zeros(len(b_kpc))

    for i in range(len(b_kpc)):
        for j in range(len(b_kpc)):
            #if((b_kpc[j]>=b_kpc[i]) & (b_kpc[j]<b_kpc[-1])):
            if((b_kpc[j]>=b_kpc_halfstep[i]) & (b_kpc[j]<b_kpc_halfstep[-1])):
                path = (b_kpc_halfstep[j+1]-b_kpc_halfstep[j])*b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
                N_OVI[i] += 2*path*cmperkpc*nOVI_kpc[j]
                N_OVII[i] += 2*path*cmperkpc*nOVII_kpc[j]
                N_OVIII[i] += 2*path*cmperkpc*nOVIII_kpc[j]
                #print "PATHLENGTH= ", i,j, path, b_kpc_halfstep[j+1], b_kpc_halfstep[j], b_kpc[i], b_kpc_halfstep[j+1]/np.sqrt(b_kpc_halfstep[j+1]**2-b_kpc[i]**2)
        if(b_kpc[i]<=300): print "radial_profile= ", b_kpc[i], np.log10(N_OVI[i]),np.log10(N_OVII[i]),np.log10(N_OVIII[i])

    N_OVI_150_sum = 0.0
    N_OVII_150_sum = 0.0
    N_OVIII_150_sum = 0.0
    area = 0.0

    for i in range(len(b_kpc)):
        if(b_kpc[i]<=150):
            N_OVI_150_sum += N_OVI[i]*np.pi*(b_kpc_halfstep[i+1]**2-b_kpc_halfstep[i]**2)
            N_OVII_150_sum += N_OVII[i]*np.pi*(b_kpc_halfstep[i+1]**2-b_kpc_halfstep[i]**2)
            N_OVIII_150_sum += N_OVIII[i]*np.pi*(b_kpc_halfstep[i+1]**2-b_kpc_halfstep[i]**2)
            area += np.pi*(b_kpc_halfstep[i+1]**2-b_kpc_halfstep[i]**2)

    N_OVI_150_sum /= area
    N_OVII_150_sum /= area
    N_OVIII_150_sum /= area

    print "<N_150kpc>= ", np.log10(N_OVI_150_sum), np.log10(N_OVII_150_sum), np.log10(N_OVIII_150_sum)

    ax_NOVI_all.plot(b_kpc, np.log10(N_OVI), color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label="%s log[$N_{OVI}$] = %4.1f"%(labelwrite,np.log10(N_OVI_150_sum)))
    ax_NOVII_all.plot(b_kpc, np.log10(N_OVII), color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label="%s log[$N_{OVII}$] = %4.1f"%(labelwrite,np.log10(N_OVII_150_sum)))
    ax_NOVIII_all.plot(b_kpc, np.log10(N_OVIII), color=colour, lw=4,alpha=1.0,ls='-',zorder=20-h,label="%s log[$N_{OVIII}$] = %4.1f"%(labelwrite,np.log10(N_OVIII_150_sum)))
    
    #Do Euler Equation:              
    for i in range(len(lR)):
        if((i > 0) & (i < len(lR)-1)):
            Mgrav = 10**Mtot[i]*M_Solar
        
            Mtherm = -1/(4*np.pi*G_Grav)*(4*np.pi*R_cgs[i]**2)/(nH_cgs[i]/XH*M_P)*(P_cgs[i-1]-P_cgs[i+1])*K_Boltz/(R_cgs[i-1]-R_cgs[i+1])
            #print R_cgs[i], nH_cgs[i], P_cgs[i+1],P_cgs[i-1], nH[i], 10**nH[i]
            #print np.log10(R_cgs[i]/R_200), np.log10(Mgrav/M_Solar), Mtherm, Mtherm/Mgrav
            ###print "EulerRatio: ", lR[i], np.log10(Mgrav/M_Solar), np.log10(T_cgs[i]) , Mtherm, Mtherm/Mgrav

if(fit == 0):
    lsfilename = lsfilename 
else:
    lsfilename = lsfilename + "_fit"

    
ax_Phot_norm.set_xscale('log')
ax_Phot_norm.set_xlim(10**lRlow,10**0.33)
ax_Phot_norm.set_ylim(-3.0,3.0)
ax_Phot_norm.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Phot_norm.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Phot_norm.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Phot_norm.set_ylabel('log $P/P_{200}$',fontsize=16)
ax_Phot_norm.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Phot_norm.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Phot_norm.legend(loc="upper right", ncol=1,fontsize=10)
fig_Phot_norm.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_Phot_norm.savefig('P_hot_norm.' + lsfilename + '.png')

ax_nHhot_norm.set_xscale('log')
ax_nHhot_norm.set_xlim(10**lRlow,10**0.33)
ax_nHhot_norm.set_ylim(-2.0,2.0)
ax_nHhot_norm.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_nHhot_norm.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_nHhot_norm.set_xlabel('$R/R_{200}$',fontsize=16)
ax_nHhot_norm.set_ylabel('log $n_\mathrm{H}/n_{\mathrm{H},200}$',fontsize=16)
ax_nHhot_norm.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_nHhot_norm.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_nHhot_norm.legend(loc="upper right", ncol=1,fontsize=10)
fig_nHhot_norm.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_nHhot_norm.savefig('nH_hot_norm.' + lsfilename + '.png')

ax_Phot_all.set_xscale('log')
ax_Phot_all.set_xlim(10**lRlow,10**0.33)
ax_Phot_all.set_ylim(1.0,8.5)
ax_Phot_all.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Phot_all.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Phot_all.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Phot_all.set_ylabel('log $P$ [cm$^{-3}$ K]',fontsize=16)
ax_Phot_all.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Phot_all.legend(loc="upper right", ncol=1,fontsize=10)
fig_Phot_all.subplots_adjust(left=0.17, bottom=0.10,top=0.98,right=0.98)
fig_Phot_all.savefig('P_hot_all.' + lsfilename + '.png')

ax_nHhot_all.set_xscale('log')
ax_nHhot_all.set_xlim(10**lRlow,10**0.33)
ax_nHhot_all.set_ylim(-5.0,4.0)
ax_nHhot_all.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_nHhot_all.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_nHhot_all.set_xlabel('$R/R_{200}$',fontsize=16)
ax_nHhot_all.set_ylabel(r'log $n_\mathrm{H} \times M_{200}/10^{12}$ [cm$^{-2}$]',fontsize=16)
ax_nHhot_all.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_nHhot_all.legend(loc="upper right", ncol=1,fontsize=10)
fig_nHhot_all.subplots_adjust(left=0.17, bottom=0.10,top=0.98,right=0.98)
fig_nHhot_all.savefig('nH_hot_all.' + lsfilename + '.png')


ax_NOVI_all.set_xlim(0,300)
ax_NOVI_all.set_ylim(12,18)
ax_NOVI_all.set_xlabel('$b$ [kpc]',fontsize=16)
ax_NOVI_all.set_ylabel('log $N_{OVI}$ [cm$^{-2}$]',fontsize=16)
ax_NOVI_all.legend(loc="upper right", ncol=1,fontsize=10)
fig_NOVI_all.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_NOVI_all.savefig('NOVI_all.' + lsfilename + '.png')

ax_NOVII_all.set_xlim(0,300)
ax_NOVII_all.set_ylim(12,18)
ax_NOVII_all.set_xlabel('$b$ [kpc]',fontsize=16)
ax_NOVII_all.set_ylabel('log $N_{OVII}$ [cm$^{-2}$]',fontsize=16)
ax_NOVII_all.legend(loc="upper right", ncol=1,fontsize=10)
fig_NOVII_all.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_NOVII_all.savefig('NOVII_all.' + lsfilename + '.png')

ax_NOVIII_all.set_xlim(0,300)
ax_NOVIII_all.set_ylim(12,18)
ax_NOVIII_all.set_xlabel('$b$ [kpc]',fontsize=16)
ax_NOVIII_all.set_ylabel('log $N_{OVIII}$ [cm$^{-2}$]',fontsize=16)
ax_NOVIII_all.legend(loc="upper right", ncol=1,fontsize=10)
fig_NOVIII_all.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_NOVIII_all.savefig('NOVIII_all.' + lsfilename + '.png')

ax_fhot_halo.set_xlim(11.7,15.3)
ax_fhot_halo.set_ylim(0,0.17)
#ax_fhot_halo.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_fhot_halo.set_xlabel('log $M_{200}$ [$M_{\odot}$]',fontsize=16)
ax_fhot_halo.set_ylabel('$f_\mathrm{hot}$',fontsize=16)
first_legend = ax_fhot_halo.legend(loc="upper left", ncol=1,fontsize=12)
second_legend = ax_fhot_halo.legend((fhot_200_label[0],fhot_500_label[0]),("$M_\mathrm{hot,200}/M_{200}$",r"$M_\mathrm{hot,500}/M_{500}$"),loc="lower right", ncol=1,fontsize=12)
fig_fhot_halo.gca().add_artist(first_legend)
fig_fhot_halo.gca().add_artist(second_legend)
#ax_fhot_halo.plot([1,1],[-10,10],color='k',lw=1,ls=':')
#ax_fhot_halo.plot([-10,10],[0,0],color='k',lw=1,ls=':')
#ax_fhot_halo.legend(loc="upper right", ncol=1,fontsize=10)
fig_fhot_halo.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_fhot_halo.savefig('fhot_halo.' + lsfilename + '.png')
