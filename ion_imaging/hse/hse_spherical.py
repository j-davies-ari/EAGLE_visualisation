import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24

#def calc_Euler_radial(rbins, totmass, P, density, vrad, vrad_mag, vtan, vtan_mag, accel, galidstr, label, M200, R200):
def calc_Euler_radial(rbins, Totmass_R, P, density, vradial, vtan, sigmaradial, accrad, galidstr, label, M200, R200):

    nr = len(rbins)-1

    rbins_c = np.zeros(len(rbins)-1)  # center of rbins
    rbins_M = np.zeros(len(rbins)-1)  # center of rbins
    for i in range(nr):
        rbins_c[i] = (rbins[i]+rbins[i+1])/2.
        rbins_M[i] = rbins[i+1]


    nr = len(rbins)-1
    Mtot_R = np.zeros(nr)
    Mtherm_R = np.zeros(nr)
    Mrot_R = np.zeros(nr)
    Mstream_R = np.zeros(nr)
    Macc_R = np.zeros(nr)
    Msum_R = np.zeros(nr)
    Mcurl_R = np.zeros(nr)
    Mraddisp_R = np.zeros(nr)
    Mrand_R = np.zeros(nr)

    fhse = open("hse_radial_%s.%s.dat"%(galidstr,label),"w")
    fhse.write("#rbins Mtot Mtherm Mrot Mstream Macc Msum Mcurl Mraddisp\n")
                       
    for i in range(nr-1):
        if(i>0):
            Mtot_R[i] = (Totmass_R[i-1]+Totmass_R[i])/2.
        else:
            Mtot_R[i] = Totmass_R[i]            
        if((i==0) | (i == nr-1)): continue
        if(np.isnan(P[i])): continue
        SurfaceArea = rbins_c[i]**2* 4*np.pi
        #SurfaceCovered[i] += 4*np.pi
        M_prefix = 1/(4*np.pi*G_Grav)*SurfaceArea
        Mtherm_R[i] = -M_prefix*(K_Boltz/density[i])*(0.5*(P[i-1]-P[i])/(rbins_c[i-1]-rbins_c[i])+0.5*(P[i]-P[i+1])/(rbins_c[i]-rbins_c[i+1]))
        
        Mrand_R[i] = -M_prefix*(1/density[i])*(0.5*(density[i-1]*sigmaradial[i-1]**2-density[i]*sigmaradial[i]**2)/(rbins_c[i-1]-rbins_c[i])+0.5*(density[i]*sigmaradial[i]**2-density[i+1]*sigmaradial[i+1]**2)/(rbins_c[i]-rbins_c[i+1]))

        Mrot_R[i] = M_prefix*(vtan[i]**2)/(rbins_c[i])

        Macc_R[i] = -M_prefix*(accrad[i])        
        if(M200<13.8): #if((accrad[0] != 0) & (accrad[nr-1] != 0)): # BDO 12/4/17
            print "Subtracting Centripetal from Acceleration (Radial)"
            Macc_R[i] -= Mrot_R[i]  # BDO 11/10/17

        Mstream_radial = -M_prefix*(vradial[i])*(0.5*(vradial[i-1]-vradial[i])/(rbins_c[i-1]-rbins_c[i])+0.5*(vradial[i]-vradial[i+1])/(rbins_c[i]-rbins_c[i+1]))

        Mstream_R[i] += Mstream_radial

        Msum_R[i] += (Mtherm_R[i] + Mrot_R[i] + Macc_R[i] + Mstream_R[i])# + Mrand_R[i])
                #if((j%5==0) & (k%5==0)): print "SurfaceArea = ", SurfaceArea, " Mtherm_3D[i] = ", Mtherm_3D[i], " i,j,k= ", i, j, k, "Mstream = ", Mstream_radial, Mstream_theta, Mstream_phi #P_3D[i], density_3D[i],rbins[i], M_prefix, (K_Boltz/density_3D[i]), (P_3D[i-1]-P_3D[i])/(rbins[i-1]-rbins[i]), (P_3D[i]-P_3D[i+1])/(rbins[i]-rbins[i+1])
                
        print "EULER_RADIAL= ", i, " Ratio= ", Msum_R[i]/Mtot_R[i], "Mtot= ", Mtot_R[i]/M_Solar, "Mtherm= ", Mtherm_R[i]/M_Solar, "Mrot= ", Mrot_R[i]/M_Solar, "Mstream= ", Mstream_R[i]/M_Solar, "Macc= ", Macc_R[i]/M_Solar
        fhse.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e\n"%(np.log10(rbins_c[i]/cmpermpc), Mtot_R[i]/M_Solar, Mtherm_R[i]/M_Solar, Mrot_R[i]/M_Solar, Mstream_R[i]/M_Solar, Macc_R[i]/M_Solar, Msum_R[i]/M_Solar))

    Euler_mass_radial_plot(rbins/cmpermpc, Mtot_R/M_Solar, Mtherm_R/M_Solar, Mrot_R/M_Solar, Mstream_R/M_Solar, Macc_R/M_Solar, M200, R200, "hse_radialsum_%s.%s"%(galidstr,label),label,Mrand_R/M_Solar)

    fhse.close()

    return (Mtot_R, Mtherm_R)


    

def calc_Euler_spherical(rbins, thetabins, phibins, Totmass_R, P_3D, density_3D, vradial_3D, vtheta_3D, vphi_3D, sigmaradial_3D, accrad_3D, galidstr, label, M200, R200):

    nr = len(rbins)-1
    ntheta = len(thetabins)-1
    nphi = len(phibins)-1
    print "rbins = ", rbins
    print "thetabins = ", thetabins
    print "phibins = ", phibins

    print "LENGTH of rbins, Totmass_R = ", len(rbins), len(Totmass_R)
    print "rbins= ", rbins
    print "Totmass_R= ", Totmass_R
    rbins_c = np.zeros(len(rbins)-1)  # center of rbins
    rbins_M = np.zeros(len(rbins)-1)  # center of rbins
    for i in range(nr):
        rbins_c[i] = (rbins[i]+rbins[i+1])/2.
        rbins_M[i] = rbins[i+1]
        print "LENGTH of thetabins, P_3D-shape, ntheta = ", len(thetabins), P_3D.shape, ntheta

    thetabins_c = np.zeros(ntheta) 
    for j in range(ntheta):
        thetabins_c[j] = (thetabins[j]+thetabins[j+1])/2.
        
        # rbins[i] -- 0.5 -- rbins_c[i] -- 0.5 -- rbins[i+1]
        
    print "nr,ntheta, nphi = ",nr, ntheta, nphi 
    SurfaceCovered = np.zeros(nr)
    SolidAngle = np.zeros((ntheta,nphi))
    Mtot_R = np.zeros(nr)
    Mtherm_R = np.zeros(nr)
    Mrot_R = np.zeros(nr)
    Mstream_R = np.zeros(nr)
    Macc_R = np.zeros(nr)
    Msum_R = np.zeros(nr)
    Mtherm_3D = np.zeros((nr, ntheta, nphi))
    Mrot_3D = np.zeros((nr, ntheta, nphi))
    Mstream_3D = np.zeros((nr, ntheta, nphi))
    Macc_3D = np.zeros((nr, ntheta, nphi))
    Msum_3D = np.zeros((nr, ntheta, nphi))    
    Mrand_R = np.zeros(nr)
    Mrand_3D = np.zeros((nr, ntheta, nphi))

    fhse = open("hse_sphere_%s.%s.dat"%(galidstr,label),"w")
    fhse.write("#rbins Mtot Mtherm Mrot Mstream Macc Msum\n")

    
    for i in range(nr-1):
        if(i>0):
            Mtot_R[i] = (Totmass_R[i-1]+Totmass_R[i])/2.
        else:
            Mtot_R[i] = Totmass_R[i]            
        for j in range(ntheta):
            for k in range(nphi):
                if((i==0) | (i == nr-1)): continue
                ###if((j==0) | (j == ntheta-1)): continue # not doing poles.  # No we want to do poles alll except for Mstream theta.  
                if(np.isnan(P_3D[i,j,k])): continue
                kleft = k-1
                kright = k+1
                if(kleft<0): kleft = nphi-1 # Wrap around
                if(kright>nphi-1): kright = 0 # Wrap around
                if(np.isnan(P_3D[i-1,j,k]) | np.isnan(P_3D[i+1,j,k])): continue
                if((P_3D[i-1,j,k] < 0) | (P_3D[i,j,k] < 0) | (P_3D[i+1,j,k] < 0)): continue
                SurfaceArea = rbins_c[i]**2*(thetabins[j+1]-thetabins[j])*(phibins[1]-phibins[0])*np.sin(thetabins_c[j])
                SurfaceCovered[i] += (thetabins[j+1]-thetabins[j])*(phibins[1]-phibins[0])*np.sin(thetabins_c[j]) #/(4*np.pi)
                SolidAngle[j,k] = (thetabins[j+1]-thetabins[j])*(phibins[1]-phibins[0])*np.sin(thetabins_c[j])/(4*np.pi)
                M_prefix = 1/(4*np.pi*G_Grav)*SurfaceArea
                Mtherm_3D[i,j,k] = -M_prefix*(K_Boltz/density_3D[i,j,k])*(0.5*(P_3D[i-1,j,k]-P_3D[i,j,k])/(rbins_c[i-1]-rbins_c[i])+0.5*(P_3D[i,j,k]-P_3D[i+1,j,k])/(rbins_c[i]-rbins_c[i+1]))

                #Mtherm_3D[i,j,k] = -M_prefix*(K_Boltz/density_3D[i,j,k])*((P_3D[i,j,k]-P_3D[i+1,j,k])/(rbins_C[i]-rbins_C[i+1]))

                ###Mtherm_3D[i,j,k] = -M_prefix*K_Boltz*((P_3D[i-1,j,k]/density_3D[i-1,j,k]-P_3D[i,j,k]/density_3D[i,j,k])/(rbins_c[i-1]-rbins_c[i])+(P_3D[i,j,k]/density_3D[i,j,k]-P_3D[i+1,j,k]/density_3D[i+1,j,k])/(rbins_c[i]-rbins_c[i+1])) # This is way wrong, because it is the partial with respect to T only!                  
                Mtherm_R[i] += Mtherm_3D[i,j,k]
                #if (i==20):
                    #print "I20= ", Mtherm_3D[i,j,k],i,j,k, "sum= ", Mtherm_R[i] 

                Mrand_3D[i,j,k] = -M_prefix*(1/density_3D[i,j,k])*(0.5*(density_3D[i-1,j,k]*sigmaradial_3D[i-1,j,k]**2-density_3D[i,j,k]*sigmaradial_3D[i,j,k]**2)/(rbins_c[i-1]-rbins_c[i])+0.5*(density_3D[i,j,k]*sigmaradial_3D[i,j,k]**2-density_3D[i+1,j,k]*sigmaradial_3D[i+1,j,k]**2)/(rbins_c[i]-rbins_c[i+1]))
                Mrand_R[i] += Mrand_3D[i,j,k]

                Mrot_3D[i,j,k] = M_prefix*(vtheta_3D[i,j,k]**2+vphi_3D[i,j,k]**2)/(rbins_c[i])
                Mrot_R[i] += Mrot_3D[i,j,k]

                Macc_3D[i,j,k] = -M_prefix*(accrad_3D[i,j,k])
                if(M200<13.8): #if((accrad_3D[0,0,0] != 0) & (accrad_3D[nr-1,0,0] != 0)): # BDO 12/4/17
                    print "Subtracting Centripetal from Acceleration (Spherical)"
                    Macc_3D[i,j,k] -= Mrot_3D[i,j,k] # BDO 11/10/17

                Macc_R[i] += Macc_3D[i,j,k]                

                Mstream_radial = -M_prefix*(vradial_3D[i,j,k])*(0.5*(vradial_3D[i-1,j,k]-vradial_3D[i,j,k])/(rbins_c[i-1]-rbins_c[i])+0.5*(vradial_3D[i,j,k]-vradial_3D[i+1,j,k])/(rbins_c[i]-rbins_c[i+1]))

                if(np.isnan(P_3D[i,j,kleft]) | np.isnan(P_3D[i,j,kright])):
                    Mstream_phi = 0
                else:
                    Mstream_phi = -M_prefix*(vphi_3D[i,j,k]/rbins_c[i]/np.sin(thetabins_c[j])) * (0.5*(vradial_3D[i,j,kleft]-vradial_3D[i,j,k])/(phibins[kleft]-phibins[k])+0.5*(vradial_3D[i,j,k]-vradial_3D[i,j,kright])/(phibins[k]-phibins[kright]))

                Mstream_theta = 0
                if(j>0):
                    if not np.isnan(P_3D[i,j-1,k]):
                        Mstream_theta += -M_prefix*(vtheta_3D[i,j,k]/rbins_c[i]) * (0.5*(vradial_3D[i,j-1,k]-vradial_3D[i,j,k])/(thetabins[j-1]-thetabins[j]))
                if(j<ntheta-1):
                    if not np.isnan(P_3D[i,j+1,k]):
                        Mstream_theta += -M_prefix*(vtheta_3D[i,j,k]/rbins_c[i]) * (0.5*(vradial_3D[i,j,k]-vradial_3D[i,j+1,k])/(thetabins[j]-thetabins[j+1]))

                Mstream_3D[i,j,k] = Mstream_radial + Mstream_theta + Mstream_phi
                if((np.isnan(Mstream_3D[i,j,k])) | (np.isinf(Mstream_3D[i,j,k]))):
                    Mstream_3D[i,j,k] = 0.0

                #Mstream_3D += -M_prefix*(vphi_3D[i,j,k]/rbins[i]) // /simtheta dvrad/dphi
                Mstream_R[i] += Mstream_3D[i,j,k]

                Msum_R[i] += (Mtherm_3D[i,j,k] + Mrot_3D[i,j,k] + Macc_3D[i,j,k] + Mstream_3D[i,j,k])# + Mrand_3D[i,j,k])
                #if((j%5==0) & (k%5==0)): print "SurfaceArea = ", SurfaceArea, " Mtherm_3D[i,j,k] = ", Mtherm_3D[i,j,k], " i,j,k= ", i, j, k, "Mstream = ", Mstream_radial, Mstream_theta, Mstream_phi #P_3D[i,j,k], density_3D[i,j,k],rbins[i], M_prefix, (K_Boltz/density_3D[i,j,k]), (P_3D[i-1,j,k]-P_3D[i,j,k])/(rbins[i-1]-rbins[i]), (P_3D[i,j,k]-P_3D[i+1,j,k])/(rbins[i]-rbins[i+1])
                
        Mtherm_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        Mrand_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        Mstream_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        Macc_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        Mrot_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        Msum_R[i] /= (SurfaceCovered[i]/(4*np.pi))
        print i, " Ratio= ", Msum_R[i]/Mtot_R[i], "Mtot= ", Mtot_R[i]/M_Solar, "Mtherm= ", Mtherm_R[i]/M_Solar, "Mrot= ", Mrot_R[i]/M_Solar, "Mstream= ", Mstream_R[i]/M_Solar, "Macc= ", Macc_R[i]/M_Solar, " SurfaceCovered = ", SurfaceCovered[i]/(4*np.pi)
        fhse.write("% 5.3f % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e % 5.3e %5.3f %5.3e\n"%(np.log10(rbins_c[i]/cmpermpc), Mtot_R[i]/M_Solar, Mtherm_R[i]/M_Solar, Mrot_R[i]/M_Solar, Mstream_R[i]/M_Solar, Macc_R[i]/M_Solar, Msum_R[i]/M_Solar, SurfaceCovered[i]/(4*np.pi),Mrand_R[i]/M_Solar))


    Euler_mass_radial_plot(rbins/cmpermpc, Mtot_R/M_Solar, Mtherm_R/M_Solar, Mrot_R/M_Solar, Mstream_R/M_Solar, Macc_R/M_Solar, M200, R200, "hse_sphericalsum_%s.%s"%(galidstr,label),label, Mrand_R/M_Solar)

    Euler_mass_spherical_plot(rbins_c/cmpermpc, thetabins, phibins, Mtot_R/M_Solar, Mtherm_3D/M_Solar, Mrot_3D/M_Solar, Mstream_3D/M_Solar, Macc_3D/M_Solar, SolidAngle, M200, R200, "hse_sphericalmap_%s.%s"%(galidstr,label),label)

    fhse.close()

    return (Mtot_R, Mtherm_R)

def Euler_mass_spherical_plot(rbins, thetabins, phibins, Mtot, Mtherm, Mrot, Mstream, Macc, SolidAngle, M200, R200, figname_base,label):

    nr = len(rbins)-1
    ntheta = len(thetabins)-1
    nphi = len(phibins)-1
    phi_otherside = nphi/2 #-1
    theta_top = ntheta
    radial_max = nr-7
    
    Msum = Mtherm+Mrot+Mstream+Macc
    rbins_M = np.zeros(len(rbins)-1)

    for i in xrange(len(rbins_M)):
        rbins_M[i] = (rbins[i]+rbins[i+1])/2. ### XXX- BDO Check, because it seems it may be rbins[i].  

    fig_sum_hse = plt.figure(figsize=(5.0,4.5))
    ax_sum_hse = fig_sum_hse.add_subplot(111, polar=True)

    fig_therm_hse = plt.figure(figsize=(5.0,4.5))
    ax_therm_hse = fig_therm_hse.add_subplot(111, polar=True)

    fig_rot_hse = plt.figure(figsize=(5.0,4.5))
    ax_rot_hse = fig_rot_hse.add_subplot(111, polar=True)

    fig_acc_hse = plt.figure(figsize=(5.0,4.5))
    ax_acc_hse = fig_acc_hse.add_subplot(111, polar=True)
    
    #theta_ = np.linspace(0+2*np.pi/(4*ntheta),2*np.pi-2*np.pi/(4*ntheta),ntheta*2)
    theta_ = np.linspace(0,2*np.pi,ntheta*2+1)
    r_ = rbins_M[0:radial_max]/R200
    pol_sum, B = np.meshgrid(theta_,r_)
    pol_therm, B = np.meshgrid(theta_,r_)
    pol_rot, B = np.meshgrid(theta_,r_)
    pol_stream, B = np.meshgrid(theta_,r_)
    pol_acc, B = np.meshgrid(theta_,r_)
 
    for i in xrange(radial_max-1):
        for j in xrange(ntheta):
            pol_sum[i,j] = Msum[i,j,0]/Mtot[i]/SolidAngle[j,0]
            pol_sum[i,theta_top+j] = Msum[i,theta_top-1-j,phi_otherside]/Mtot[i]/SolidAngle[j,phi_otherside]
            if(i==40):
                print "POL_SUM= ", pol_sum[i,j], Msum[i,j,0], Mtot[i], SolidAngle[j,0], j

                
            pol_therm[i,j] = Mtherm[i,j,0]/Mtot[i]/SolidAngle[j,0]
            pol_therm[i,theta_top+j] = Mtherm[i,theta_top-1-j,phi_otherside]/Mtot[i]/SolidAngle[j,phi_otherside]
            pol_rot[i,j] = Mrot[i,j,0]/Mtot[i]/SolidAngle[j,0]
            pol_rot[i,theta_top+j] = Mrot[i,theta_top-1-j,phi_otherside]/Mtot[i]/SolidAngle[j,phi_otherside]
            pol_acc[i,j] = Macc[i,j,0]/Mtot[i]/SolidAngle[j,0]
            pol_acc[i,theta_top+j] = Macc[i,theta_top-1-j,phi_otherside]/Mtot[i]/SolidAngle[j,phi_otherside]
            if(j==0):
                pol_sum[i,2*theta_top] = Msum[i,j,0]/Mtot[i]/SolidAngle[j,0]
                pol_therm[i,2*theta_top] = Mtherm[i,j,0]/Mtot[i]/SolidAngle[j,0]
                pol_rot[i,2*theta_top] = Mrot[i,j,0]/Mtot[i]/SolidAngle[j,0]
                pol_acc[i,2*theta_top] = Macc[i,j,0]/Mtot[i]/SolidAngle[j,0]


    levels_sum = [0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
    polsum = ax_sum_hse.contourf(theta_,r_,pol_sum, levels_sum, cmap='RdYlBu_r', extend='both')
    #ax_sum_hse.set_xticklabels([])
    #ax_sum_hse.set_yticklabels([])
    ax_sum_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    cb = fig_sum_hse.colorbar(polsum)
    cb.set_label("M$_{\mathrm{Sum}}/$M$_{\mathrm{Tot}}$")
    fig_sum_hse.savefig(figname_base + ".polsum.png")                

    levels_therm = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
    poltherm = ax_therm_hse.contourf(theta_,r_,pol_therm, levels_therm, cmap='RdYlBu_r',extend='both')    
    #ax_therm_hse.set_xticklabels([])
    #ax_therm_hse.set_yticklabels([])
    ax_therm_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    cb = fig_therm_hse.colorbar(poltherm)  
    cb.set_label("M$_{\mathrm{Therm}}/$M$_{\mathrm{Tot}}$")
    fig_therm_hse.savefig(figname_base + ".poltherm.png")                
    
    levels_rot = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
    polrot = ax_rot_hse.contourf(theta_,r_,pol_rot, levels_rot, cmap='RdYlBu_r',extend='both')    
    #ax_rot_hse.set_xticklabels([])
    #ax_rot_hse.set_yticklabels([])
    ax_rot_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    cb = fig_rot_hse.colorbar(polrot)  
    cb.set_label("M$_{\mathrm{Rot}}/$M$_{\mathrm{Tot}}$")
    fig_rot_hse.savefig(figname_base + ".polrot.png")                

    levels_acc = [-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5]
    polacc = ax_acc_hse.contourf(theta_,r_,pol_acc, levels_acc, cmap='RdYlBu_r',extend='both')    
    #ax_acc_hse.set_xticklabels([])
    #ax_acc_hse.set_yticklabels([])
    ax_acc_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    cb = fig_acc_hse.colorbar(polacc)  
    cb.set_label("M$_{\mathrm{Acc}}/$M$_{\mathrm{Tot}}$")
    fig_acc_hse.savefig(figname_base + ".polacc.png")                
    
    #a = np.linspace(0,2*np.pi,50)
    #b = np.linspace(0,1,50)
    #A, B = np.meshgrid(a, b)
    #c = np.random.random(A.shape)

    #actual plotting
    #import matplotlib.cm as cm
    #ax = plt.subplot(111, polar=True)
    #ax.set_yticklabels([])
    #ctf = ax.contourf(a, b, c, cmap=cm.jet)
    #plt.colorbar(ctf)

    
    
def Euler_mass_radial_plot(rbins, Mtot, Mtherm, Mrot, Mstream, Macc, M200, R200, figname_base,label,Mrand):

    R200 = R200*1e+03
    rbins = rbins*1e+03 # kpc as of 5/10/18
    Msum = Mtherm+Mrot+Mstream+Macc
    rbins_M = np.zeros(len(rbins)-1)
    
    for i in xrange(len(rbins_M)):
        rbins_M[i] = (rbins[i]+rbins[i+1])/2. ### XXX- BDO Check, because it seems it may be rbins[i].  

    
    fig_log_hse = plt.figure(figsize=(4.5,3.375))
    ax_log_hse = fig_log_hse.add_subplot(111)

    ax_log_hse.plot(rbins_M, np.log10(Mtot), color = 'black', lw=4,label="$M_{\mathrm{tot}}$",zorder=-10)
    ax_log_hse.plot(rbins_M, np.log10(Mtherm), color = 'red', lw=1,label="$M_{\mathrm{therm}}$",zorder=5)
    ax_log_hse.plot(rbins_M, np.log10(-Mtherm), color = 'red', lw=1, ls=':',zorder=5)
    ax_log_hse.plot(rbins_M, np.log10(Mrot), color = 'lime', lw=1,label="$M_{\mathrm{rot}}$",zorder=5)
    if(M200<13.8): # BDO 12/4/17
        ax_log_hse.plot(rbins_M, np.log10(Macc), color = 'magenta', lw=1,label="$M_{\mathrm{acc}}$",zorder=5)
        ax_log_hse.plot(rbins_M, np.log10(-Macc), color = 'magenta', lw=1,ls=":",zorder=5)
    ax_log_hse.plot(rbins_M, np.log10(Mstream), color = 'blue', lw=1,label="$M_{\mathrm{stream}}$",zorder=5)
    ax_log_hse.plot(rbins_M, np.log10(-Mstream), color = 'blue', lw=1,ls=":",zorder=5)
    ax_log_hse.plot(rbins_M, np.log10(Msum), color = 'orange', lw=2, label="$M_{\mathrm{sum}}$",zorder=0)

    ax_log_hse.plot([-1000,10000], [M200,M200], color = 'black', lw=2,ls=':')
    print "R200 hse_spherical = ", R200
    ax_log_hse.plot([R200,R200], [-10,20], color = 'black', lw=2, ls=':')

    ax_log_hse.set_xscale('log')
    #ax_log_hse.set_xlim(-2.0,0.5)
    ax_log_hse.set_xlim(10,3000)
    ax_log_hse.set_xticks([10,20,50,100,200,500,1000,2000])
    ax_log_hse.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_log_hse.set_ylim(10.0,16.0)
    ax_log_hse.set_xlabel('R [kpc]', fontsize=16)
    ax_log_hse.set_ylabel('log M [$M_{\odot}$]', fontsize=16)
    ax_log_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    ax_log_hse.legend(loc="upper center", ncol=3,fontsize=11)

    fig_log_hse.subplots_adjust(left=0.17, bottom=0.17,top=0.90,right=0.97)
    fig_log_hse.savefig(figname_base + ".log.png")
                

    fig_lin_hse = plt.figure(figsize=(4.5,3.375))
    ax_lin_hse = fig_lin_hse.add_subplot(111)

    ax_lin_hse.plot(rbins_M, Mtot/Mtot, color = 'black', lw=4,zorder=-10)#,label="$M_{\mathrm{tot}}$")
    ax_lin_hse.plot(rbins_M, Mtherm/Mtot, color = 'red', lw=1,label="$\mathcal{S}_{\mathrm{therm}}$",zorder=5)
    ax_lin_hse.plot(rbins_M, Mrot/Mtot, color = 'lime', lw=1,label="$\mathcal{S}_{\mathrm{rot}}$",zorder=5)
    if(M200<13.8): # BDO 12/4/17
        ax_lin_hse.plot(rbins_M, Macc/Mtot, color = 'magenta', lw=1,label="$\mathcal{S}_{\mathrm{acc}}$",zorder=5)
    ax_lin_hse.plot(rbins_M, Mstream/Mtot, color = 'blue', lw=1,label="$\mathcal{S}_{\mathrm{stream}}$",zorder=5)
    #ax_lin_hse.plot(rbins_M, Mrand/Mtot, color = 'yellow', lw=1,label="$\mathcal{S}_{\mathrm{rand}}$",zorder=-10)
    ax_lin_hse.plot(rbins_M, Msum/Mtot, color = 'orange', lw=2, label="$\mathcal{S}_{\mathrm{sum}}$",zorder=0)

    ax_lin_hse.plot([-1000,10000], [0,0], color = 'black', lw=2,ls=':')

    ax_lin_hse.plot([R200,R200], [-10,20], color = 'black', lw=2, ls=':')

    ax_lin_hse.set_xscale('log')
    #ax_lin_hse.set_xlim(-2.0,0.5)
    ax_lin_hse.set_xlim(10,3000)
    ax_lin_hse.set_xticks([10,20,50,100,200,500,1000,2000])
    ax_lin_hse.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    #ax_lin_hse.set_xlim(-2.0,0.5)
    ax_lin_hse.set_ylim(-0.5,2.0)
    ax_lin_hse.set_xlabel('R [kpc]', fontsize=16)
    #ax_lin_hse.set_ylabel('$M/M_{\mathrm{tot}}$', fontsize=16)
    ax_lin_hse.set_ylabel('$\mathcal{S}_{term}$', fontsize=16)
    ax_lin_hse.set_title('log $M_{\mathrm{200}}=$ %s'%(M200),fontsize=16)
    ax_lin_hse.legend(loc="upper center", ncol=3,fontsize=12)

    fig_lin_hse.subplots_adjust(left=0.17, bottom=0.17,top=0.90,right=0.97)
    fig_lin_hse.savefig(figname_base + ".lin.png")                
