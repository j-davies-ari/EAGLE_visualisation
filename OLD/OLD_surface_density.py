import numpy as np
import h5py as h5
from sys import exit
from python_tools import *
import matplotlib.pyplot as plt
import random

def anarchy_kernel(particle_xyz,gridpoint_xys,h): # For a given particle, returns the SPH kernel on a 2D grid
    grid_size = int(np.sqrt(len(gridpoint_xys[:,0])))
    no_pixels = int(len(gridpoint_xys[:,0]))
    kernel = np.zeros(no_pixels)
    ones = np.reshape(np.ones(no_pixels),(no_pixels,1))
    particle_xys = ones*particle_xyz
    r = np.sqrt((abs(gridpoint_xys[:,0]-particle_xys[:,0]))**2 + (abs(gridpoint_xys[:,1]-particle_xys[:,1]))**2)
    mask = np.where(r<=h)
    kernel[mask] += (21/(2*np.pi*h**3)) * (1-r[mask]/h)**4 * (1+4*r[mask]/h)
    kernel_out = np.reshape(kernel,(grid_size,grid_size))
    return kernel_out

def twod_kernel(particle_xyz,gridpoint_xys,h): # For a given particle, returns the SPH kernel on a 2D grid
    grid_size = int(np.sqrt(len(gridpoint_xys[:,0])))
    no_pixels = int(len(gridpoint_xys[:,0]))
    kernel = np.zeros(no_pixels)
    ones = np.reshape(np.ones(no_pixels),(no_pixels,1))
    particle_xys = ones*particle_xyz
    r = np.sqrt((abs(gridpoint_xys[:,0]-particle_xys[:,0]))**2 + (abs(gridpoint_xys[:,1]-particle_xys[:,1]))**2)
    mask = np.where(r<=h)
    kernel[mask] += (7/(np.pi*h**2)) * (1-r[mask]/h)**4 * (1+4*r[mask]/h)
    kernel_out = np.reshape(kernel,(grid_size,grid_size))
    return kernel_out

def dimensionless_sph_weight(mass,density,h,n_dims):
    return mass/(density*(h**n_dims))

def circle_kernel(centre,gridpoint_xys,radius):
    no_pixels = int(len(gridpoint_xys[:,0]))
    grid_size = int(np.sqrt(np.float(no_pixels)))
    kernel = np.zeros(no_pixels)
    ones = np.reshape(np.ones(no_pixels),(no_pixels,1))
    centre_xys = ones*centre
    
    r = np.sqrt((abs(gridpoint_xys[:,0]-centre_xys[:,0]))**2 + (abs(gridpoint_xys[:,1]-centre_xys[:,1]))**2)
    mask = np.where(r<=radius)
    kernel[mask] += 1
    kernel_out = np.reshape(kernel,(grid_size,grid_size))
    return kernel_out,len(kernel[kernel==1])


kpc = 3.0857e19
msun = 1.989e30
h_0 = 0.6777

##############################################################################################################
# SET RESOLUTION
res = 255
sim = 'L0025N0752'
run = 'REFERENCE'
frame_size = 36 # In kpc
pixel_size = 2.*np.float(frame_size)/np.float(res)
dz = 5.
##############################################################################################################
ct = 0
sdensities = []
sfr_sdensities = []
for n in np.arange(5,50):
    ct += 1
    '''
    if n != 14:
        continue
    '''
    try:
	    data = h5.File('/home/arijdav1/Desktop/'+sim+'_'+run+'/028_z000p000/group'+str(n)+'.hdf5','r')
    except IOError:
	    continue

    print 'Group no. ',n

    #print data['Gas'].keys()

    star_coords = np.array(data['Stars/Coordinates'])
    star_h = np.array(data['Stars/SmoothingLength'])
    star_mass = np.array(data['Stars/Mass'])
    star_density = np.array(data['Stars/Density'])

    gas_coords = np.array(data['Gas/Coordinates'])
    gas_h = np.array(data['Gas/SmoothingLength'])
    gas_mass = np.array(data['Gas/Mass'])
    gas_density = np.array(data['Gas/Density'])
    sfr = np.array(data['Gas/StarFormationRate'])

    #print sfr

    '''
    #for testing
    gas_coords = np.array([[0.,-0.01,0.],[0.,0.01,0.]])
    gas_h = np.array([[0.006],[0.006]])
    gas_mass = np.array([[1.],[1.]])
    gas_density = np.array([[1.],[1.]])
    '''
    
    # Conversions from whatever is in the data
    star_coords *= 1e3
    gas_coords *= 1e3
    star_h *= 1e3
    gas_h *= 1e3

    star_mass *= 1e10
    gas_mass *= 1e10
    star_density *= 10
    gas_density *= 10
    

    '''
    # Convert from CGS into kpc and Msun units
    gas_coords /= 3.08568e20
    gas_h /= 3.08568e20
    gas_mass /= 1.989e33
    gas_density /= 6.770e-29
    '''

    # Set minimum smoothing length to half a pixel
    star_h[star_h<pixel_size/2.] = pixel_size/2.    
    gas_h[gas_h<pixel_size/2.] = pixel_size/2.
    
    # Remove all particles that aren't within h of the frame or within dz range
    mask1 = (abs(gas_coords[:,0])<=frame_size)&(abs(gas_coords[:,1])<=frame_size)&(abs(gas_coords[:,2])<=dz)
    gas_coords = gas_coords[mask1]
    gas_h = gas_h[mask1]
    gas_mass = gas_mass[mask1]
    gas_density = gas_density[mask1]
    sfr = sfr[mask1]

    mask2 = (abs(star_coords[:,0])<=frame_size)&(abs(star_coords[:,1])<=frame_size)
    star_coords = star_coords[mask2]

    no_star = len(star_coords[:,0])
    no_gas = len(gas_coords[:,0])
    xedges = np.linspace(-frame_size,frame_size,num=res+1)
    yedges = np.linspace(-frame_size,frame_size,num=res+1)
    bs = get_binsizes(xedges)[0]

    xcents = get_bincentres(xedges)
    ycents = get_bincentres(yedges)
    gridx, gridy = np.meshgrid(xcents, ycents)
    xlist = np.reshape(np.ravel(gridx),(len(np.ravel(gridx)),1))
    ylist = np.reshape(np.ravel(gridy),(len(np.ravel(gridy)),1))
    gxys = np.hstack((xlist,ylist))

    gasmap_sph = np.zeros((res,res))
    sfrmap = np.zeros((res,res))
    for p in range(no_gas):
        progress_bar(p,no_gas)
        kern = twod_kernel(gas_coords[p,:],gxys,gas_h[p]) #######

        #print np.sum(kern)

        if np.sum(kern)<1e-10:
            continue

        gasmap_sph += gas_mass[p]*kern/np.sum(kern)
        sfrmap += sfr[p]*kern/np.sum(kern)

        #gasmap_sph += gas_mass[p]*kern
        #sfrmap += sfr[p]*kern

    gasmap_sph /= pixel_size**2         #############

    gmassmap = gasmap_sph*pixel_size**2

    gasmap_sph[gasmap_sph<1e-4]=1e-3
    gasmap_sph /= 1e6 # Convert to Msun pc-2
    image_total = np.sum(gmassmap)
    #image_total = np.sum(gasmap_sph)
    log_gasmap = np.log10(gasmap_sph)

    # Normalisation test
    print 'Total mass of gas particles: ',np.sum(gas_mass)
    print 'Total mass smoothed onto image: ',image_total

    # Try a surface density profile
    profile1 = np.sum(gasmap_sph,axis=1)/np.float(res)
    profile2 = np.sum(gasmap_sph,axis=0)/np.float(res)
    
    '''
    radii = np.arange(1.,50.,1.)
    for ii,rad in enumerate(radii):
        circle,npix = circle_kernel([0.,0.,0.],gxys,rad)
        masked_mass = gmassmap*circle
        masked_sfr = sfrmap*circle
        sdensities.append(np.sum(masked_mass)/(npix*pixel_size**2))
        sfr_sdensities.append(np.sum(masked_sfr)/(npix*pixel_size**2))
    '''
    
    # Generate random aperture centres
    # Find surface values as soon as the aperture expands to enclose 1e7 Msun mass
    for ii in range(50):
        xcoord = random.uniform(-frame_size,frame_size)
        ycoord = random.uniform(-frame_size,frame_size)
        print 'Centre at ',xcoord,',',ycoord
        ri = 0.1
        count = 0
        while True:
            count +=1
            if frame_size-ri < abs(xcoord):
                break
            elif frame_size-ri < abs(ycoord):
                break
            circle,npix = circle_kernel([xcoord,ycoord,0.],gxys,ri)

            masked_mass = gmassmap*circle
            #print 'Enclosed mass: ',np.sum(masked_mass)
            if np.sum(masked_mass)<1e7:
                ri +=0.1
                continue
            else:
                masked_sfr = sfrmap*circle
                sdensities.append(np.sum(masked_mass)/(npix*pixel_size**2))
                sfr_sdensities.append(np.sum(masked_sfr)/(npix*pixel_size**2))
                break
        print count,' iterations'
    
    # Enclosed surface gas and sfr density as a fn of radius
    '''
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212,sharex=ax1)
    ax1.plot(radii,sdensities)
    ax2.plot(radii,sfr_sdensities)
    #ax1.set_xlim(0,90)
    #ax1.set_ylim(0,0.4)
    ax1.set_ylabel('Gas surface density (Msun pc-2)')
    ax2.set_ylabel('SFR surface density (Msun pc-2)')
    ax2.set_xlabel('Aperture radius (kpc)')
    plt.show()
    '''
    
    # Show the galaxy
    '''
    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_gasmap,cmap='hot',extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.colorbar()
    #plt.scatter(star_coords[:,0],star_coords[:,1],color='c',marker='*')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_gas.png')
    plt.show()
    '''


    # Attempt at a profile
    '''
    plt.figure(figsize=(10,6))
    plt.plot(xcents,np.log10((profile1+profile2)/2.),color='k')
    plt.plot(xcents,np.log10(profile1),color='b',ls='--')
    plt.plot(xcents,np.log10(profile2),color='r',ls='--')
    plt.xlabel('r')
    plt.ylabel('log stacked sigma_gas')
    plt.show()
    '''

sdensities = np.array(sdensities)/1e6 # Convert to Msun pc-2

KS_x = np.logspace(-2,3,num=100)
KS_y = 1.4*np.log10(KS_x) + np.log10(1.515e-4)

log_sdens = np.log10(sdensities)
log_sfrdens = np.log10(sfr_sdensities)

log_sfrdens[log_sfrdens<-6] = -6

plt.figure()
plt.scatter(log_sdens,log_sfrdens)
plt.plot(np.log10(KS_x),KS_y,ls='--',c='r')
plt.ylabel(r'$\dot{Sigma}_*$ ($M_{\odot}$ $\mathrm{yr}^{-1}$ $\mathrm{kpc}^{-2}$')
plt.xlabel(r'$Sigma_g$ ($M_{\odot}$ $\mathrm{pc}^{-2}$')
plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/KS_TEST.png')
plt.show()


