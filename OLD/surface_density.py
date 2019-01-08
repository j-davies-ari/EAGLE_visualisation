import numpy as np
import h5py as h5
from sys import exit
from python_tools import *
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from multiprocessing import Pool
from scipy import stats
#from multiprocessing import Manager
#from multiprocessing import Array
#import ctypes as c

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
res = 256
sim = 'L0025N0752'
run = 'REFERENCE'
frame_size = 70 # In kpc
pixel_size = 2.*np.float(frame_size)/np.float(res)
dz = 10.
ap_enclosed = 300
do_stars = False
show_gals = False
##############################################################################################################
ct = 0
sdensities = []
sfr_sdensities = []
av_metals = []
for n in np.arange(10,45):
    ct += 1
    '''
    if n != 31:
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
    gas_metals = np.array(data['Gas/Metallicity'])
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
    gas_metals = gas_metals[mask1]

    mask2 = (abs(star_coords[:,0])<=frame_size)&(abs(star_coords[:,1])<=frame_size)&(abs(star_coords[:,2])<=dz)
    star_coords = star_coords[mask2]
    star_h = star_h[mask2]
    star_mass = star_mass[mask2]
    star_density = star_density[mask2]

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
    
    
    gasmap = np.zeros((res,res))
    sfrmap = np.zeros((res,res))
    metalmap = np.zeros((res,res))
    starmap = np.zeros((res,res))
    
    '''
    manager = Manager()
    gasmap_sph = manager.Array('f', res*res)
    sfrmap = manager.Array('f', np.zeros((res,res)))
    metalmap = manager.Array('f', np.zeros((res,res)))
    
    
    def create_gasmap(gas_params):
        gas_xs,gas_ys,gas_zs,gas_mass,gas_h,gas_metals,sfr = gas_params
        kern = twod_kernel([gas_xs,gas_ys,gas_zs],gxys,gas_h)
        gasmap_sph = gas_mass*kern
        return gasmap_sph

    def create_sfrmap(gas_params):
        gas_xs,gas_ys,gas_zs,gas_mass,gas_h,gas_metals,sfr = gas_params
        kern = twod_kernel([gas_xs,gas_ys,gas_zs],gxys,gas_h)
        sfrmap = sfr*kern
        return sfrmap

    def create_metalmap(gas_params):
        gas_xs,gas_ys,gas_zs,gas_mass,gas_h,gas_metals,sfr = gas_params
        kern = twod_kernel([gas_xs,gas_ys,gas_zs],gxys,gas_h)
        metalmap = gas_metals*gas_mass*kern
        return metalmap

    def create_starmap(star_params):
        star_xs,star_ys,star_zs,star_mass,star_h = star_params
        kern = twod_kernel([star_xs,star_ys,star_zs],gxys,star_h)
        starmap_sph = star_mass*kern
        return starmap_sph

    def parallelise(f,params,num_cores):
        pool = Pool(num_cores)
        result = pool.map(f,params)
        #cleaned = [x for x in result if not x is None]
        cleaned = np.asarray(result)
        pool.close()
        pool.join()
        return cleaned

    num_cores = 16

    if __name__ == '__main__':  
        gas_params = zip(gas_coords[:,0],gas_coords[:,1],gas_coords[:,2],gas_mass,gas_h,gas_metals,sfr)
        print 'Creating a gas map using ',num_cores,' cores...'
        indv_maps = parallelise(create_gasmap,gas_params,num_cores)
        gasmap = np.sum(indv_maps,axis=0)
        del indv_maps
        print 'Creating an SFR map using ',num_cores,' cores...'
        indv_maps = parallelise(create_sfrmap,gas_params,num_cores)
        sfrmap = np.sum(indv_maps,axis=0)
        del indv_maps
        print 'Creating a metallicity map using ',num_cores,' cores...'
        indv_maps = parallelise(create_metalmap,gas_params,num_cores)
        metalmap = np.sum(indv_maps,axis=0)
        del indv_maps
        print 'Creating a star map using ',num_cores,' cores...'
        star_params = zip(star_coords[:,0],star_coords[:,1],star_coords[:,2],star_mass,star_h)
        indv_maps = parallelise(create_starmap,star_params,num_cores)
        starmap = np.sum(indv_maps,axis=0)
        del indv_maps


    '''
    print 'Creating gas-based maps...'
    for p in tqdm(range(no_gas)):
        kern = twod_kernel(gas_coords[p,:],gxys,gas_h[p]) #######
        if np.sum(kern)<1e-10:
            continue
        gasmap += gas_mass[p]*kern
        sfrmap += sfr[p]*kern
        metalmap += gas_metals[p]*gas_mass[p]*kern
    
    if do_stars == True:
        print 'Creating star-based maps...'
        for p in tqdm(range(no_star)):
            star_kern = twod_kernel(star_coords[p,:],gxys,star_h[p])
            if np.sum(star_kern)<1e-10:
                continue
            starmap += star_mass[p]*star_kern
    else:
        starmap += np.ones((res,res))
    
    #gasmap_sph *= 2*dz         #############
    #sfrmap *= 2*dz

    gmassmap = gasmap*pixel_size**2
    rawsfrmap = sfrmap*pixel_size**2
    starmassmap = starmap*pixel_size**2

    gasmap /= 1e6 # Convert to Msun pc-2
    starmap /= 1e6 # Convert to Msun pc-2
    metalmap /= 1e6 # Convert to Msun pc-2

    gasmap[gasmap<1e-4]=1e-3
    sfrmap[sfrmap<1e-7]=1e-7
    metalmap[metalmap<1e-7]=1e-7

    starmap[starmap<1e-1]=1e-1
    image_total = np.sum(gmassmap)
    imagesfr_total = np.sum(rawsfrmap)
    imagestar_total = np.sum(starmassmap)
    
    log_gasmap = np.log10(gasmap)
    log_sfrmap = np.log10(sfrmap)
    log_metalmap = np.log10(metalmap)
    log_starmap = np.log10(starmap)    

    # Normalisation test
    print 'Total mass of gas particles: ',np.sum(gas_mass)
    print 'Total mass smoothed onto image: ',image_total
    print 'Total SFR of gas particles: ',np.sum(sfr)
    print 'Total SFR smoothed onto image: ',imagesfr_total
    print 'Total mass of star particles: ',np.sum(star_mass)
    print 'Total stellar mass smoothed onto image: ',imagestar_total

    # Try a surface density profile
    profile1 = np.sum(gasmap,axis=1)/np.float(res)
    profile2 = np.sum(gasmap,axis=0)/np.float(res)
    
    '''
    radii = np.arange(1.,50.,1.)
    for ii,rad in enumerate(radii):
        circle,npix = circle_kernel([0.,0.,0.],gxys,rad)
        masked_mass = gmassmap*circle
        masked_sfr = sfrmap*circle
        sdensities.append(np.sum(masked_mass)/(npix*pixel_size**2))
        sfr_sdensities.append(np.sum(masked_sfr)/(npix*pixel_size**2))
    '''
    
    nos,xe,ye = np.histogram2d(gas_coords[:,0],gas_coords[:,1],[xedges,yedges])

    # Generate random aperture centres
    # Find surface values as soon as the aperture expands to enclose 1e7 Msun mass
    print 'Calculating surface densities...'
    for ii in tqdm(range(100)):
        xcoord = random.uniform(-frame_size,frame_size)
        ycoord = random.uniform(-frame_size,frame_size)
        #print 'Centre at ',xcoord,',',ycoord
        ri = 0.1
        count = 0
        while True:
            count +=1
            if frame_size-ri < abs(xcoord):
                break
            elif frame_size-ri < abs(ycoord):
                break
            circle,npix = circle_kernel([xcoord,ycoord,0.],gxys,ri)
            masked_number = nos*circle
            if np.sum(masked_number)<ap_enclosed:
                ri +=0.1
                continue
            else:
                masked_mass = gmassmap*circle
                masked_sfr = sfrmap*circle
                masked_metals = metalmap*circle
                sdensities.append(np.sum(masked_mass)/(npix*pixel_size**2))
                sfr_sdensities.append(np.sum(masked_sfr)/(npix*pixel_size**2))
                av_metals.append(np.sum(masked_metals)/npix)
                #print 'Final radius: ',ri,' kpc'
                break
        #print count,' iterations'
    
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
    
    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_gasmap,cmap='hot',extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$')
    plt.ylabel(r'$y$ $\mathrm{kpc}$')
    col = plt.colorbar()
    col.set_label(r'$\log\Sigma_g$ ($M_{\odot}$ $\mathrm{pc}^{-2})$')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_gas.png')
    if show_gals:
        plt.show()

    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_sfrmap,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$')
    plt.ylabel(r'$y$ $\mathrm{kpc}$')
    col = plt.colorbar()
    col.set_label(r'$\log\dot{\Sigma}_*$ ($M_{\odot}$ $\mathrm{yr}^{-1}$ $\mathrm{kpc}^{-2})$')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_sfr.png')
    if show_gals:
        plt.show()

    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_metalmap,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$')
    plt.ylabel(r'$y$ $\mathrm{kpc}$')
    col = plt.colorbar()
    col.set_label(r'$\log\Sigma_Z$ ($M_{\odot}$ $\mathrm{pc}^{-2})$')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_metals.png')
    if show_gals:
        plt.show()
    
    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_starmap,cmap='bone',extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$')
    plt.ylabel(r'$y$ $\mathrm{kpc}$')
    col = plt.colorbar()
    col.set_label(r'$\log\Sigma_*$ ($M_{\odot}$ $\mathrm{pc}^{-2})$')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_star.png')
    if show_gals:
        plt.show()


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

slope, intercept, r_value, p_value, std_err = stats.linregress(log_sdens[log_sfrdens>-6],log_sfrdens[log_sfrdens>-6])

bestfit_y = slope*np.log10(KS_x) + intercept

fig = plt.figure()
plt.scatter(log_sdens,log_sfrdens,c=np.log10(av_metals))
col = plt.colorbar()
col.set_label(r'$\Sigma_Z$ ($M_{\odot}$ $\mathrm{pc}^{-2})$')
plt.plot(np.log10(KS_x),KS_y,ls='--',c='r')
plt.plot(np.log10(KS_x),bestfit_y,ls='--',c='b')
plt.ylabel(r'$\dot{\Sigma}_*$ ($M_{\odot}$ $\mathrm{yr}^{-1}$ $\mathrm{kpc}^{-2})$')
plt.xlabel(r'$\Sigma_g$ ($M_{\odot}$ $\mathrm{pc}^{-2})$')
plt.xlim(-2,2)
plt.ylim(-6.5,1)
plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/KS_TEST_apN'+str(ap_enclosed)+'_dz'+str(dz)+'.png')
fig.show()


