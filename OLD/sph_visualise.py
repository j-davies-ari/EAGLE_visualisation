import numpy as np
import h5py as h5
from sys import exit
from python_tools import *
import matplotlib.pyplot as plt

def dimensionless_anarchy_kernel(particle_xy,gridpoint_xys,h): # For a given particle, returns the SPH kernel on a 2D grid
    grid_size = int(np.sqrt(len(gridpoint_xys[:,0])))
    no_pixels = int(len(gridpoint_xys[:,0]))
    kernel = np.zeros(no_pixels)
    ones = np.reshape(np.ones(no_pixels),(no_pixels,1))
    particle_xys = ones*particle_xy
    r = np.sqrt((abs(gridpoint_xys[:,0]-particle_xys[:,0]))**2 + (abs(gridpoint_xys[:,1]-particle_xys[:,1]))**2)
    mask = np.where(r<=h)
    kernel[mask] += (21/(2*np.pi*h**2)) * (1-r[mask]/h)**4 * (1+4*r[mask]/h)
    kernel_out = np.reshape(kernel,(grid_size,grid_size))
    return kernel_out

def dimensionless_sph_weight(mass,density,h,n_dims):
    return mass/(density*(h**n_dims))

kpc = 3.0857e19
msun = 1.989e30
h_0 = 0.6777

##############################################################################################################
# SET RESOLUTION
res = 256
sim = 'L0025N0752'
run = 'REFERENCE'
frame_size = 36 # In kpc
pixel_size = 2.*np.float(frame_size)/np.float(res)
##############################################################################################################

for n in np.arange(5,100):
    
    if n != 14:
        continue
    
    try:
	    data = h5.File('/home/arijdav1/Desktop/'+sim+'_'+run+'/028_z000p000/group'+str(n)+'.hdf5','r')
    except IOError:
	    continue

    print 'Group no. ',n

    star_coords = np.array(data['Stars/Coordinates'])
    star_h = np.array(data['Stars/SmoothingLength'])
    star_mass = np.array(data['Stars/Mass'])
    star_density = np.array(data['Stars/Density'])

    gas_coords = np.array(data['Gas/Coordinates'])
    gas_h = np.array(data['Gas/SmoothingLength'])
    gas_mass = np.array(data['Gas/Mass'])
    gas_density = np.array(data['Gas/Density'])
    #sfr = np.array(data['Gas/StarFormationRate'])

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
    
    # Remove all particles that aren't within h of the frame
    mask1 = (abs(gas_coords[:,0])<=frame_size+gas_h)&(abs(gas_coords[:,1])<=frame_size+gas_h)
    gas_coords = gas_coords[mask1]
    gas_h = gas_h[mask1]
    gas_mass = gas_mass[mask1]
    gas_density = gas_density[mask1]

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
    for p in range(no_gas):
        progress_bar(p,no_gas)
        kern = dimensionless_anarchy_kernel(gas_coords[p,:2],gxys,gas_h[p]) 

        if np.sum(kern)<1e-10:
            continue

        gasmap_sph += gas_mass[p]*kern/np.sum(kern)
        #gasmap_sph += gas_mass[p]*kern

    gasmap_sph[gasmap_sph<1e-4]=1e-3
    image_total = np.sum(gasmap_sph)
    log_gasmap = np.log10(gasmap_sph)

    # Normalisation test
    print 'Total mass of gas particles: ',np.sum(gas_mass)
    print 'Total mass smoothed onto image: ',image_total

    # Rough surface density calculation
    chunk_sidelen = 20.
    px_req = chunk_sidelen/pixel_size
    inds = int(np.round(px_req/2.))
    start = res/2-inds
    end = res/2+inds
    chunk = gasmap_sph[start:end,start:end]
    surface_density = np.sum(chunk)/chunk_sidelen**2

    print 'Surface Density: ',surface_density,' Msun kpc^-2'



    '''
    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(starmap_sph,cmap='bone')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_stars.png')
    plt.show()
    '''    

    fig = plt.figure(figsize=(10,10))
    im1 = plt.imshow(log_gasmap,cmap='hot',extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.colorbar()
    #plt.scatter(star_coords[:,0],star_coords[:,1],color='c',marker='*')
    plt.savefig('/home/arijdav1/Desktop/figures/sph_pics/'+sim+'_'+run+'_'+'group'+str(n)+'_gas.png')
    plt.show()


