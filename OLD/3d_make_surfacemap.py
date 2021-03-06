import numpy as np
import h5py as h5
from sys import exit
from sys import argv
from python_tools import *
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from scipy import stats
import pickle
import os

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

def full_kernel(particle_xyz,gridpoint_xyz,h): # For a given particle, returns the SPH kernel on a 2D grid
    grid_size = int(np.power(len(gridpoint_xyz[:,0]),1./3.))+1
    no_pixels = int(len(gridpoint_xyz[:,0]))
    kernel = np.zeros(no_pixels)
    ones = np.reshape(np.ones(no_pixels),(no_pixels,1))
    particle_xyzs = ones*particle_xyz
    r = np.sqrt((abs(gridpoint_xyz[:,0]-particle_xyzs[:,0]))**2 + (abs(gridpoint_xyz[:,1]-particle_xyzs[:,1]))**2 + (abs(gridpoint_xyz[:,2]-particle_xyzs[:,2]))**2)
    mask = np.where(r<=h)
    kernel[mask] += (21/(np.pi*h**3)) * (1-r[mask]/h)**4 * (1+4*r[mask]/h)
    kernel_out = np.reshape(kernel,(grid_size,grid_size,grid_size))
    return kernel_out

def dimensionless_sph_weight(mass,density,h,n_dims):
    return mass/(density*(h**n_dims))

##############################################################################################################
# SET RESOLUTION
res = 256
sim = str(argv[1])
run = str(argv[2])
#frame_size = 70 # In kpc

print sim
print run

do_stars = True
show_gals = False

halotype = str(argv[3])

if 'nosfr' in argv:
    nosfr = True    # Pick only SFR=0 gas for X ray maps?
else:
    nosfr = False

if 'onlysfr' in argv:
    onlysfr = True    # Pick only SFR=0 gas for X ray maps?
else:
    onlysfr = False

if 'edgeon' in argv:
    edgeon = True
    edgeflag = '_edgeon'
    edgefold = '/edgeon'
else:
    edgeon = False
    edgeflag = ''
    edgefold = '/faceon'

if 'zoom' in argv:
    zoom = True
    zoomflag = 'zoom/'
else:
    zoom = False
    zoomflag = 'full/'

##############################################################################################################
ct = 0
for n in np.arange(0,500):
    ct += 1
    
    if n != 32:
        continue
    
    try:
	    data = h5.File('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/028_z000p000/'+halotype+'/group'+str(n)+'.hdf5','r+')
    except IOError:
	    continue

    print 'Group no. ',n

    # Establish frame size
    r200 = np.float64(data['Volume/r200'])*1e3 # Convert r200 to kpc

    if zoom:
        dz = r200/4.
        frame_size = r200/4.
    else:
        dz = r200*2.
        frame_size = r200*2. # go out to log(r/r200) = 0.5
    pixel_size = 2.*np.float(frame_size)/np.float(res)

    save_loc = '/data5/arijdav1/saved_maps/'+sim+'_'+run+'/'+halotype+edgefold+'/'+zoomflag

    if not os.path.exists(save_loc):
        os.makedirs(save_loc)

    of = save_loc+'group'+str(n)+'.pkl'
    outfile = open(of, 'w')

    star_coords = np.array(data['Stars/Coordinates'])
    star_h = np.array(data['Stars/SmoothingLength'])
    star_mass = np.array(data['Stars/Mass'])
    star_density = np.array(data['Stars/Density'])

    gas_coords = np.array(data['Gas/Coordinates'])
    gas_h = np.array(data['Gas/SmoothingLength'])
    gas_mass = np.array(data['Gas/Mass'])
    gas_density = np.array(data['Gas/Density'])
    gas_metals = np.array(data['Gas/Metallicity'])
    gas_temp = np.array(data['Gas/Temperature'])
    sfr = np.array(data['Gas/StarFormationRate'])
    gas_Lx = np.array(data['Gas/Xray_luminosity'])
    
    # Conversions from whatever is in the data
    star_coords *= 1e3
    gas_coords *= 1e3
    star_h *= 1e3
    gas_h *= 1e3

    star_mass *= 1e10
    gas_mass *= 1e10
    star_density *= 10
    gas_density *= 10
    

    n_p = len(gas_coords[:,0])
    n_s = len(star_coords[:,0])

    if edgeon:
        gas_coords = np.hstack((np.reshape(gas_coords[:,0],(n_p,1)),np.reshape(gas_coords[:,2],(n_p,1)),np.reshape(gas_coords[:,1],(n_p,1))))
        star_coords = np.hstack((np.reshape(star_coords[:,0],(n_s,1)),np.reshape(star_coords[:,2],(n_s,1)),np.reshape(star_coords[:,1],(n_s,1))))

    # Set minimum smoothing length to half a pixel
    star_h[star_h<pixel_size/2.] = pixel_size/2.    
    gas_h[gas_h<pixel_size/2.] = pixel_size/2.
    
    # Remove all particles that aren't within h of the frame or within dz range
    mask1 = (abs(gas_coords[:,0])<=frame_size)&(abs(gas_coords[:,1])<=frame_size)&(abs(gas_coords[:,2])<=dz)
    gas_coords = gas_coords[mask1]
    gas_h = gas_h[mask1]
    gas_mass = gas_mass[mask1]
    gas_density = gas_density[mask1]
    gas_temp = gas_temp[mask1]
    sfr = sfr[mask1]
    gas_metals = gas_metals[mask1]
    gas_Lx = gas_Lx[mask1]

    mask2 = (abs(star_coords[:,0])<=frame_size)&(abs(star_coords[:,1])<=frame_size)&(abs(star_coords[:,2])<=dz)
    star_coords = star_coords[mask2]
    star_h = star_h[mask2]
    star_mass = star_mass[mask2]
    star_density = star_density[mask2]

    if nosfr: # Take only NSF gas?
        mask3 = np.where(sfr==0)[0]
        gas_coords = gas_coords[mask3]
        gas_h = gas_h[mask3]
        gas_mass = gas_mass[mask3]
        gas_density = gas_density[mask3]
        gas_temp = gas_temp[mask3]
        sfr = sfr[mask3]
        gas_metals = gas_metals[mask3]
        gas_Lx = gas_Lx[mask3]

    if onlysfr: # Take only SF gas?
        mask3 = np.where(sfr>0.)[0]
        gas_coords = gas_coords[mask3]
        gas_h = gas_h[mask3]
        gas_mass = gas_mass[mask3]
        gas_density = gas_density[mask3]
        gas_temp = gas_temp[mask3]
        sfr = sfr[mask3]
        gas_metals = gas_metals[mask3]
        gas_Lx = gas_Lx[mask3]

    no_star = len(star_coords[:,0])
    no_gas = len(gas_coords[:,0])
    xedges = np.linspace(-frame_size,frame_size,num=res+1)
    yedges = np.linspace(-frame_size,frame_size,num=res+1)
    zedges = np.linspace(-dz,dz,num=res+1)
    bs = get_binsizes(xedges)[0]
    zbs = get_binsizes(zedges)[0]

    xcents = get_bincentres(xedges)
    ycents = get_bincentres(yedges)
    zcents = get_bincentres(zedges)
    gridx, gridy, gridz = np.meshgrid(xcents, ycents, zcents)
    xlist = np.reshape(np.ravel(gridx),(len(np.ravel(gridx)),1))
    ylist = np.reshape(np.ravel(gridy),(len(np.ravel(gridy)),1))
    zlist = np.reshape(np.ravel(gridz),(len(np.ravel(gridz)),1))
    gxys = np.hstack((xlist,ylist,zlist))

    gasmap = np.zeros((res,res,res))
    sfrmap = np.zeros((res,res,res))
    metalmap = np.zeros((res,res,res))
    tempmap = np.zeros((res,res,res))
    starmap = np.zeros((res,res,res))
    xraymap = np.zeros((res,res,res))
    
    print 'Creating gas-based maps...'

    for p in tqdm(range(no_gas)):
        kern = full_kernel(gas_coords[p,:],gxys,gas_h[p]) #######
        if np.sum(kern)<1e-10:
            continue
        gasmap += gas_mass[p]*kern # surface density
        tempmap += gas_mass[p]*gas_temp[p]*kern # mass-weighted mean
        sfrmap += gas_mass[p]*sfr[p]*kern # mass-weighted mean
        metalmap += gas_mass[p]*gas_metals[p]*gas_mass[p]*kern # mass-weighted mean
        xraymap += gas_Lx[p]*kern # surface density
    
    if do_stars == True:
        print 'Creating star-based maps...'
        for p in tqdm(range(no_star)):
            star_kern = full_kernel(star_coords[p,:],gxys,star_h[p])
            if np.sum(star_kern)<1e-10:
                continue
            starmap += star_mass[p]*star_kern # surface density
    else:
        starmap += np.ones((res,res))

    gasmap = np.sum(gasmap,axis=0)
    tempmap = np.sum(tempmap,axis=0)
    sfrmap = np.sum(sfrmap,axis=0)
    metalmap = np.sum(metalmap,axis=0)
    xraymap = np.sum(xraymap,axis=0)
    starmap = np.sum(starmap,axis=0)

    gmassmap = gasmap*pixel_size**2
    starmassmap = starmap*pixel_size**2
    rawLxmap = xraymap*pixel_size**2
    
    # Finish off mass-weighting certain quantities
    tempmap *= pixel_size**2 / gmassmap
    sfrmap *= pixel_size**2 / gmassmap
    metalmap *= pixel_size**2 / gmassmap

    gasmap /= 1e6 # Convert to Msun pc-2
    starmap /= 1e6 # Convert to Msun pc-2
    metalmap /= 1e6 # Convert to Msun pc-2
    xraymap /= 1e6 # Convert to erg s-1 pc-2

    image_total = np.sum(gmassmap)
    imagestar_total = np.sum(starmassmap)
    imagexray_total = np.sum(rawLxmap)  

    # Normalisation test
    print 'Total mass of gas particles: ',np.sum(gas_mass)
    print 'Total mass smoothed onto image: ',image_total
    print 'Total mass of star particles: ',np.sum(star_mass)
    print 'Total stellar mass smoothed onto image: ',imagestar_total
    print 'Total L_x of particles: ',np.sum(gas_Lx)
    print 'Total L_x smoothed onto image: ',imagexray_total
    print '\n'
    print 'Saving...'

    nos,xe,ye = np.histogram2d(gas_coords[:,0],gas_coords[:,1],[xedges,yedges])
    
    map_out = dict([('sim',sim),('run',run),('frame_size',frame_size),('res',res),('dz',dz),('pixel_size',pixel_size),('group',n),('gxys',gxys),('gasmap',gasmap),('starmap',starmap),('metalmap',metalmap),('sfrmap',sfrmap),('tempmap',tempmap),('xraymap',xraymap),('hist',nos)])

    pickle.dump(map_out,outfile)
    outfile.close()
    print 'Done! Data saved to ',of

    
    
