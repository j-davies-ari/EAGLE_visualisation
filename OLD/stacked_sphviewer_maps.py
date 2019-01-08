import numpy as np
import h5py as h5
from sys import exit
from python_tools import *
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy import stats
import os
import eagle as E
from find_closest_galaxy import groupnum_table
from sys import argv
import sphviewer
from mpl_toolkits.axes_grid1 import make_axes_locatable

def load(field,dictionary):
    return np.array(dictionary[field])

sim = 'L0025N0376'
run = 'REFERENCE'
snap = '028_z000p000'
infotype = 'SUBFIND'
table = 'Subhalo'
halotype = 'all'
galtype = 'Lstar'

prop = argv[1]

save_loc = '/home/arijdav1/Dropbox/phd/figures/stacked_images/'+prop

res = 512
ext = 250.

'''

# Pick a mass range and select the group numbers to load in
simstr = '/data5/simulations/EAGLE/'+sim+'/'+run+'/data/'
num_subs = np.array(E.readArray("SUBFIND_GROUP", simstr, snap, "/FOF/NumOfSubhalos"))
masslist = np.array(E.readArray("SUBFIND_GROUP",simstr,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10

if isolate_galaxies:
    numbers = np.arange(len(masslist)) + 1
    groupmass = np.array(E.readArray("SUBFIND_GROUP",simstr,snap,'FOF/GroupMass')[num_subs>0])
    subhalo_masses = np.array(E.readArray("SUBFIND", simstr, snap, "/Subhalo/Mass"))
    subhalo_groupnums = np.array(E.readArray("SUBFIND", simstr, snap, "/Subhalo/GroupNumber"))

    halo_indices = []
    print 'Picking isolated galaxies...'
    for g in tqdm(range(len(numbers))):
        sub_m = subhalo_masses[subhalo_groupnums==numbers[g]]
        if len(sub_m)==0:
            continue
        if sub_m[0]/groupmass[g] > 0.9:
            halo_indices.append(g)
    halo_indices = np.array(halo_indices)

    masslist = masslist[halo_indices]
    numbers = numbers[halo_indices]

    del groupmass
    del subhalo_masses
    del subhalo_groupnums

if galtype=='dwarf':
    mask = np.where((masslist>=10**11.3)&(masslist<=10**11.5))[0]
elif galtype=='Lstar':
    mask = np.where((masslist>=10**11.9)&(masslist<=10**12.25))[0]
elif galtype=='group':
    mask = np.where((masslist>=10**12.6)&(masslist<=10**13.2))[0]
elif galtype=='all':
    mask = np.where((masslist>=10**10)&(masslist<=10**15))[0]
if isolate_galaxies:
    gns = numbers[mask]
    del numbers
else:
    gns = mask+1
del num_subs
del masslist
'''

# Load in pre-processed halo data
halofile = h5.File('/data5/arijdav1/halo_data/'+sim+'_'+run+'_'+snap+'.hdf5','r')
gasdata = halofile['All_CGM']

# Load halo properties
Mstar_30kpc = load('Mstar_30kpc',halofile)*1e10
M200 = load('M200',halofile)*1e10
gns = load('GroupNumber',halofile)

isolated_flag = np.array(halofile['IsolatedFlag'])
mask = np.where((isolated_flag>0.)&(M200>np.power(10.,11.5)))[0]









px_size = ext/float(res) # in kpc
px_size_pc = px_size * 1e3

stacked_img = np.zeros((res,res),dtype=np.float64)

print 'Making and stacking maps...'
# Make each surface map and stack it onto the rest at each step
for n in tqdm(range(len(gns))):
    try:
        data = h5.File('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/028_z000p000/'+halotype+'/group'+str(gns[n])+'.hdf5','r')
    except IOError:
        continue

    if prop == 'stars':
        column_pos = np.array(data['Stars/Coordinates'])
        smoothing_length = np.array(data['Stars/SmoothingLength'])
        quantity = np.array(data['Stars/Mass'])*1e10
    else:

        # Define hot gas as T > 2x10^4 K
        temp = np.array(data['Gas/Temperature'])
        L_x = np.array(data['Gas/Xray_luminosity'])
        hot = np.where(L_x>0.)[0]

        column_pos = np.array(data['Gas/Coordinates'])[hot,:]
        smoothing_length = np.array(data['Gas/SmoothingLength'])[hot]
        if prop == 'gasmass':
            quantity = np.array(data['Gas/Mass'])[hot]*1e10
        elif prop == 'xrays':
            quantity = (L_x[hot].astype(np.float64))/1e30


    N = len(quantity)
    
    pos = np.zeros((3,N))
    for i in range(3): # Arrange co-ordinates into 3 rows rather than 3 colums for sph-viewer
        pos[i,:] = column_pos[:,i]

    pos *= 1e3 # convert to kpc
    smoothing_length *= 1e3

    if prop == 'gasmass':
        cmap = 'hot'
        axlab = r'$\log\Sigma_g$ [$M_{\odot}$ $\mathrm{kpc}^{-2}]$'
    elif prop == 'xrays':
        cmap = 'gnuplot'
        axlab = r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{kpc}^{-2}]$'
    elif prop == 'stars':
        cmap = 'bone'
        axlab = r'$\log\Sigma_*$ [$M_{\odot}$ $\mathrm{kpc}^{-2}]$'
    else:
        raise IOError('Plot options are "gasmass", "stars", "temperature" or "xrays"')

    Particles = sphviewer.Particles(pos,quantity,hsml=smoothing_length)
    Scene = sphviewer.Scene(Particles)
    
    Scene.update_camera(x=0.,y=0.,z=0.,r='infinity',t=90.,extent=[-ext,ext,-ext,ext],xsize=res,ysize=res)
    Render = sphviewer.Render(Scene)
    #Render.set_logscale()

    img = Render.get_image()
    sdensity_img = img/px_size**2

    '''
    # For quickly viewing the individual images
    extent = Render.get_extent()
    plt.figure()
    plt.imshow(sdensity_img,cmap='hot',extent=extent, origin='lower')
    plt.show()
    '''

    '''
    fig = plt.figure(1,figsize=(5,5))    ax1 = fig.add_subplot(111)
    Render.set_logscale()    Render.histogram(bins=100, log=True)
    plt.show()
    '''
    stacked_img += sdensity_img    

if prop == 'xrays':
    stacked_img *= 1e30
    #stacked_img[stacked_img>1e33] = 1e33 # to set a max


stacked_img = np.log10(stacked_img/len(gns)) # Divide through by number of halos stacked and set a logscale
extent = Render.get_extent()fig = plt.figure(1,figsize=(8,8))ax1 = fig.add_subplot(111)
if prop == 'stars':    im1 = plt.imshow(stacked_img, extent=extent, origin='lower', cmap=cmap, vmin=0.)
else:
    im1 = plt.imshow(stacked_img, extent=extent, origin='lower', cmap=cmap)
divider = make_axes_locatable(ax1)
cax = divider.append_axes("right", size="5%", pad=0.15)
col = plt.colorbar(im1, cax=cax)
col.set_label(axlab,fontsize=16)
plt.savefig(save_loc+'_map.png')
plt.close()

plt.figure(figsize=(8,6))
plt.hist(np.ravel(stacked_img),bins=50,histtype='step',log=True)
plt.ylim(0.9,1e6)
plt.savefig(save_loc+'_hist.png')
plt.show()



def radial_sdensity_profile(image,extent,resolution,nbins=50,plot_to=250.,axis_type='linear'):
    # Find the radial surface density profile

    # Make a grid of pixel radii from the centre (in kpc)
    gridx, gridy = np.meshgrid(np.linspace(-extent,extent,num=resolution),np.linspace(-extent,extent,num=resolution))
    radial_grid = np.sqrt(gridx**2+gridy**2)

    # Flatten the grid and corresponding image values and sort them by radius
    flat_radial_grid = np.ravel(radial_grid)
    flat_img = np.ravel(np.float64(10**image))
    sort_indices = np.argsort(flat_radial_grid)
    flat_radial_grid = flat_radial_grid[sort_indices]
    flat_img = flat_img[sort_indices]

    # Define the radial bins
    if axis_type == 'linear':
        radial_binedges = np.linspace(0.,plot_to,nbins+1)
        radial_bincentres = get_bincentres(radial_binedges)
    elif axis_type == 'loglog':
        radial_binedges = np.logspace(0.,np.log10(plot_to),nbins+1)
        radial_bincentres = get_bincentres(np.log10(radial_binedges))
    else:
        raise ValueError('Please pick a valid axis type ("linear" or "loglog")')

    median_sdensities = np.zeros(nbins)
    for b in range(nbins):
        startstop = np.searchsorted(flat_radial_grid,[radial_binedges[b],radial_binedges[b+1]])
        median_sdensities[b] = np.median(flat_img[startstop[0]:startstop[1]])

        if axis_type == 'loglog' and b<10:
            print flat_img[startstop[0]:startstop[1]]

    norm = np.amax(median_sdensities[~np.isnan(median_sdensities)])
    '''
    if axis_type == 'loglog':
        print np.log10(median_sdensities)
        print np.log10(norm)
    '''
    if axis_type == 'linear':
        return median_sdensities/norm, radial_bincentres
    else:   
        return np.array(np.log10(median_sdensities)-np.log10(norm)), np.array(radial_bincentres)

profile,bcs = radial_sdensity_profile(stacked_img,ext,res,axis_type='linear')

logprofile,logbcs = radial_sdensity_profile(stacked_img,ext,res,axis_type='loglog')


fig,ax1 = plt.subplots()
left,bottom,width,height = [0.6,0.6,0.25,0.25]
ax2 = fig.add_axes([left, bottom, width, height])

ax1.plot(bcs,profile,c='r')
ax2.plot(10**logbcs,10**logprofile,c='r')
ax2.set_xscale('log')
ax2.set_yscale('log')

ax1.set_xlabel('$r$ $[\mathrm{kpc}]$')
ax2.set_xlabel('$r$ $[\mathrm{kpc}]$')
ax1.set_ylabel('$\mathrm{Normalised}$ $\Sigma_{\mathrm{median}}$')
ax2.set_ylabel('$\mathrm{Norm.}$ $\Sigma_{\mathrm{med.}}$')

plt.savefig(save_loc+'_betaprofile.png')
plt.show()





