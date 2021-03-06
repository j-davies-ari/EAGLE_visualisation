import numpy as np
import h5py as h5
from sys import exit
from sys import argv
from python_tools import *
import matplotlib.pyplot as plt
import random
from tqdm import tqdm
from scipy import stats
import eagle as E
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab

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

##############################################################################################################
# SET RESOLUTION
res = 256
sim = str(argv[1])
run = str(argv[2])

simloc = '/data5/simulations/EAGLE/'+sim+'/'+run+'/data/'
snap = '028_z000p000'

ap_enclosed = 300
show_gals = False

halotype = str(argv[3])

if 'edgeon' in argv:
    edgeon=True
    edgeflag='_edgeon'
    edgefold='/edgeon'
else:
    edgeon=False
    edgeflag=''
    edgefold='/faceon'

if 'zoom' in argv:
    zoom = True
    zoomflag = 'zoom/'
else:
    zoom = False
    zoomflag = 'full/'

find_median_sdensities = False

#loc = '/data5/arijdav1/saved_maps/'+halotype+'/'
loc = '/data5/arijdav1/saved_maps/'+sim+'_'+run+'/'+halotype+edgefold+'/'+zoomflag
save_loc = '/home/arijdav1/Dropbox/phd/figures/sph_pics/'+sim+'_'+run+'/'+halotype+edgefold+'/'+zoomflag

if not os.path.exists(save_loc):
    os.makedirs(save_loc)
##############################################################################################################

num_subs = np.array(E.readArray("SUBFIND_GROUP",simloc,snap, "/FOF/NumOfSubhalos"))
masslist = np.array(E.readArray("SUBFIND_GROUP",simloc,snap,'FOF/Group_M_Crit200')[num_subs>0])*1e10
gns = np.where(num_subs>0)[0]+1

subfind_subgn = E.readArray('SUBFIND',simloc,snap,'/Subhalo/SubGroupNumber')
subfind_gn = E.readArray('SUBFIND',simloc,snap,'/Subhalo/GroupNumber')[subfind_subgn==0]
sfrs = E.readArray('SUBFIND',simloc,snap,'/Subhalo/ApertureMeasurements/SFR/030kpc')[subfind_subgn==0]

ct = 0
sdensities = []
sfr_sdensities = []
metals_sdensities = []
stars_sdensities = []
xray_sdensities = []
for n in np.arange(0,500):
    ct += 1
    
    if n != 32:
        continue
    
    try:
        alldata = h5.File('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/028_z000p000/'+halotype+'/group'+str(n)+'.hdf5','r+')
    except IOError:
	    continue

    print 'Group no. ',n

    # Get M_200 and SFR_30kpc
    m200 = np.log10(masslist[gns==n])
    sfr_30kpc = sfrs[subfind_gn==n]

    # Establish frame size
    r200 = np.float64(alldata['Volume/r200'])*1e3 # Convert r200 to kpc
    frame_size = r200*2. # go out to 2r200
    pixel_size = 2.*np.float(frame_size)/np.float(res)
    dz = r200*2.

    if zoom:
        frame_size = r200/4. # go out to 1/2 r200
        pixel_size = 2.*np.float(frame_size)/np.float(res)
        dz = r200/4.

    del alldata

    try:
        data = mapload('group'+str(n),loc)
    except IOError:
	    continue

    

    #pixel_size = data['pixel_size']
    gxys = data['gxys']
    gasmap = data['gasmap']
    starmap = data['starmap']
    metalmap = data['metalmap']
    tempmap = data['tempmap']
    sfrmap = data['sfrmap']
    xraymap = data['xraymap']
    nos = data['hist']

    pixel_size_pc = pixel_size*1000

    gmassmap = gasmap*pixel_size_pc**2  ### Note that pixel sizes in pc must be used for everything but SFR
    rawsfrmap = sfrmap*pixel_size**2
    starmassmap = starmap*pixel_size_pc**2
    rawmetalmap = metalmap*pixel_size_pc**2
    rawxraymap = xraymap*pixel_size_pc**2

    circle,npix = circle_kernel([0.,0.,0.],gxys,r200)
    Lx_r200 = np.log10(np.sum(rawxraymap*circle))

    gasmap[gasmap<1e-4]=1e-3
    #sfrmap[sfrmap<1e-7]=1e-7
    #metalmap[metalmap<1e-4]=1e-4
    starmap[starmap<1e-1]=1e-1
    xraymap[xraymap<1e20]=1e20
    
    log_gasmap = np.log10(gasmap)
    log_sfrmap = np.log10(sfrmap)
    log_metalmap = np.log10(metalmap)
    log_tempmap = np.log10(tempmap)
    log_starmap = np.log10(starmap)
    log_xraymap = np.log10(xraymap)

    plt.figure()
    plt.hist(np.ravel(log_xraymap))
    plt.show()
    exit()


    if find_median_sdensities:
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
                    masked_sfr = rawsfrmap*circle
                    masked_metals = rawmetalmap*circle
                    masked_stars = starmassmap*circle
                    masked_xray = rawxraymap*circle
                    sdensities.append(np.sum(masked_mass)/(npix*pixel_size_pc**2))
                    if np.sum(masked_sfr)<1e-8:
                        sfr_sdensities.append(1e-8)
                    else:
                        sfr_sdensities.append(np.sum(masked_sfr)/(npix*pixel_size**2))
                    stars_sdensities.append(np.sum(masked_stars)/(npix*pixel_size_pc**2))
                    metals_sdensities.append(np.sum(masked_metals)/(npix*pixel_size_pc**2))
                    xray_sdensities.append(np.sum(masked_xray)/(npix*pixel_size_pc**2))
                    #print 'Final radius: ',ri,' kpc'
                    break
            #print count,' iterations'
        
        print 'Median surface density in gas: ',np.median(sdensities)
        print 'Median surface density in SFR: ',np.median(sfr_sdensities)
        print 'Median surface density in metals: ',np.median(metals_sdensities)
        print 'Median surface density in stars: ',np.median(stars_sdensities)
        print 'Median surface density in X-rays: ',np.median(xray_sdensities)

    # Show the galaxy

    textstr = '$\logM_{200}=%.2f$ $[\mathrm{M}_{\odot}]$\n$\logM_*=%.2f$ $[\mathrm{M}_{\odot}]$\n$\dot{M}_{*,30\mathrm{kpc}}=%.2f$ $[\mathrm{M}_{\odot}$ $\mathrm{yr}^{-1}]$\n$\logL_X(r<r_{200})=%.2f$ $[\mathrm{erg}$ $\mathrm{s}^{-1}]$'%(m200,np.log10(np.sum(starmassmap)),sfr_30kpc,Lx_r200)
    props = dict(boxstyle='round', facecolor='w', alpha=0.5)
    
    r200_circle = pylab.Circle((0,0), radius=r200, color='w',ls='--', fill=False)

    fig0 = plt.figure(figsize=(10,10))
    ax0 = plt.gca()
    im1 = plt.imshow(log_gasmap,cmap='gist_heat',vmin=-2.5,vmax=2.6,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax0)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\log\Sigma_g$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax0.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_gas.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig0)
    
    fig1 = plt.figure(figsize=(10,10))
    ax1 = plt.gca()
    im1 = plt.imshow(sfrmap,cmap='gist_rainbow',vmin=0.,vmax=0.01,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\langle \mathrm{SFR}\rangle$ $[M_{\odot}$ $\mathrm{yr}^{-1}]$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax1.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_sfr.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig1)

    fig2 = plt.figure(figsize=(10,10))
    ax2 = plt.gca()
    im1 = plt.imshow(log_tempmap,cmap='hot',vmin=3.9,vmax=7.,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\log\langle T\rangle$ $[\mathrm{K}])$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax2.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_temp.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig2)
    
    fig3 = plt.figure(figsize=(10,10))
    ax3 = plt.gca()
    im1 = plt.imshow(log_metalmap,cmap='gist_rainbow',vmin=-4.5,vmax=1.3,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\log\langle Z\rangle$ $[Z_{\odot}]$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax3.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_metals.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig3)
    
    fig4 = plt.figure(figsize=(10,10))
    ax4 = plt.gca()
    im1 = plt.imshow(log_starmap,cmap='bone',vmin=-1.,vmax=3.,extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\log\Sigma_*$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax4.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_star.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig4)

    fig5 = plt.figure(figsize=(10,10))
    ax5 = plt.gca()
    im1 = plt.imshow(log_xraymap,cmap='gnuplot',extent=[-frame_size,frame_size,-frame_size,frame_size])
    plt.xlabel(r'$x$ $\mathrm{kpc}$',fontsize=16)
    plt.ylabel(r'$y$ $\mathrm{kpc}$',fontsize=16)
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes("right", size="5%", pad=0.15)
    col = plt.colorbar(im1, cax=cax)
    col.set_label(r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{pc}^{-2}]$',fontsize=16)
    plt.text(0.05, 0.95,textstr,transform=ax5.transAxes,fontsize=14,verticalalignment='top',bbox=props)
    plt.savefig(save_loc+'group'+str(n)+'_Xray.png',bbox_inches='tight')
    if show_gals:
        plt.show()
    plt.close(fig5)

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

exit()

# KENNICUTT SCHMITT STUFF


#sdensities = np.array(sdensities)/1e6 # Convert to Msun pc-2

KS_x = np.logspace(-2,3,num=100)
KS_y = 1.4*np.log10(KS_x) + np.log10(1.515e-4)

log_sdens = np.log10(sdensities)
log_sfrdens = np.log10(sfr_sdensities)

#log_sfrdens[log_sfrdens<-6] = -6

slope, intercept, r_value, p_value, std_err = stats.linregress(log_sdens[log_sfrdens>-8],log_sfrdens[log_sfrdens>-8])

#slope, intercept, r_value, p_value, std_err = stats.linregress(log_sdens,log_sfrdens)

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
plt.ylim(-8.5,1)
plt.savefig(save_loc+str(res)+'_KS_TEST_apN'+str(ap_enclosed)+'_dz'+str(dz)+'.png')
if show_gals:
    plt.show()
