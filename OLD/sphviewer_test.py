import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import os
from sys import exit
from sys import argv
from tqdm import tqdm
import sphviewer
from mpl_toolkits.axes_grid1 import make_axes_locatable

sim = 'L0025N0376'
run = 'REFERENCE'
halotype = 'all'

prop = str(argv[1])

if prop == 'temperature':
    mean_mass_weighted = True
else:
    mean_mass_weighted = False

save_loc = '/home/arijdav1/Dropbox/phd/figures/sph_movies/'+prop+'/'

if not os.path.exists(save_loc):
    os.makedirs(save_loc)

res = 512
ext = 400.

px_size = ext/float(res) # in kpc
px_size_pc = px_size * 1e3

for n in np.arange(0,5000):

    if n != 32:
        continue    

    try:
        data = h5.File('/data5/arijdav1/saved_regions/'+sim+'_'+run+'/028_z000p000/'+halotype+'/group'+str(n)+'.hdf5','r')
    except IOError:
        continue

    if prop == 'stars':
        column_pos = np.array(data['Stars/Coordinates'])
        smoothing_length = np.array(data['Stars/SmoothingLength'])
        quantity = np.array(data['Stars/Mass'])*1e10
    else:
        column_pos = np.array(data['Gas/Coordinates'])
        smoothing_length = np.array(data['Gas/SmoothingLength'])
        if prop == 'gasmass':
            quantity = np.array(data['Gas/Mass'])*1e10
        elif prop == 'xrays':
            quantity = np.array(data['Gas/Xray_luminosity'])
        elif prop == 'temperature':
            mass = np.array(data['Gas/Mass'])*1e10
            quantity = np.float32(mass*np.array(data['Gas/Temperature'])) # recast just in case

    N = len(quantity)
    
    pos = np.zeros((3,N))
    for i in range(3): # Arrange co-ordinates into 3 rows rather than 3 colums for sph-viewer
        pos[i,:] = column_pos[:,i]

    pos *= 1e3 # convert to kpc
    smoothing_length *= 1e3
    
    if prop == 'gasmass':
        cmap = 'hot'
        axlab = r'$\log\Sigma_g$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$'
    elif prop == 'xrays':
        cmap = 'gnuplot'
        axlab = r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{pc}^{-2}]$'
    elif prop == 'stars':
        cmap = 'bone'
        axlab = r'$\log\Sigma_*$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$'
    elif prop == 'temperature':
        cmap = 'seismic'
        axlab = r'$\langle\log(T)\rangle$ $[\mathrm{K}]$'
    else:
        raise IOError('Plot options are "gasmass", "stars", "temperature" or "xrays"')

    if not mean_mass_weighted:

        Particles = sphviewer.Particles(pos,quantity,hsml=smoothing_length)
        Scene = sphviewer.Scene(Particles)
        
        phi_list = np.arange(0.,359.,1)
        for p, phi in tqdm(enumerate(phi_list)):
            Scene.update_camera(x=0.,y=0.,z=0.,r='infinity',t=45.,p=phi,extent=[-ext,ext,-ext,ext],xsize=res,ysize=res)
            Render = sphviewer.Render(Scene)
            Render.set_logscale()

            img = Render.get_image()
            img -= np.log10(px_size_pc**2)

            
            if p == 0:
                vmin = np.amin(img)
                vmax = np.amax(img)

                if prop == 'stars':
                    vmin = -0.5
                        extent = Render.get_extent()            fig = plt.figure(1,figsize=(8,8))            ax1 = fig.add_subplot(111)            im1 = plt.imshow(img, extent=extent, origin='lower', cmap=cmap,vmin=vmin,vmax=vmax)
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.15)
            col = plt.colorbar(im1, cax=cax)
            col.set_label(axlab,fontsize=16)
            plt.savefig(save_loc+'group'+str(n)+'_phi%03d'%(p)+'.png')
            plt.close()

    else:
        mass_Particles = sphviewer.Particles(pos,mass,hsml=smoothing_length)
        quantity_Particles = sphviewer.Particles(pos,quantity,hsml=smoothing_length)
        mass_Scene = sphviewer.Scene(mass_Particles)
        quantity_Scene = sphviewer.Scene(quantity_Particles)
    
        phi_list = np.arange(0.,359.,1)
        for p, phi in tqdm(enumerate(phi_list)):
            mass_Scene.update_camera(x=0.,y=0.,z=0.,r='infinity',t=45.,p=phi,extent=[-ext,ext,-ext,ext],xsize=res,ysize=res)
            quantity_Scene.update_camera(x=0.,y=0.,z=0.,r='infinity',t=45.,p=phi,extent=[-ext,ext,-ext,ext],xsize=res,ysize=res)
            mass_Render = sphviewer.Render(mass_Scene)
            #mass_Render.set_logscale()
            quantity_Render = sphviewer.Render(quantity_Scene)
            #quantity_Render.set_logscale()

            mass_img = mass_Render.get_image()
            quantity_img = quantity_Render.get_image()

            img = np.log10(quantity_img/mass_img)

            if p == 0:
                vmin = np.amin(img)
                vmax = np.amax(img)
                
                if prop == 'temperature':
                    vmin = 4.
                    vmax = 6.5
                            extent = quantity_Render.get_extent()            fig = plt.figure(1,figsize=(8,8))            ax1 = fig.add_subplot(111)            im1 = plt.imshow(img, extent=extent, origin='lower', cmap=cmap,vmin=vmin,vmax=vmax)
            divider = make_axes_locatable(ax1)
            cax = divider.append_axes("right", size="5%", pad=0.15)
            col = plt.colorbar(im1, cax=cax)
            col.set_label(axlab,fontsize=16)
            plt.savefig(save_loc+'group'+str(n)+'_phi%03d'%(p)+'.png')
            plt.close()



