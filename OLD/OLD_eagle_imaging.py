import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
import os
from sys import exit
from sys import argv
from tqdm import tqdm
import sphviewer
import eagle as E
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpi4py import MPI
import read_eagle as read
from copy import deepcopy
import safe_colours
import numexpr as ne
import matplotlib.offsetbox
from matplotlib.lines import Line2D
safe_colours = safe_colours.initialise()

class AnchoredHScaleBar(matplotlib.offsetbox.AnchoredOffsetbox):
    """ size: length of bar in data units
        extent : height of bar ends in axes units """
    def __init__(self, size=1, extent = 0.03, label="", loc=2, ax=None,
                 pad=0.4, borderpad=0.5, ppad = 0, sep=2, prop=None, 
                 frameon=True, **kwargs):
        if not ax:
            ax = plt.gca()
        trans = ax.get_xaxis_transform()
        size_bar = matplotlib.offsetbox.AuxTransformBox(trans)
        line = Line2D([0,size],[0,0], **kwargs)
        vline1 = Line2D([0,0],[-extent/2.,extent/2.], **kwargs)
        vline2 = Line2D([size,size],[-extent/2.,extent/2.], **kwargs)
        size_bar.add_artist(line)
        size_bar.add_artist(vline1)
        size_bar.add_artist(vline2)
        txt = matplotlib.offsetbox.TextArea(label, minimumdescent=False)
        self.vpac = matplotlib.offsetbox.VPacker(children=[size_bar,txt],align="center", pad=ppad, sep=sep) 
        matplotlib.offsetbox.AnchoredOffsetbox.__init__(self, loc, pad=pad, borderpad=borderpad, child=self.vpac, prop=prop, frameon=frameon)

def value_locate(refx, x): # Give this an array of x's and it returns the closest indices where they appear in refx
	refx = np.array(refx)
	x = np.array(x)
	loc = np.zeros(len(x), dtype='int')
	for i in xrange(len(x)):
		ix = x[i]
		ind = ((refx - ix) <= 0).nonzero()[0]
		if len(ind) == 0:
			loc[i] = 0.
		else: loc[i] = ind[-1]
	return loc

##################################################################################################################################

snapdict = dict([('28', ['028_z000p000', 0.]), \
                 ('27', ['027_z000p101', 0.101]), \
                 ('26', ['026_z000p183', 0.183]), \
                 ('25', ['025_z000p271', 0.271]), \
                 ('24', ['024_z000p366', 0.366]), \
                 ('23', ['023_z000p503', 0.503]), \
                 ('22', ['022_z000p615', 0.615]), \
                 ('21', ['021_z000p736', 0.736]), \
                 ('20', ['020_z000p865', 0.865]), \
                 ('19', ['019_z001p004', 1.004]), \
                 ('18', ['018_z001p259', 1.259]), \
                 ('17', ['017_z001p487', 1.487]), \
                 ('16', ['016_z001p737', 1.737]), \
                 ('15', ['015_z002p012', 2.012]), \
                 ('14', ['014_z002p237', 2.237]), \
                 ('13', ['013_z002p478', 2.478]), \
                 ('12', ['012_z003p017', 3.017]), \
                 ('11', ['011_z003p528', 3.528]), \
                 ('10', ['010_z003p984', 3.984]), \
                 ('9', ['009_z004p485', 4.485]), \
                 ('8', ['008_z005p037', 5.037]), \
                 ('7', ['007_z005p487', 5.487]), \
                 ('6', ['006_z005p971', 5.971]), \
                 ('5', ['005_z007p050', 7.05]), \
                 ('4', ['004_z008p075', 8.075]), \
                 ('3', ['003_z008p988', 8.988]), \
                 ('2', ['002_z009p993', 9.993]), \
                 ('1', ['001_z015p132', 15.132]), \
                 ('0', ['000_z020p000', 20.])])


##################################################################################################################################

class singlehalo(object):
    def __init__(self, gn, prop,
                 sim='L0100N1504',
                 run='REFERENCE',
                 snapnum=28,
                 halotype='all'):

        self.groupnumber = gn
        self.property = prop
        tag = snapdict[str(snapnum)][0]

        try:
            data = h5.File(
                '/data5/arijdav1/saved_regions/' + sim + '_' + run + '/' + tag + '/' + halotype + '/group' + str(
                    gn) + '.hdf5', 'r')
        except IOError:
            print 'Invalid group number for this box.'
            exit()

        if prop == 'stars':
            pos = np.array(data['Stars/Coordinates']).T
            smoothing_length = np.array(data['Stars/SmoothingLength'])
            quantity = np.array(data['Stars/Mass']) * 1e10
        else:
            pos = np.array(data['Gas/Coordinates']).T
            smoothing_length = np.array(data['Gas/SmoothingLength'])
            if prop == 'gas':
                quantity = np.array(data['Gas/Mass']) * 1e10
            elif prop == 'hotgas':
                quantity = np.array(data['Gas/Mass']) * 1e10
                temp = np.array(data['Gas/Temperature'])
                mask = np.where(temp > np.power(10., 5.5))[0]
                column_pos = column_pos[mask]
                smoothing_length = smoothing_length[mask]
                quantity = quantity[mask]
            elif prop == 'xrays':
                quantity = np.array(
                    data['Gas/Xray_luminosity']) / 1e30  # make the numbers smaller so sphviewer can deal with them
            else:
                raise IOError('Plot options are "gas", "hotgas", "stars" or "xrays"')


        N = len(quantity)

        pos *= 1e3  # convert to kpc
        smoothing_length *= 1e3

        if prop == 'gas':
            self.cmap = 'hot'
            self.axlab = r'$\log\Sigma_g$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$'
        elif prop == 'hotgas':
            self.cmap = 'hot'
            self.axlab = r'$\log\Sigma_{\mathrm{g,hot}}$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$'
        elif prop == 'xrays':
            self.cmap = 'gnuplot'
            self.axlab = r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{pc}^{-2}]$'
        elif prop == 'stars':
            self.cmap = 'bone'
            self.axlab = r'$\log\Sigma_*$ [$M_{\odot}$ $\mathrm{pc}^{-2}]$'
        else:
            raise IOError('Plot options are "gas", "hotgas", "stars" or "xrays"')

        Particles = sphviewer.Particles(pos, quantity, hsml=smoothing_length)
        self.Scene = sphviewer.Scene(Particles)

    def image(self, x=0.,
              y=0.,
              z=0.,
              theta=45.,
              phi=0.,
              extent=400.,
              resolution=1024,
              camera_distance='infinity',
              save=True,
              show=True,
              path='/home/arijdav1/Dropbox/phd/figures/sph_pics/'):

        px_size = extent / float(resolution)  # in kpc
        px_size_pc = px_size * 1e3

        self.Scene.update_camera(x=x, y=y, z=z, r=camera_distance, t=theta, p=phi,
                                 extent=[-extent, extent, -extent, extent], xsize=resolution, ysize=resolution)

        Render = sphviewer.Render(self.Scene)
        Render.set_logscale()
        extent = Render.get_extent()
        img = Render.get_image()

        if self.property == 'xrays':
            img += 30.

        img -= np.log10(px_size_pc ** 2)

        vmin = np.amin(img[np.isfinite(img)])
        vmax = np.amax(img[np.isfinite(img)])

        if self.property == 'stars':
            vmin = 0.

        if self.property == 'gas':
            vmin = -0.3
        if self.property == 'hotgas':
            vmin = 0.

            # if self.property == 'xrays':
            #   vmin = 0.

        fig = plt.figure(1, figsize=(20, 20))
        ax1 = fig.add_subplot(111)
        im1 = plt.imshow(img, extent=extent, origin='lower', cmap=self.cmap, vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.15)
        col = plt.colorbar(im1, cax=cax)
        col.set_label(self.axlab, fontsize=16)
        if save:
            plt.savefig(path + 'group' + str(self.groupnumber) + '_' + self.property + '_t%03d' % (theta) + '_p%03d' % (
            phi) + '.png')
        if show:
            plt.show()
        plt.close()


##################################################################################################################################

def load_array(quantity, ptype, array_type='SNAP', sim='L0100N1504', run='REFERENCE', tag='028_z000p000'):
    path = '/data5/simulations/EAGLE/' + sim + '/' + run + '/data'
    return np.array(E.readArray(array_type, path, tag, '/PartType' + str(ptype) + '/' + quantity),dtype=np.float32)


class wholebox(object):
    def __init__(self, prop,
                 sim='L0100N1504',
                 run='REFERENCE',
                 snapnum=28):
        print 'Initialising box for imaging...'

        tag = snapdict[str(snapnum)][0]

        sim_path = '/data5/simulations/EAGLE/' + sim + '/' + run + '/data/'

        # Get volume information
        boxsize = E.readAttribute('SNAP', sim_path, tag, "/Header/BoxSize")
        h = E.readAttribute('SNAP', sim_path, tag, "/Header/HubbleParam")
        a_0 = E.readAttribute('SNAP', sim_path, tag, "/Header/ExpansionFactor")

        # Point read_eagle to the data
        snapfile = sim_path + 'snapshot_' + tag + '/snap_' + tag + '.0.hdf5'
        comm = MPI.COMM_WORLD
        comm_rank = comm.Get_rank()
        comm_size = comm.Get_size()
        # Open snapshot
        snap = read.EagleSnapshot(snapfile)
        # Select region of interest
        snap.select_region(0.,boxsize,0.,boxsize,0.,boxsize)
        # Split selection between processors
        # This assigns an equal number of hash cells to each processor.
        snap.split_selection(comm_rank,comm_size)

        if prop == 'stars':
            #pos = load_array('Coordinates', 4, sim=sim, run=run, tag=tag).T
            #smoothing_length = load_array('SmoothingLength', 4, sim=sim, run=run, tag=tag)
            #quantity = load_array('Mass', 4, sim=sim, run=run, tag=tag) * 1e10

            pos = snap.read_dataset(4,'Coordinates') * a_0/h
            pos = pos.T
            smoothing_length = snap.read_dataset(4, 'SmoothingLength') * a_0 / h
            quantity = snap.read_dataset(4, 'Mass') / h * 1e10
        else:
            #pos = load_array('Coordinates', 0, sim=sim, run=run, tag=tag).T
            #smoothing_length = load_array('SmoothingLength', 0, sim=sim, run=run, tag=tag)

            pos = snap.read_dataset(0, 'Coordinates') * a_0 / h
            print pos
            pos = pos.T
            smoothing_length = snap.read_dataset(0, 'SmoothingLength') * a_0 / h

            if prop == 'gas':
                quantity = snap.read_dataset(0, 'Mass') / h / 1e10
                print quantity
            elif prop == 'xrays':
                pids = snap.read_dataset(0, 'ParticleIDs')
                print 'Matching x-rays to particles'

                xray_data = h5.File('/data6/arijdav1/Lx_matching/'+sim+'_'+run+'/'+tag+'.hdf5','r')
                xrays = np.array(xray_data['Xray_luminosity']) / 1e30
                xray_pids = np.array(xray_data['ParticleIDs'])
                #match_sort = np.argsort(xray_pids)
                #xrays = xrays[match_sort]
                #xray_pids = xray_pids[match_sort]

                quantity = xrays[np.searchsorted(xray_pids,pids)]


            else:
                raise IOError('Plot options are "gas","stars" or "xrays"')

        N = len(quantity)

        pos *= 1e3  # convert to kpc
        smoothing_length *= 1e3

        print N

        Particles = sphviewer.Particles(pos, quantity, hsml=smoothing_length)

        print Particles.get_pos()
        print Particles.get_mass()
        print Particles.get_hsml()

        self.Scene = sphviewer.Scene(Particles)

        self.sim = sim
        self.run = run
        self.tag = tag
        self.property = prop
        self.boxsize = boxsize / h

    def get_xyz(self, gns):  # returns an array of the centre coordinates (in kpc) for a list of group numbers
        gn_inds = np.asarray(gns) - 1
        simpath = '/data5/simulations/EAGLE/' + self.sim + '/' + self.run + '/data/'
        num_subs = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, "/FOF/NumOfSubhalos"))
        first_subhalo = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, 'FOF/FirstSubhaloID')[num_subs > 0])
        subfind_centres = np.array(E.readArray('SUBFIND', simpath, self.tag, 'Subhalo/CentreOfPotential'))[
                          first_subhalo, :]
        return subfind_centres[gn_inds] * 1e3

    def image(self, gn, x, y, z,  # coordinates must be in kpc!
              theta=0.,
              phi=0.,
              extent=400.,
              vmax = None,
              vmin = None,
              resolution=1024,
              camera_distance='infinity',
              cmap = None,
              save=True,
              show=True,
              set_contrast = True,
              path='/home/arijdav1/Dropbox/phd/figures/sph_pics/',
              movieflag=False,
              framenumber=0):

        if os.path.exists(path) == False:
            os.makedirs(path)


        if camera_distance == 'infinity':
            px_size = 2*extent / float(resolution)  # in kpc
            px_size *= 1e3
            mass_units = '$M_{\odot}$ $\mathrm{pc}^{-2}$'
            lum_units = '$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{pc}^{-2}$'
        else:
            px_size = (90./float(resolution))*60. # in arcmin
            mass_units = '$M_{\odot}$ $\mathrm{arcmin}^{-2}$'
            lum_units = '$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{arcmin}^{-2}$'


        if self.property == 'gas':
            if cmap==None:
                cmap = 'afmhot'
            axlab = r'$\log\Sigma_g$ [%s]'%(mass_units)
        elif self.property == 'xrays':
            if cmap==None:
                cmap = 'inferno'
            axlab = r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [%s]'%(lum_units)
        elif self.property == 'stars':
            if cmap==None:
                cmap = 'bone'
            axlab = r'$\log\Sigma_*$ [%s]'%(mass_units)
        else:
            raise IOError('Plot options are "gas", "stars" or "xrays"')


        
        if camera_distance == 'infinity':
            self.Scene.update_camera(x=x, y=y, z=z, r='infinity', t=theta, p=phi,
                                     extent=[-extent, extent, -extent, extent], xsize=resolution, ysize=resolution)

        else:
            self.Scene.update_camera(x=x, y=y, z=z, r=camera_distance, t=theta, p=phi,
                                     xsize=resolution, ysize=resolution)

        Render = sphviewer.Render(self.Scene)
        #Render.set_logscale()
        extent = Render.get_extent()
        rendered_img = Render.get_image()
        img = deepcopy(rendered_img)
        img[img==0.] = 1e-10
        img = np.log10(img)

        if self.property == 'xrays':
            img += 30.

        if self.property == 'gas':
            img += 20.

        img -= np.log10(px_size ** 2)


        if camera_distance == 'infinity':
            # If you don't know the dynamic range, set the defaults and make them callable
            if set_contrast:
                if self.property == 'stars':
                    self.vmin = -1.5
                elif self.property == 'gas':
                    self.vmin = -1.5
                elif self.property == 'xrays':
                    self.vmin = 22.5

                self.vmax = np.amax(img[np.isfinite(img)])

            elif not set_contrast and vmax==None or vmin==None:
                print 'Please use set_contrast or specify a vmax and vmin'
                exit()

        else:
            if set_contrast:
                minim = np.amin(img[np.isfinite(img)])
                maxim = np.amax(img[np.isfinite(img)])
                self.vmax = maxim
                
                if self.property == 'stars':
                    self.vmin = np.percentile(np.linspace(minim,maxim,100),72)
                elif self.property == 'gas':
                    self.vmin = np.percentile(np.linspace(minim,maxim,100),25)
                elif self.property == 'xrays':
                    self.vmin = 22.5

            elif not set_contrast and vmax==None or vmin==None:
                print 'Please use set_contrast or specify a vmax and vmin'
                exit()

        if vmax == None:
            vmax = self.vmax
        if vmin == None:
            vmin = self.vmin

        inches = np.floor((resolution+500.)/100.)

        fig = plt.figure(1, figsize=(inches,inches))
        ax1 = fig.add_subplot(111)
        im1 = plt.imshow(img, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.15)
        col = plt.colorbar(im1, cax=cax)
        col.set_label(axlab, fontsize=16)
        if save:
            if movieflag:
                plt.savefig(path+'frame%s'%(str(framenumber))+'.png',bbox_inches='tight')

            else:

                plt.savefig(path + 'group' + str(gn) + '_' +\
                        self.property + '_ext%0d' % (extent[1]) +\
                        '_t%03d' % (theta) + '_p%03d' % (phi) + '.png',bbox_inches='tight')
        if show:
            plt.show()
        plt.close()


    def flythrough(self,stopping_points, extent,
                        movie_length = 30,
                        framerate = 24,
                        resolution = 1024):

        num_frames_per_stop = int(float(movie_length*framerate)/float(len(stopping_points)-1))


        xs = []
        ys = []
        zs = []
        for s in range(len(stopping_points)-1):
            xs.extend(np.linspace(stopping_points[s][0],stopping_points[s+1][0],num_frames_per_stop))
            ys.extend(np.linspace(stopping_points[s][1],stopping_points[s+1][1],num_frames_per_stop))
            zs.extend(np.linspace(stopping_points[s][2],stopping_points[s+1][2],num_frames_per_stop))
            
        print 'Making flythrough movie...'
        for i in tqdm(range(len(xs))):
            
            self.image(0,xs[i],ys[i],zs[i],extent = extent, resolution = resolution,
                    path = '/data6/arijdav1/flythrough_movies/'+self.sim+'_'+self.run+'/'+self.property+'/',
                    movieflag = True,
                    framenumber = i,
                    show=False)


    def rotate(self,centre, extent, 
                        rotate_by = 360.,
                        movie_length = 30,
                        framerate = 24,
                        resolution = 1024):

        num_frames = movie_length*framerate

        phi_values = np.linspace(0.,rotate_by,num_frames)
            
        print 'Making rotating movie...'
        for i in tqdm(range(len(phi_values))):
            if i == 0:
                set_contrast = True
                vmax = None
                vmin = None
            else:
                set_contrast = False
                vmax = self.vmax
                vmin = self.vmin
            
            self.image(0,centre[0],centre[1],centre[2],extent = extent, resolution = resolution, phi = phi_values[i],
                    path = '/data6/arijdav1/rotating_movies/'+self.sim+'_'+self.run+'/'+self.property+'/',
                    movieflag = True,
                    framenumber = i,
                    show=False,
                    set_contrast = set_contrast,
                    vmax = vmax,
                    vmin = vmin)


class region(object): # For loading in a region around your chosen galaxy
    def __init__(self, prop,
                 sim='L0025N0376',
                 run='REFERENCE',
                 snapnum=28):
        print 'Initialising box for imaging...'

        tag = snapdict[str(snapnum)][0]

        if prop == 'xrays':
            print 'Loading particle X-ray luminosities'
            xray_data = h5.File('/hpcdata7/arijdav1/Lx_matching/'+sim+'_'+run+'/'+tag+'.hdf5','r')
            self.xrays = np.array(xray_data['Xray_luminosity']) / 1e30
            self.xray_pids = np.array(xray_data['ParticleIDs'])

        self.sim_path = '/data5/simulations/EAGLE/' + sim + '/' + run + '/data/'

        # Get volume information
        boxsize = E.readAttribute('SNAP', self.sim_path, tag, "/Header/BoxSize")
        self.h = E.readAttribute('SNAP', self.sim_path, tag, "/Header/HubbleParam")
        self.a_0 = E.readAttribute('SNAP', self.sim_path, tag, "/Header/ExpansionFactor")

        self.sim = sim
        self.run = run
        self.tag = tag
        self.property = prop
        self.boxsize = boxsize * self.a_0/self.h
        self.snapnum = snapnum



    def get_xyz(self, gns):  # returns an array of the centre coordinates (in kpc) for a list of group numbers
        print 'Finding centre coordinates...'
        gn_inds = np.asarray(gns) - 1
        simpath = '/data5/simulations/EAGLE/' + self.sim + '/' + self.run + '/data/'
        num_subs = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, "/FOF/NumOfSubhalos"))
        first_subhalo = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, 'FOF/FirstSubhaloID')[num_subs > 0])
        subfind_centres = np.array(E.readArray('SUBFIND', simpath, self.tag, 'Subhalo/CentreOfPotential'))[first_subhalo, :]
        return subfind_centres[gn_inds] * 1e3


    def get_alignment(self,gns): # Relies on a galaxy having been saved by halo_isolator
        # gives an alignment vector in degrees, based on the angular momentum of the stars
        try:
            ngrps = len(gns)
        except TypeError:
            ngrps = 1
            
        angles = np.zeros((ngrps,3))

        for g in range(ngrps):

            if ngrps == 1:
                grnum = gns
            else:
                grnum = gns[g]

            try:
                data = h5.File(
                    '/hpcdata7/arijdav1/saved_regions/' + self.sim + '_' + self.run + '/' + self.tag + '/all/group' + str(
                        grnum) + '.hdf5', 'r')
            except IOError:
                print 'Group number not saved by halo isolator'
                exit()

            s_r200 = np.array(data['Volume/r200'])
            sc = np.array(data['Stars/Coordinates'])
            sv = np.array(data['Stars/Velocity'])


            sm = np.array(data['Stars/Mass'])
            s_r2 = np.einsum('...j,...j->...',sc,sc)
            mask = np.where(s_r2<s_r200**2)[0]

            sc = sc[mask] * 3.08567758149137e16 * 1e6 # convert coords to km as velocity is in km s-1
            sv = sv[mask]

            sm = sm.reshape((len(sm),1))

            J = sm * np.cross(sc,sv)
            Jtot = np.sum(J,axis=0)
            Jmag = np.sqrt(Jtot[0]**2 + Jtot[1]**2 + Jtot[2]**2)
            # Theta and phi are reversed from the normal definitions because the images are in x-y plane
            theta = np.arccos(Jtot[2]/Jmag) # in degrees
            phi = np.arcsin(Jtot[1]/(Jmag*np.sin(theta)))

            angles[g,0] = np.arctan(Jtot[2]/(Jmag*np.sin(theta)*np.sin(phi))) * 180./np.pi
            angles[g,1] = np.arctan(Jtot[2]/(Jmag*np.sin(theta)*np.cos(phi))) * 180./np.pi
            angles[g,2] = phi * 180./np.pi


        return angles



    def select(self,centre,region_size): # Region size in Mpc
        
        print 'Loading region...'
        code_centre = centre * self.h/(self.a_0*1e3) # convert to h-less comoving code units
        region_size *= self.h/self.a_0

        # Point read_eagle to the data
        snapfile = self.sim_path + 'snapshot_' + self.tag + '/snap_' + self.tag + '.0.hdf5'

        # Open snapshot
        snap = read.EagleSnapshot(snapfile)
        # Select region of interest
        snap.select_region(code_centre[0]-region_size/2.,
                            code_centre[0]+region_size/2.,
                            code_centre[1]-region_size/2.,
                            code_centre[1]+region_size/2.,
                            code_centre[2]-region_size/2.,
                            code_centre[2]+region_size/2.)

        if self.property == 'stars':
            pos = snap.read_dataset(4,'Coordinates') * self.a_0/self.h
            smoothing_length = snap.read_dataset(4, 'SmoothingLength') * self.a_0 / self.h
            quantity = snap.read_dataset(4, 'Mass') / self.h * 1e10
        else:
            pos = snap.read_dataset(0, 'Coordinates') * self.a_0 / self.h
            smoothing_length = snap.read_dataset(0, 'SmoothingLength') * self.a_0 / self.h

            if self.property == 'gas':
                quantity = snap.read_dataset(0, 'Mass') / self.h / 1e10
            elif self.property == 'xrays':
                pids = snap.read_dataset(0, 'ParticleIDs')
                quantity = self.xrays[np.searchsorted(self.xray_pids,pids)]
            else:
                raise IOError('Plot options are "gas","stars" or "xrays"')

        
        print 'Wrapping box...'
        region_size /= self.h/self.a_0
        centre_mpc = centre /1e3
        pos = ne.evaluate("pos-centre_mpc")
        pos[pos[:,0]<(-1.*self.boxsize/2.),0] += self.boxsize
        pos[pos[:,1]<(-1.*self.boxsize/2.),1] += self.boxsize
        pos[pos[:,2]<(-1.*self.boxsize/2.),2] += self.boxsize
        pos[pos[:,0]>self.boxsize/2.,0] -= self.boxsize
        pos[pos[:,1]>self.boxsize/2.,1] -= self.boxsize
        pos[pos[:,2]>self.boxsize/2.,2] -= self.boxsize
        pos = ne.evaluate("pos+centre_mpc")
        



        pos = pos.T

        N = len(quantity)

        pos *= 1e3  # convert to kpc
        smoothing_length *= 1e3

        print 'Creating scene...'
        Particles = sphviewer.Particles(pos, quantity, hsml=smoothing_length)
        self.Scene = sphviewer.Scene(Particles)

    

    def image(self, gn, centre,  # coordinates must be in kpc!
              theta=0.,
              phi=0.,
              align = None,
              extent=400.,
              vmax = None,
              vmin = None,
              resolution=1024,
              camera_distance='infinity',
              cmap = None,
              save=True,
              show=True,
              set_contrast = True,
              path='/home/arijdav1/Dropbox/phd/figures/sph_pics/',
              movieflag=False,
              circleradius = None,
              scalebar_size = None,
              scalebar_colour = None,
              framenumber=0,
              returnarray=False):

        if align != None: # make sure the alignment isn't ruined by too many inputs
            theta = 0.
            phi = 0.

        if os.path.exists(path) == False:
            os.makedirs(path)

        if camera_distance == 'infinity':
            px_size = 2*extent / float(resolution)  # in kpc
            px_size *= 1e3
            mass_units = '$M_{\odot}$ $\mathrm{pc}^{-2}$'
            lum_units = '$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{pc}^{-2}$'
        else:
            px_size = (90./float(resolution))*60. # in arcmin
            mass_units = '$M_{\odot}$ $\mathrm{arcmin}^{-2}$'
            lum_units = '$\mathrm{erg}$ $\mathrm{s}^{-1}$ $\mathrm{arcmin}^{-2}$'


        if self.property == 'gas':
            if cmap==None:
                cmap = 'afmhot'
            axlab = r'$\log\Sigma_g$ [%s]'%(mass_units)
        elif self.property == 'xrays':
            if cmap==None:
                cmap = 'inferno'
            axlab = r'$\log\Sigma_{X,\mathrm{0.5-2keV}}$ [%s]'%(lum_units)
        elif self.property == 'stars':
            if cmap==None:
                cmap = 'bone'
            axlab = r'$\log\Sigma_*$ [%s]'%(mass_units)
        else:
            raise IOError('Plot options are "gas", "stars" or "xrays"')


        
        if camera_distance == 'infinity':
            self.Scene.update_camera(x=centre[0], y=centre[1], z=centre[2], r='infinity', t=theta, p=phi,
                                     extent=[-extent, extent, -extent, extent], xsize=resolution, ysize=resolution)

        else:
            self.Scene.update_camera(x=x, y=y, z=z, r=camera_distance, t=theta, p=phi,
                                     xsize=resolution, ysize=resolution)

        
        if align == 'face' or align == 'edge':
            angles = self.get_alignment(gn)
            rot_x = angles[0][0]
            rot_y = angles[0][1]
            rot_z = angles[0][2]

            print rot_x
            print rot_y
            print rot_z

            #self.Scene.update_camera(p=-1.*(90.+phi_J))
            self.Scene.update_camera(t=rot_x,p=rot_y,roll=rot_z)
        

        print 'Rendering image...'
        Render = sphviewer.Render(self.Scene)
        #Render.set_logscale()
        extent = Render.get_extent()
        rendered_img = Render.get_image()


        if returnarray == True:
            return rendered_img 

        img = deepcopy(rendered_img)
        img[img==0.] = 1e-10
        img = np.log10(img)

        if self.property == 'xrays':
            img += 30.

        if self.property == 'gas':
            img += 20.

        img -= np.log10(px_size ** 2)


        if camera_distance == 'infinity':
            # If you don't know the dynamic range, set the defaults and make them callable
            if set_contrast:
                if self.property == 'stars':
                    self.vmin = -1.5
                elif self.property == 'gas':
                    self.vmin = -1.5
                elif self.property == 'xrays':
                    self.vmin = 22.5

                self.vmax = np.amax(img[np.isfinite(img)])

            elif not set_contrast and vmax==None or vmin==None:
                print 'Please use set_contrast or specify a vmax and vmin'
                exit()

        else:
            if set_contrast:
                minim = np.amin(img[np.isfinite(img)])
                maxim = np.amax(img[np.isfinite(img)])
                self.vmax = maxim
                
                if self.property == 'stars':
                    self.vmin = np.percentile(np.linspace(minim,maxim,100),72)
                elif self.property == 'gas':
                    self.vmin = np.percentile(np.linspace(minim,maxim,100),25)
                elif self.property == 'xrays':
                    self.vmin = 22.5

            elif not set_contrast and vmax==None or vmin==None:
                print 'Please use set_contrast or specify a vmax and vmin'
                exit()

        if vmax == None:
            vmax = self.vmax
        if vmin == None:
            vmin = self.vmin

        inches = np.floor((resolution+500.)/100.)

        fig = plt.figure(1, figsize=(inches,inches))
        ax1 = fig.add_subplot(111)
        im1 = plt.imshow(img, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax1)
        cax = divider.append_axes("right", size="5%", pad=0.15)
        col = plt.colorbar(im1, cax=cax)
        col.set_label(axlab, fontsize=16)

        if circleradius != None:
            ax1.add_artist(plt.Circle((0,0),circleradius,fc='none',edgecolor='green',lw=2,ls='--'))

        if scalebar_size != None:
            ob = AnchoredHScaleBar(size=scalebar_size, label="", ax=ax1, loc=4, frameon=False, pad=0.6,sep=4,color=scalebar_colour) 
            ax1.add_artist(ob)
            if scalebar_size>=1000.:
                sb_label = '%.0f Mpc'%(scalebar_size/1000.)
            else:
                sb_label = '%.0f kpc'%(scalebar_size)

            #ax1.text(0.63*extent[1],-0.82*extent[1],sb_label,color=scalebar_colour,fontsize=16)

        if save:
            if movieflag:
                plt.savefig(path+'frame%s'%(str(framenumber))+'.png',bbox_inches='tight')

            else:

                plt.savefig(path + self.sim + '_snap' + str(self.snapnum) + '_group' + str(gn) + '_' +\
                        self.property + '_ext%0d' % (extent[1]) +\
                        '_t%03d' % (theta) + '_p%03d' % (phi) + '.png',bbox_inches='tight',dpi=400)
        if show:
            plt.show()
        plt.close()




if __name__ == '__main__':
    '''
    stars = region('stars')
    
    centre = stars.get_xyz(14)
    stars.select(centre,28.)
    stars.image(14,centre,extent=2000,save=False,theta=90,phi=10,scalebar_size=1000,circleradius=311.25)
    '''
    gas = region('gas')
    
    centre = gas.get_xyz(14)
    gas.select(centre,28.)
    gas.image(14,centre,extent=2000,save=False,theta=90,phi=10,scalebar_size=1000,scalebar_colour='white',circleradius=311.25)


    '''
    # User guide DO NOT DELETE
    gas = region('gas',sim='L0100N1504')
    centre = gas.get_xyz(2012)
    gas.select(centre,20.)
    gas.image(2012,centre,extent=10000.,path='/data6/arijdav1/v_vs_t_images/')
    '''


