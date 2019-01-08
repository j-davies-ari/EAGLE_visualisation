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

class region(object): # For loading in a region around your chosen galaxy
    def __init__(self, prop,
                 sim='L0025N0376',
                 run='REFERENCE',
                 snapnum=28,
                 tag = None,
                 quiet = False):

        if not quiet:
            print 'Initialising box for imaging...'

        if sim == 'L0100N1504':
            storage_loc = '/hpcdata6'
        else:
            storage_loc = '/hpcdata5'

        if tag == None:
            tag = snapdict[str(snapnum)][0]
        else:
            snapnum = int(tag[:3])

        if prop == 'xrays':
            if not quiet:
                print 'Loading particle X-ray luminosities'
            xray_data = h5.File('/hpcdata7/arijdav1/Lx_matching/'+sim+'_'+run+'/'+tag+'.hdf5','r')
            self.xrays = np.array(xray_data['Xray_luminosity']) / 1e30
            self.xray_pids = np.array(xray_data['ParticleIDs'])

        self.sim_path = storage_loc+'/simulations/EAGLE/' + sim + '/' + run + '/data/'

        # Get volume information
        boxsize = E.readAttribute('SNAP', self.sim_path, tag, "/Header/BoxSize")
        self.h = E.readAttribute('SNAP', self.sim_path, tag, "/Header/HubbleParam")
        self.a_0 = E.readAttribute('SNAP', self.sim_path, tag, "/Header/ExpansionFactor")

        self.z = (1./self.a_0) - 1.

        self.sim = sim
        self.run = run
        self.tag = tag
        self.property = prop
        self.boxsize = boxsize * self.a_0/self.h
        self.snapnum = snapnum
        self.quiet = quiet
        self.storage_loc = storage_loc


    def get_xyz(self, gns):  # returns an array of the centre coordinates (in kpc) for a list of group numbers
        if not self.quiet:
            print 'Finding centre coordinates...'
        gn_inds = np.asarray(gns) - 1
        #if self.run == 'REFERENCE_ApogeeRun':
        #    simpath = self.storage_loc+'/simulations/EAGLE/' + self.sim + '/' + self.run + '/data2/'
        #else:      
        simpath = self.storage_loc+'/simulations/EAGLE/' + self.sim + '/' + self.run + '/data/'
        #num_subs = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, "/FOF/NumOfSubhalos"))
        first_subhalo = np.array(E.readArray("SUBFIND_GROUP", simpath, self.tag, 'FOF/FirstSubhaloID'))
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
                data = h5.File('/hpcdata7/arijdav1/saved_regions/' + self.sim + '_' + self.run + '/' + self.tag + '/all/group' + str(grnum) + '.hdf5', 'r')
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
        if not self.quiet:
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

        if not self.quiet:
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
        if not self.quiet:
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
              cbar_loc = 'right',
              cbarpad = -80.,
              showaxes = False,
              save=True,
              show=True,
              set_contrast = True,
              path='/home/arijdav1/Dropbox/phd/figures/sph_pics/',
              fname = None,
              imageonly=False,
              circleradius = None,
              scalebar_size = None,
              scalebar_colour = None,
              redshift_label = False,
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
        
        if not self.quiet:
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

        if imageonly: # save the img array as a picture with no colourbars etc, at right resolution
            if fname == None:
                savepath = path + self.sim + '_snap%03d'%(self.snapnum) + '_group' + str(gn) + '_' +\
                        self.property + '_ext%0d' % (extent[1]) +\
                        '_t%03d' % (theta) + '_p%03d' % (phi) + '_noborder.png'
            else:
                savepath = path + fname




            #plt.imsave(fname=savepath,arr=img, vmin=vmin, vmax=vmax, cmap=cmap, origin='lower')

            fig = plt.figure(figsize=(1,1))
            ax = fig.add_axes([0,0,1,1])
            plt.tick_params(axis='x',which='both',bottom='off',top='off')
            plt.tick_params(axis='y',which='both',left='off',right='off')
            im1 = ax.imshow(img,vmin=vmin, vmax=vmax, cmap=cmap, origin='lower')
            if redshift_label:
                ax.text(0.95, 0.9, r'$z=%.1f$'%(self.z),
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax.transAxes,
                    color='white', fontsize=3)
            if circleradius != None:
                ax.add_artist(plt.Circle((0,0),circleradius,fc='none',edgecolor='w',lw=2,ls='--'))
            if scalebar_size != None:
                ob = AnchoredHScaleBar(size=scalebar_size, label="", ax=ax, loc=4, frameon=False, pad=0.6,sep=4,color=scalebar_colour) 
                ax.add_artist(ob)

            if not showaxes:
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)

            #divider = make_axes_locatable(ax)
            #cax = divider.append_axes("top",size="5%",pad=0.)
            #col = plt.colorbar(im1, cax=cax, orientation='horizontal')
            #col.ax.xaxis.set_ticks_position('top')
            #caxticks = col.ax.xaxis.get_major_ticks()
            #caxticks[0].label2.set_visible(False)
            #caxticks[-1].label2.set_visible(False)
            #col.set_label(axlab, labelpad=cbarpad, fontsize=16)

            plt.savefig(savepath,dpi=resolution)
            if show:
                plt.show()
            plt.close(fig)
            return

        #inches = np.floor((resolution+500.)/100.)
        #print inches

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        im1 = plt.imshow(img, extent=extent, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)
        divider = make_axes_locatable(ax1)
        if cbar_loc == 'top':
            cax = divider.append_axes("top",size="5%",pad=0.)
            col = plt.colorbar(im1, cax=cax, orientation='horizontal')
            col.ax.xaxis.set_ticks_position('top')
            caxticks = col.ax.xaxis.get_major_ticks()
            caxticks[0].label2.set_visible(False)
            caxticks[-1].label2.set_visible(False)
            col.set_label(axlab, labelpad=cbarpad, fontsize=24)
            cax.tick_params(axis='both', which='major', labelsize=16)
        else:
            cax = divider.append_axes("right", size="5%", pad=0.15)
            col = plt.colorbar(im1, cax=cax)
            col.set_label(axlab, fontsize=16)


        if not showaxes:
            ax1.axes.get_xaxis().set_visible(False)
            ax1.axes.get_yaxis().set_visible(False)


        if circleradius != None:
            ax1.add_artist(plt.Circle((0,0),circleradius,fc='none',edgecolor='w',lw=2,ls='--'))

        if scalebar_size != None:
            ob = AnchoredHScaleBar(size=scalebar_size, label="", ax=ax1, loc=4, frameon=False, pad=0.6,sep=4,color=scalebar_colour) 
            ax1.add_artist(ob)
            if scalebar_size>=1000.:
                sb_label = '%.0f Mpc'%(scalebar_size/1000.)
            else:
                sb_label = '%.0f kpc'%(scalebar_size)

            #ax1.text(0.63*extent[1],-0.82*extent[1],sb_label,color=scalebar_colour,fontsize=16)

        if save:

            DPI = fig.get_dpi()
            #DPI = 1024.
            fig.set_size_inches(1024.0/float(DPI),1080.0/float(DPI)) # for colorbar

            if fname == None:
                savepath = path + self.sim + '_snap%03d'%(self.snapnum) + '_group' + str(gn) + '_' +\
                        self.property + '_ext%0d' % (extent[1]) +\
                        '_t%03d' % (theta) + '_p%03d' % (phi) + '_noborder.png'
            else:
                savepath = path + fname

            plt.savefig(savepath,pad_inches=0.0,dpi='figure')
        if show:
            plt.show()
        plt.close()




if __name__ == '__main__':

    
    gas = region('gas',sim='L0100N1504')
    
    centre = gas.get_xyz(511)
    gas.select(centre,5.)
    gas.image(511,centre,extent=2000,save=False,theta=20,phi=200,scalebar_size=1000,scalebar_colour='white',circleradius=317.82544,cbar_loc='top')
    

    '''

    xrays = region('xrays',sim='L0025N0376',run='REFERENCE_ApogeeRun',tag='999_z000p000')

    for n in np.arange(9,100):
        print 'Group number ',n
        centre = xrays.get_xyz(n)
        xrays.select(centre,28.)
        xrays.image(n,centre,extent=5000.,save=False)
    
    '''

