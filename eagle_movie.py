import numpy as np
#from eagle_imaging_ions import region
from eagle_imaging import region
import os
import glob
from sys import argv
import eagle as E
from tqdm import tqdm
from scipy.interpolate import interp1d
import centre_tracking
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


class movie(object):

    def __init__(self, groupnum, prop, moviename, # moviename is the folder in /hpcdata7/arijdav1/sph_images/ which will store everything
                        sim='L0025N0376',
                        run = 'REFERENCE_ApogeeRun',
                        resolution=1024, # resolution of the movie
                        extent=2000.,
                        ion=None): # side length of the movie in kpc

        print 'Initialising movie...'

        # Establish data location
        root_dir = '/hpcdata5/simulations/EAGLE/'+sim+'/'+run+'/data/'
        simpaths = sorted(glob.glob(root_dir+'snapshot*'))
        snaps = [f[-12:] for f in simpaths]

        # Set save destination
        path = '/hpcdata7/arijdav1/sph_images/'+moviename+'/'+prop+'/'
        if not os.path.exists(path):
            os.makedirs(path)

        regionsize = 8 * extent/1000. # how much to load in and image, in Mpc (default 8x the movie extent)

        self.groupnum = groupnum
        self.prop = prop
        self.moviename = moviename

        self.path = path
        self.sim = sim
        self.run = run
        self.resolution = resolution
        self.extent = extent
        self.regionsize = regionsize
        self.root_dir = root_dir
        self.snaps = snaps

        self.interpolation_snaps = [0,249,499,749,999]

        # Set these to none in case you want to run set_dynamic_range
        self.vmax_vals = None
        self.vmin_vals = None
        self.ion=ion

      
    def get_camera_positions(self):
        
        '''
        This function MUST be run first.
        Uses the merger tree to find the main branch of the group central we're interested in.
        Interpolates the centres of potential for smooth camera movement.
        Sets a (physical) camera_positions variable, as well as a variable for the SnapNums over which the galaxy exists.
        It also sets 5 SnapNums for interpolation of the dynamic range
        '''
        
        print 'Establishing camera positions...'
        tree = centre_tracking.mergertree()
        tree.find_mainbranch(self.groupnum)
        COP_comov = tree.centre_of_potential(comoving=True) * 1e3 # in kpc for imaging
        interp_comov_position = centre_tracking.interpolate_COP(COP_comov) # interpolate the comoving position
        self.tree_aexp = tree.expansion_factors()
        self.camera_positions = interp_comov_position * self.tree_aexp.reshape(len(self.tree_aexp),1) # then convert it to physical units
        tsns = tree.mainbranch_snapnums()
        L = len(tsns)
        self.interpolation_snaps = [tsns[0],tsns[int(L/4.)],tsns[int(L/2.)],tsns[int(3.*L/4.)],tsns[-1]]
        #self.tree_snapnums = tsns

        # Chop the movie down to the merger tree
        self.snaps = np.array(self.snaps)[tsns]
        self.snapnumbers = tsns


    def set_dynamic_range(self):  

        '''
        Creates a series of test images allowing you to set the dynamic range the movie interpolates over.
        Sets the vmin and vmax interpolation points at the end.
        '''

        vmax_vals = []
        vmin_vals = []

        print 'Running dynamic range tests... '

        for ts in self.interpolation_snaps:

            tag = self.snaps[ts]
            snapnum = int(tag[:3])
            a_exp = np.float(self.tree_aexp[self.snapnumbers==snapnum])
            camera_position = self.camera_positions[self.snapnumbers==snapnum]

            gas = region(self.prop,sim=self.sim,run=self.run,tag=tag,quiet=True,ion=self.ion)
            gas.select(camera_position[0],self.regionsize*a_exp)

            gas.image(self.groupnum,camera_position[0],
                        extent=self.extent*a_exp,
                        resolution=self.resolution,
                        save=False,
                        show=True)


            while True:
                vmax_temp = np.float32(raw_input('Enter vmax: '))
                vmin_temp = np.float32(raw_input('Enter vmin: '))

                print 'Checking range... '
                gas.image(self.groupnum,camera_position[0],
                            extent=self.extent*a_exp,
                            resolution=self.resolution,
                            vmin = vmin_temp,
                            vmax = vmax_temp,
                            save=False,
                            show=True)

                happy = raw_input('Happy with range? (y/n) ')
                if happy == 'y':
                    vmax_vals.append(vmax_temp)
                    vmin_vals.append(vmin_temp)
                    break
                else:
                    continue

        print 'Vmax interpolation values: ',vmax_vals
        print 'Vmin interpolation values: ',vmin_vals

        self.vmax_vals = vmax_vals
        self.vmin_vals = vmin_vals


    def generate_colour_scale(self):

        '''
        Interpolates the dynamic range across the whole movie, based either on the results of set_dynamic_range,
        or on the default values set in this function.
        '''

        print 'Establishing dynamic range...'
        if self.vmax_vals == None or self.vmin_vals == None:
            if self.prop == 'stars':
                self.vmax_vals = [3.2,3.2,3.2,3.2,3.2]
                self.vmin_vals = [0.,0.,0.,0.,0.]

            if self.prop == 'gas':
                self.vmax_vals = [1.8, 2.2, 2.5, 1.6, 1.0]
                self.vmin_vals = [1.0, 0.0, -0.75, -1.3, -1.5]

            elif self.prop == 'xrays':
                self.vmax_vals = [30.0, 35.0, 35.0, 35.0, 32.0]
                self.vmin_vals = [20.0, 22.5, 22.5, 22.5, 22.5]

            elif self.prop == 'entropy':
                self.vmax_vals = [-11.7, -10.8, -10.8, 35.0, 32.0]
                self.vmin_vals = [-14.5, -14.2, -14.2, 22.5, 22.5]

            elif self.prop == 'ion' and self.ion == 'c4':
                self.vmax_vals = [15., 15., 15., 15., 15.]
                self.vmin_vals = [11., 11., 11., 11., 11.]



        fmax = interp1d(self.interpolation_snaps,self.vmax_vals)
        fmin = interp1d(self.interpolation_snaps,self.vmin_vals)

        self.vmax_interpolated = fmax(self.snapnumbers)
        self.vmin_interpolated = fmin(self.snapnumbers)

    def make_movie(self,theta=0.,phi=0.,rotate_by=0.,imageonly=False,start_at=0):

        print 'Making movie...'

        num_frames = len(self.snaps)

        rotation = np.linspace(phi,phi+rotate_by,num_frames)

        for n in tqdm(range(num_frames)):
            
            if self.snapnumbers[n] < start_at:
                continue

            tag = self.snaps[n]
            snapnum = self.snapnumbers[n]
            a_exp = self.tree_aexp[n]

            fname = 'snap%03d.png'%(snapnum)


            if self.prop == 'ion' and self.snapnumbers[n] < 174:

                img = np.zeros((self.resolution,self.resolution))

                fig = plt.figure(figsize=(1,1))
                ax = fig.add_axes([0,0,1,1])
                plt.tick_params(axis='x',which='both',bottom='off',top='off')
                plt.tick_params(axis='y',which='both',left='off',right='off')
                im1 = plt.imshow(img, extent=[-1, 1, -1, 1], origin='lower', cmap='viridis',vmin=0.,vmax=1.)
                ax.axes.get_xaxis().set_visible(False)
                ax.axes.get_yaxis().set_visible(False)

                plt.savefig(self.path+fname,dpi=self.resolution)
                plt.close(fig)
                
                continue



            gas = region(self.prop,sim=self.sim,run=self.run,tag=tag,quiet=True,ion=self.ion)
            gas.select(self.camera_positions[n,:],self.regionsize*a_exp)
            gas.image(self.groupnum,
                        self.camera_positions[n,:],
                        theta=theta,
                        phi=rotation[n],
                        extent=self.extent*a_exp,
                        resolution=self.resolution,
                        vmin = self.vmin_interpolated[n],
                        vmax = self.vmax_interpolated[n],
                        imageonly=imageonly,
                        path = self.path,
                        fname = fname,
                        redshift_label=True,
                        show=False,
                        cbar_loc='top')

        self.final_phi = rotation[-1]
        self.final_theta = theta
        self.imageonly = imageonly

    def zoom_rotate_at_end(self,final_extent, num_zoom_frames = 100, num_rot_frames = 400):

        print 'Zooming in...'

        tag = self.snaps[-1]
        snapnum = self.snapnumbers[-1]
        a_exp = self.tree_aexp[-1]

        extent_list = np.linspace(self.extent*a_exp,final_extent*a_exp,num_zoom_frames)

        gas = region(self.prop,sim=self.sim,run=self.run,tag=tag,quiet=True,ion=self.ion)
        gas.select(self.camera_positions[-1,:],self.regionsize*a_exp)

        for n in tqdm(range(num_zoom_frames)):

            fname = 'snap%03d_zoom%03d.png'%(snapnum,n)

            gas.image(self.groupnum,
                        self.camera_positions[-1,:],
                        theta=self.final_theta,
                        phi=self.final_phi,
                        extent=extent_list[n],
                        resolution=self.resolution,
                        vmin = self.vmin_interpolated[-1],
                        vmax = self.vmax_interpolated[-1],
                        imageonly=self.imageonly,
                        path = self.path,
                        fname = fname,
                        redshift_label=True,
                        show=False,
                        cbar_loc='top')

        rotation = np.linspace(self.final_phi,self.final_phi+360.,num_rot_frames)        

        for n in tqdm(range(num_rot_frames)):

            fname = 'snap%03d_zoom%03d_rotate%03d.png'%(snapnum,num_zoom_frames-1,n)

            gas.image(self.groupnum,
                        self.camera_positions[-1,:],
                        theta=self.final_theta,
                        phi=rotation[n],
                        extent=extent_list[-1],
                        resolution=self.resolution,
                        vmin = self.vmin_interpolated[-1],
                        vmax = self.vmax_interpolated[-1],
                        imageonly=self.imageonly,
                        path = self.path,
                        fname = fname,
                        redshift_label=True,
                        cbar_loc='top',
                        show=False)

        

def main(grnum,image_property,initial_extent=1000.,final_extent=500.,start_at=0,set_dr=False,ion=None):

    mov = movie(grnum,image_property,'group%i_zoom_movie'%(grnum),extent=initial_extent)
    mov.get_camera_positions()

    if set_dr:
        mov.set_dynamic_range()
    #exit()

    mov.generate_colour_scale()
    mov.make_movie(phi=90.,rotate_by=90.,start_at=start_at,imageonly=True)
    mov.zoom_rotate_at_end(final_extent=final_extent)


if __name__ == '__main__':
    #main('ion',ion='c4')
    main(10,'stars',initial_extent=100.,final_extent=100.,set_dr=False,start_at=544)

