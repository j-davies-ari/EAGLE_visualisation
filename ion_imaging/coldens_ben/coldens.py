import numpy as np
import eagle_constants_and_units as c
reload(c)
import ctypes as ct
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
import turtle

    # --------------------------------------------------------------------------------
    #
    #   NAME:
    #       COLDENS
    #
    #   PURPOSE:
    #       This function uses the SPH interpolation scheme to make column density
    #       maps of simulation output. It calls the C routines in the shared library
    #       'HsmlAndProject.so' to project the element or ion abundances of the SPH
    #       particles onto a grid and returns the 2D array containing the column
    #       density at each pixel position. It optionally plots the result as a
    #       colour-coded figure. Note that the projection is always done along the
    #       z axis, but it possibly to rotate the particle positions so as to
    #       essentially project along any given axis.
    #
    #   CALLING SEQUENCE:
    #       import coldens
    #       result = coldens.main(coords, l_smooth, mass_ion, centre, L_x, L_y, L_z,
    #                             npix_x, npix_y, desngb, BoxSize [, fig_name=None,
    #                             Vmin=None, Vmax=None, ion=None, theta=0.0, phi=0.0,
    #                             psi=0.0])
    #
    #   INPUTS:
    #       coords:   Array (N, 3) containing the particle positions [Mpc].
    #       l_smooth: SPH smoothing lengts [Mpc].
    #       mass_ion: Total mass [Msun] per SPH particle of the ion or element one
    #                 wants to make a column density map of.
    #       centre:   Coordinates [Mpc] of the centre of the map.
    #       L_x:      Size of the map [Mpc] in the horizontal direction.
    #       L_y:      Size of the map [Mpc] in the vertical direction.
    #       L_z:      Thickness of the slice [Mpc] one wants the project.
    #       npix_x:   Number of pixels in the horizontal direction.
    #       npix_y:   Number of pixels in the vertical direction.
    #       desngb:   Number of SPH neighbours to be used in the projection.
    #       BoxSize:  Size of the simulation box [Mpc].
    #
    #   OPTIONAL INPUTS:
    #       fig_name: Name of the (pdf) figure where the column density map should
    #                 be saved. If not specified, the result in not saved as a
    #                 colour-coded map.
    #       Vmin:     Minimum log column density [cm^-2] to be plotted. Default is
    #                 the minimum of ResultW.
    #       Vmax:     Maximum log column density [cm^-2] to be plotted. Default is
    #                 the maximum of ResultW.
    #       ion:      String containing the name of the ion, which will be used to
    #                 label the colour bar in the figure as 'N_ion'. If not
    #                 specified, the label will be just 'N'.
    #       theta:    Angle [degree] of rotation about the x axis. Positive
    #                 direction is from y to z.
    #       phi:      Angle [degree] of rotation about the y axis. Positive
    #                 direction is from z to x.
    #       psi:      Angle [degree] of rotation about the z axis. Positive
    #                 direction is from x to y.
    #
    #   OUTPUTS:
    #       ResultW:  Array (npix_x, npix_y) containing the log column density [cm^-2]
    #                 at each pixel position.
    #
    #   OPTIONAL OUTPUT:
    #       Figure showing ResultW as a colour-coded column density map.
    #
    # --------------------------------------------------------------------------------


def translate(old_coords, centre, boxsize):

    print 'Translating particle positions: (%.2f, %.2f, %.2f) -> (0, 0, 0) Mpc' \
          % (centre[0], centre[1], centre[2])

    # translates old coordinates into coordinates wrt centre
    # taking into account the periodicity of the box

    new_coords = np.ndarray((len(old_coords[:,0]), 3), dtype=np.float32)

    # x coordinates
    x_normal = np.where(np.abs(old_coords[:,0] - centre[0]) <= boxsize / 2.0)
    x_periodic_neg = np.where(old_coords[:,0] - centre[0] > boxsize / 2.0)
    x_periodic_pos = np.where(old_coords[:,0] - centre[0] < -1.0 * boxsize / 2.0)
    new_coords[x_normal[0],0] = old_coords[x_normal[0],0] - centre[0]
    new_coords[x_periodic_neg[0],0] = old_coords[x_periodic_neg[0],0] - centre[0] - boxsize
    new_coords[x_periodic_pos[0],0] = old_coords[x_periodic_pos[0],0] - centre[0] + boxsize

    # y coordinates
    x_normal = np.where(np.abs(old_coords[:,1] - centre[1]) <= boxsize / 2.0)
    x_periodic_neg = np.where(old_coords[:,1] - centre[1] > boxsize / 2.0)
    x_periodic_pos = np.where(old_coords[:,1] - centre[1] < -1.0 * boxsize / 2.0)
    new_coords[x_normal[0],1] = old_coords[x_normal[0],1] - centre[1]
    new_coords[x_periodic_neg[0],1] = old_coords[x_periodic_neg[0],1] - centre[1] - boxsize
    new_coords[x_periodic_pos[0],1] = old_coords[x_periodic_pos[0],1] - centre[1] + boxsize

    # z coordinates
    x_normal = np.where(np.abs(old_coords[:,2] - centre[2]) <= boxsize / 2.0)
    x_periodic_neg = np.where(old_coords[:,2] - centre[2] > boxsize / 2.0)
    x_periodic_pos = np.where(old_coords[:,2] - centre[2] < -1.0 * boxsize / 2.0)
    new_coords[x_normal[0],2] = old_coords[x_normal[0],2] - centre[2]
    new_coords[x_periodic_neg[0],2] = old_coords[x_periodic_neg[0],2] - centre[2] - boxsize
    new_coords[x_periodic_pos[0],2] = old_coords[x_periodic_pos[0],2] - centre[2] + boxsize

    return new_coords


def rotate(old_coords, theta, phi, psi):

    print 'Rotating particles over angles of (%.1f, %.1f, %.1f) degrees about x, y and z axis' \
          % (theta, phi, psi)

    # from angles [degree] to angles [rad]
    theta_rad = theta * (np.pi / 180.0)
    phi_rad   = phi * (np.pi / 180.0)
    psi_rad   = psi * (np.pi / 180.0)

    # construct rotation matrix
    rot_matrix = np.matrix(np.zeros((3,3)))
    rot_matrix[0,0] =  np.cos(phi_rad) * np.cos(psi_rad)
    rot_matrix[0,1] = -np.cos(phi_rad) * np.sin(psi_rad)
    rot_matrix[0,2] =  np.sin(phi_rad)
    rot_matrix[1,0] =  np.cos(theta_rad) * np.sin(psi_rad) + np.sin(theta_rad) * np.sin(phi_rad) * np.cos(psi_rad)
    rot_matrix[1,1] =  np.cos(theta_rad) * np.cos(psi_rad) - np.sin(theta_rad) * np.sin(phi_rad) * np.sin(psi_rad)
    rot_matrix[1,2] = -np.sin(theta_rad) * np.cos(phi_rad)
    rot_matrix[2,0] =  np.sin(theta_rad) * np.sin(psi_rad) - np.cos(theta_rad) * np.sin(phi_rad) * np.cos(psi_rad)
    rot_matrix[2,1] =  np.sin(theta_rad) * np.cos(psi_rad) - np.cos(theta_rad) * np.sin(phi_rad) * np.sin(psi_rad)
    rot_matrix[2,2] =  np.cos(theta_rad) * np.cos(phi_rad)
    
    # perform rotations with a matrix operation
    old_matrix = np.matrix(old_coords)
    new_matrix = rot_matrix * old_matrix.T
    new_coords = np.ndarray((len(old_coords[:,0]), 3), dtype=np.float32)
    new_coords[:,0] = new_matrix[0,:]    
    new_coords[:,1] = new_matrix[1,:]
    new_coords[:,2] = new_matrix[2,:]

    return new_coords


def make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix, redshift=None, extralabel=None, haloinfostr=None, mhalo=None,docbar=None, R200=None,kpc=None):

    kpc = 1
    print '\n--- Making colour map ---\n'

    print 'Min value: log N = %.5e cm^-2' % (Vmin)
    print 'Max value: log N = %.5e cm^-2' % (Vmax)
    print 'File:', fig_name

    if kpc != None:
        Xmin *= 1000.
        Xmax *= 1000.
        Ymin *= 1000.
        Ymax *= 1000.
        
    labelcolor = "black"
    
    #fig = plt.figure(figsize = (6, 5))
    if docbar == None or int(docbar)>0 :  
        fig = plt.figure(figsize = (npix/80.+0.5, npix/80.))
    else:
        fig = plt.figure(figsize = (npix/80., npix/80.))

    ax = plt.subplot(111)
    if npix >= 550:
        labelfontsize = 32
        extrastrfontsize = 16
        adjustable_color_label_pad = 16 # was 12
    else:
        labelfontsize = 20 #24 #16
        extrastrfontsize = 13
        adjustable_color_label_pad = 16
        
        
    # axis properties
    ax.set_xlim(Xmin, Xmax)
    ax.set_ylim(Ymin, Ymax)
    if(kpc != 'None'):
        ax.set_xlabel(r'[kpc]',fontsize=14)
    else:
        ax.set_xlabel(r'[Mpc]',fontsize=14)        
    #ax.set_title(r'Before AGN Turns On',fontsize=16)
    #ax.set_ylabel(r'Y [Mpc]')
    ax.minorticks_on()


    # colour map
    if ion == "Dispersion" or ion == "sigma" or ion == "sig" :
        img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("YlGnBu"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
    else:
        if ion == "T" or ion == "temperature" or ion == "To6" or ion == "T(OVI)":
            img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("gnuplot2"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
        else:
            if Vmin < 0.0 and Vmax > 0.0: # Difference map
                img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("RdBu"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
            else:  # Density  
                if Vmin < 0.0 and Vmax <= 0.0: 
                    img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("gnuplot2_r"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
                else:
                    if Vmin > 16: # Hydrogen
                        img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("plasma"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
                    else:
                        if Vmax >= 20: # HI
                            img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("viridis"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
                            labelcolor = "white"
                        else:
                            img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("gist_ncar_r"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
                            ###img = ax.imshow(ResultW.T, extent=(Xmin,Xmax,Ymin,Ymax), cmap=cm.get_cmap("inferno"), origin='lower', interpolation='nearest', vmin=Vmin, vmax=Vmax)
                            labelcolor = "black"

    if R200 != None:
        print "Adding circle at R200= ", float(R200)
        circle1 = plt.Circle((0, 0), float(R200), color='gray', fill=False)
        ax.add_artist(circle1)
        
    if docbar == None or int(docbar)>0 :          
        # colour bar legend
        div = axgrid.make_axes_locatable(ax)
        cax = div.append_axes("right", size="5%", pad=-0.2, zorder=20)
        #cax.set_zorder(20)
        cbar = plt.colorbar(img, cax=cax)
        #cbar.set_zorder(20)
    #if Vmin < 0.0: # Difference map
    #cbar.ax.set_ylabel(r'log N$_{\mathrm{%s}}$(NEQ-Equ) [cm$^{-2}$]' % (ion))
        if ion != None:
            if ion == "vx" or ion == "vy" or ion == "vz" or ion == "vdotroverr" or ion == "Velocity" or ion == "Vel" :
                cbar.ax.set_xlabel(r'v [km/s]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                
            else:
                if ion == "Dispersion" or ion == "sigma" or ion == "sig" :
                    cbar.ax.set_xlabel(r'$\sigma$ [km/s]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                else:
                    if ((ion == "T") or (ion == "temperature")):
                        cbar.ax.set_xlabel(r'lg$T$ [K]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                    else:
                        if ion == "nH":
                            cbar.ax.set_xlabel(r'lg$n_{\mathrm{H}}$ [cm$^{-3}$]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                        else:
                            if ion == "Z":
                                cbar.ax.set_xlabel(r'lg$Z$',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                            else:
                                if ion == "To6" or ion == "T(OVI)":
                                    cbar.ax.set_xlabel(r'lg$T_{\mathrm{OVI}}$ [K]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                                else:
                                    if ion == "xOVI":
                                        cbar.ax.set_xlabel(r'lg$x_{\mathrm{OVI}}$',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                                    else:
                                        if ion == "xOVII":
                                            cbar.ax.set_xlabel(r'lg$x_{\mathrm{OVII}}$',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                                        else:
                                            if ion == "xOVIII":
                                                cbar.ax.set_xlabel(r'lg$x_{\mathrm{OVIII}}$',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                                            else:
                                                if(Vmin<0):
                                                    if(Vmax <= 0):
                                                        cbar.ax.set_xlabel(r'lg%s'%(extralabel),labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
                                                    else:
                                                        cbar.ax.set_xlabel(r'$\Delta$lg$N$ [cm$^{-2}$]',labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)                                #cbar.ax.set_xlabel(r'$\Delta$log N$_{\mathrm{%s}}$ [cm$^{-2}$]' % (extralabel))
                                                else:
                                                    cbar.ax.set_xlabel(r'lg$N_{\mathrm{%s}}$ [cm$^{-2}$]' % (extralabel),labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)
        else:
            cbar.ax.set_xlabel(r'%s'%(extralabel),labelpad=adjustable_color_label_pad,fontsize=extrastrfontsize)

    print 'labelfontsize =', labelfontsize
    if redshift != None:
        ax.text((Xmax-Xmin)*0.925+Xmin,(Ymax-Ymin)*0.95+Ymin,r'z = %4.2f'%(redshift),horizontalalignment='right', verticalalignment='top', fontsize=26, color=labelcolor)
    if mhalo != None:
        ax.text((Xmax-Xmin)*0.925+Xmin,(Ymax-Ymin)*0.95+Ymin,r'log M$_{200}$=%4.1f'%(mhalo),horizontalalignment='right', verticalalignment='top', fontsize=labelfontsize, color=labelcolor)
    if extralabel != None:
        ax.text((Xmax-Xmin)*0.075+Xmin,(Ymax-Ymin)*0.95+Ymin,r'%s'%(extralabel),horizontalalignment='left', verticalalignment='top', fontsize=26, color=labelcolor)
    if haloinfostr != None:
        ax.text((Xmax-Xmin)*0.5+Xmin,(Ymax-Ymin)*0.075+Ymin,r'%s'%(haloinfostr),fontsize=extrastrfontsize, horizontalalignment='center', verticalalignment='top',zorder=100, color=labelcolor)
    
    # save PDF figure
    plt.savefig(fig_name)
    
    # close figure
    plt.close(fig)


def main(coords, l_smooth, mass_ion, centre, L_x, L_y, L_z, npix_x, npix_y, desngb, BoxSize, fig_name=None, Vmin=None, Vmax=None, ion=None, redshift=None, mhalo=None, extralabel=None, haloinfostr=None, theta=0.0, phi=0.0, psi=0.0, docbar=None, emission=None, R200=None, kpc=None):

    if(kpc):
        l_smooth *= 1.e+03
        L_x *= 1.e+03
        L_y *= 1.e+03
        L_z *= 1.e+03
        BoxSize *= 1.e+03
        centre *= 1.e+03
        if(R200!= None):
            R200 *= 1.e+03

    if docbar == None: 
        docbar = 1

    # ==============================================
    # Input
    # ==============================================

    print '\n--- Starting column density script ---\n'

    print 'Number of particles:', len(mass_ion)
    print 'Min, max, median of x [Mpc]:', np.min(coords[:,0]), np.max(coords[:,0]), np.median(coords[:,0])
    print 'Min, max, median of y [Mpc]:', np.min(coords[:,1]), np.max(coords[:,1]), np.median(coords[:,1])
    print 'Min, max, median of z [Mpc]:', np.min(coords[:,2]), np.max(coords[:,2]), np.median(coords[:,2])
    print 'Min, max, median of smoothing lengths [Mpc]:', np.min(l_smooth), np.max(l_smooth), np.median(l_smooth)
    print 'Min, max, median of ion mass [Msun]: %.5e %.5e %.5e' \
          % (np.min(mass_ion), np.max(mass_ion), np.median(mass_ion))
    print 'Centre of map: (%.2f, %.2f, %.2f) Mpc' % (centre[0], centre[1], centre[2])
    print 'Size of map: %.2f x %.2f Mpc' % (L_x, L_y)
    print 'Thickness of slice: %.2f Mpc' % (L_z)
    print 'Grid: %i x %i pixels' % (npix_x, npix_y)
    print 'Number of SPH neighbours:', desngb
    print 'Size of simulation box: %.2f Mpc' % (BoxSize)


    # ==============================================
    # Select particles, translate and rotate
    # ==============================================
    
    print '\n--- Selecting particles in cube ---\n'
    
    # select all particles within map region + maximum smoothing length
    max_Hsml = np.max(l_smooth)
    if(max_Hsml > 0.5): max_Hsml = 0.5 # BDO added because of spurious large particles in zooms, limit is 500 kpc.
    diff_x = np.minimum((coords[:,0] - centre[0] + BoxSize) % BoxSize,
                        (centre[0] - coords[:,0] + BoxSize) % BoxSize)
    diff_y = np.minimum((coords[:,1] - centre[1] + BoxSize) % BoxSize,
                        (centre[1] - coords[:,1] + BoxSize) % BoxSize)
    diff_z = np.minimum((coords[:,2] - centre[2] + BoxSize) % BoxSize,
                        (centre[2] - coords[:,2] + BoxSize) % BoxSize)
    select = ((diff_x < L_x / 2.0 + max_Hsml) & \
              (diff_y < L_y / 2.0 + max_Hsml) & \
              (diff_z < L_z / 2.0 + max_Hsml))
    coords = coords[select,:]
    NumPart = sum(select)
    
    print 'Number of particles selected:', NumPart

    # translate particle positions (with periodic boundary conditions)
    trans_coords = translate(coords, centre, BoxSize)

    # OPTIONAL: rotate particle positions
    if theta != 0.0 or phi != 0.0 or psi != 0.0:
        old_coords = trans_coords.copy()
        trans_coords = rotate(old_coords, theta, phi, psi)

    # positions [Mpc], kernel sizes [Mpc] and masses [Msun]
    pos = np.zeros((NumPart, 3)).astype(np.float32)
    pos[:,0] = trans_coords[:,0]
    pos[:,1] = trans_coords[:,1]
    pos[:,2] = trans_coords[:,2]
    Hsml = l_smooth[select].astype(np.float32)
    mass_ion = mass_ion[select].astype(np.float32)
    
    
    # ==============================================
    # Putting everything in right format for C routine
    # ==============================================
    
    print '\n--- Calling findHsmlAndProject ---\n'
    
    # define edges of the map wrt centre coordinates [Mpc]
    Xmin = -1.0 * L_x / 2.0
    Xmax = L_x / 2.0
    Ymin = -1.0 * L_y / 2.0
    Ymax = L_y / 2.0
    Zmin = -1.0 * L_z / 2.0
    Zmax = L_z / 2.0

     # define axes (projection always along Axis3, but you can rotate the particle position matrix)
    Axis1 = 0
    Axis2 = 1
    Axis3 = 2

    # maximum kernel size [Mpc]
    Hmax = 0.5 * L_y
    
    # arrays to be filled with resulting maps
    ResultW = np.zeros((npix_x, npix_y)).astype(np.float32)
    ResultQ = np.zeros((npix_x, npix_y)).astype(np.float32)
    
    # input arrays for C routine
    c_pos = pos[:,:]
    c_Hsml = Hsml[:]
    c_mass_ion = mass_ion[:]
    c_ResultW = ResultW[:,:]
    c_ResultQ = ResultQ[:,:]
    
    # check if HsmlAndProject changes 
    print 'Mass in: %.5e Msun' % (np.sum(c_mass_ion))
    
    # path to shared library
    #lib_path = '/home/oppenheimer/pythonCode/coldens_ben/HsmlAndProject/HsmlAndProject.so'
    #lib_path = '/home/extboppe/pythonCode/coldens_ben/HsmlAndProject_OMP/HsmlAndProject_v3_notree_gadget_omp.so'
    lib_path = '/home/arijdav1/Dropbox/phd/Code/sph_visualisation/coldens_ben/HsmlAndProject_OMP/HsmlAndProject_v3_notree_gadget_omp.so'
    #lib_path = '/cosma/home/dp004/dc-oppe1/pythonCode/coldens_ben/HsmlAndProject/HsmlAndProject.so'

    # load the library
    my_library = ct.CDLL(lib_path)
    
    # set the parameter types (numbers with ctypes, arrays with ndpointers)
    my_library.findHsmlAndProject.argtypes = [ct.c_int,
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,3)), 
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)), 
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)), 
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)), 
                                  ct.c_float,
                                  ct.c_float,
                                  ct.c_float,
                                  ct.c_float,
                                  ct.c_float,
                                  ct.c_float,
                                  ct.c_int,
                                  ct.c_int,
                                  ct.c_int,
                                  ct.c_int,
                                  ct.c_int,
                                  ct.c_int,
                                  ct.c_float,
                                  ct.c_double,
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(npix_x,npix_y)), 
                                  np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(npix_x,npix_y))]
    
    # set the return type
    my_library.findHsmlAndProject.restype = None
    
    print '----------'
    
    # call the findHsmlAndProject C routine
    my_library.findHsmlAndProject(ct.c_int(NumPart),    # number of particles in map
                                  c_pos,                # positions wrt to centre (NumPart, 3)
                                  c_Hsml,               # SPH kernel
                                  c_mass_ion,           # quantity to be mapped by projection (or weighting for average)
                                  c_mass_ion,           # quantity to be mapped by averaging
                                  ct.c_float(Xmin),     # left edge of map
                                  ct.c_float(Xmax),     # right edge of map
                                  ct.c_float(Ymin),     # bottom edge of map
                                  ct.c_float(Ymax),     # top edge of map
                                  ct.c_float(Zmin),     # near edge of map
                                  ct.c_float(Zmax),     # far edge of map
                                  ct.c_int(npix_x),     # number of pixels in x direction
                                  ct.c_int(npix_y),     # number of pixels in y direction
                                  ct.c_int(desngb),     # number of neightbours for SPH interpolation
                                  ct.c_int(Axis1),      # horizontal axis (x direction)
                                  ct.c_int(Axis2),      # vertical axis (y direction)
                                  ct.c_int(Axis3),      # projection axis (z direction)
                                  ct.c_float(Hmax),     # maximum allowed smoothing kernel
                                  ct.c_double(BoxSize), # size of simulation box
                                  c_ResultW,            # RESULT: map of projected mass (npix_x, npix_y)
                                  c_ResultQ)            # RESULT: map of mass weighted projected quantity (npix_x, npix_y)
    
    print '----------'
    
    # check if mass conserved (but first one also counts some particles outside the map actually)
    print 'Mass in:  %.5e Msun' % (np.sum(c_mass_ion))
    print 'Mass out: %.5e Msun' % (np.sum(ResultW))
    
    
    # ==============================================
    # Calculate column density and make map
    # ==============================================
    
    print '\n--- Calculating column density ---\n'
    
    # grid cell area [cm^2]
    area = (L_x / np.float32(npix_x)) * (L_y / np.float32(npix_y)) * c.cm_per_mpc ** 2
    
    # log10(column density [cm^-2])
    #factor = c.solar_mass / (c.atomw_C * c.u) / area  # BDO- was always divided by carbon mass!
    factor = c.solar_mass / (c.u) / area # BDO- No atomoic weight!!! Given in atomic units in line
    if emission != None:
        factor = 1./area
    
    print "ResultW= ", ResultW
    
    if Vmin < -50: # For cases of negative velocity
        ResultW = ResultW * factor # log10([cm^-2])
    else:
        ResultW = np.log10(ResultW) + np.log10(factor) # log10([cm^-2])

    print "area = ", area, " c.u= ", c.u, " c.solar_mass= ", c.solar_mass, "factor= ", factor
    print 'Min log column density: %.5e cm^-2' % (np.min(ResultW))
    print 'Max log column density: %.5e cm^-2' % (np.max(ResultW))
    print 'Mean log column density: %.5e cm^-2' % (np.mean(ResultW))
    print 'Median log column density: %.5e cm^-2' % (np.median(ResultW))

    # OPTIONAL: create colour map
    if fig_name != None:

        # minimum and maximum log column density to plot
        if Vmin == None:
            Vmin = np.min(ResultW)
        if Vmax == None:
            Vmax = np.max(ResultW)


        if mhalo != None:
            make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, extralabel=extralabel, mhalo=mhalo, docbar=docbar, R200=R200)
        else:
            if haloinfostr != None and redshift != None and extralabel != None:
                make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, redshift, extralabel, haloinfostr, docbar=docbar)
            else:
                if redshift != None and extralabel != None:
                    make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, redshift, extralabel, docbar=docbar)
                else:
                    if redshift != None:
                        make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, redshift, extralabel, docbar=docbar)
                    else:
                        if extralabel != None:
                            make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, extralabel=extralabel, docbar=docbar)
                        else:
                            make_colourmap(fig_name, ResultW, Xmin, Xmax, Ymin, Ymax, Vmin, Vmax, ion, npix_x, extralabel=extralabel, docbar=docbar)

    
    return ResultW
