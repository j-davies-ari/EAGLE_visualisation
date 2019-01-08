import numpy as np
#import eagle_constants_and_units as c
##reload(c)

import ctypes as ct
#import matplotlib
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import mpl_toolkits.axes_grid1 as axgrid
#import ion_header as ionh
import string
#import time

#import make_maps_opts_locs as ol
#reload(ol)
import h5py

import eagle as E

simloc = '/hpcdata5/simulations/EAGLE/L0050N0752/REFERENCE/data/'

coords = np.array(E.readArray("SNAP", simloc, '028_z000p000', "/PartType0/Coordinates"))
mass = np.array(E.readArray("SNAP", simloc, '028_z000p000', "/PartType0/Mass"))

NumPart = len(mass)
Ls =

def project(NumPart,Ls,Axis1,Axis2,Axis3,box3,periodic,npix_x,npix_y,kernel,dct,tree,ompproj=True):
    '''
    dct must be a dictionary containing arrays 'coords', 'lsmooth', 'qW', 'qQ' (prevents copying of large arrays)
    '''

    # positions [Mpc / cm/s], kernel sizes [Mpc] and input quantities
    # a quirk of HsmlAndProject is that it only works for >= 100 particles. Pad with zeros if less.
    if NumPart >=100:
        pos = dct['coords'].astype(np.float32)
        Hsml = dct['lsmooth'].astype(np.float32)
        qW = dct['qW'].astype(np.float32)
        qQ = dct['qQ'].astype(np.float32)
    
    else:
        qQ = np.zeros((100,),dtype = np.float32)
        qQ[:NumPart] = dct['qQ'].astype(np.float32)  
        qW = np.zeros((100,),dtype = np.float32)
        qW[:NumPart] = dct['qW'].astype(np.float32)
        Hsml = np.zeros((100,),dtype = np.float32)
        Hsml[:NumPart] = dct['lsmooth'].astype(np.float32)
        pos = np.ones((100,3),dtype = np.float32)*1e8 #should put the particles outside any EAGLE projection region
        pos[:NumPart,:] = dct['coords'].astype(np.float32)
        NumPart = 100
    
    # ==============================================
    # Putting everything in right format for C routine
    # ==============================================
    
    print '\n--- Calling findHsmlAndProject ---\n'
    
    # define edges of the map wrt centre coordinates [Mpc]
    # in the periodic case, the C function expects all coordinates to be in the [0, BoxSize] range (though I don't think it actually reads Xmin etc. in for this)    
    # these need to be defined wrt the 'rotated' axes, e.g. Zmin, Zmax are always the min/max along the projection direction    
    if not periodic: # 0-centered
        Xmin = -1.0 * Ls[Axis1]/ 2.0
        Xmax = Ls[Axis1] / 2.0
        Ymin = -1.0 * Ls[Axis2] / 2.0
        Ymax = Ls[Axis2] / 2.0
        Zmin = -1.0 * Ls[Axis3] / 2.0
        Zmax = Ls[Axis3] / 2.0
  
    else: # half box centered (BoxSize used for x-y periodic boundary conditions)
        Xmin, Ymin = (0.,)*2
        Xmax,Ymax = (box3[Axis1],box3[Axis2])
        Zmin, Zmax = (0.5*(box3[Axis3] - Ls[Axis3]), 0.5*(box3[Axis3] + Ls[Axis3]))
    
    BoxSize = box3[Axis1]
        
    # maximum kernel size [Mpc] (modified from Marijke's version)
    Hmax = 0.5*min(Ls[Axis1],Ls[Axis2]) # Axis3 might be velocity; whole different units, so just ignore
        
    # arrays to be filled with resulting maps
    ResultW = np.zeros((npix_x, npix_y)).astype(np.float32)
    ResultQ = np.zeros((npix_x, npix_y)).astype(np.float32)
    
    # input arrays for C routine (change in c_pos <-> change in pos)
    c_pos = pos[:,:]
    c_Hsml = Hsml[:]
    c_QuantityW = qW[:]
    c_QuantityQ = qQ[:]
    c_ResultW = ResultW[:,:]
    c_ResultQ = ResultQ[:,:]
    
    # check if HsmlAndProject changes 
    print 'Total quantity W in: %.5e' % (np.sum(c_QuantityW))
    print 'Total quantity Q in: %.5e' % (np.sum(c_QuantityQ))
    
    # path to shared library
    if ompproj:
        sompproj = '_omp'
    else:
        sompproj = ''
    if tree:
        # in v3, projection can use more particles than c_int max,
        # but the tree building cannot
        if not ct.c_int(NumPart).value == NumPart:
            print(' ***         Warning         ***\n\nNumber of particles %i overflows C int type.\n This will likely cause the tree building routine in HsmlAndProjcet_v3 to fail.\nSee notes on v3 version.\n\n*****************************\n')
        if periodic:
            lib_path = ol.hsml_dir + 'HsmlAndProject_v3_%s_perbc%s.so' %(kernel, sompproj)       
        else:    
            lib_path = ol.hsml_dir + 'HsmlAndProject_v3_%s%s.so' %(kernel, sompproj)       
    else:
        if periodic:
            lib_path = ol.hsml_dir + 'HsmlAndProject_v3_notree_%s_perbc%s.so' %(kernel, sompproj)       
        else:    
            lib_path = ol.hsml_dir + 'HsmlAndProject_v3_notree_%s%s.so' %(kernel, sompproj)
    
    print('Using projection file: %s \n' % lib_path)
    # load the library
    my_library = ct.CDLL(lib_path)
    
    # set the parameter types (numbers with ctypes, arrays with ndpointers)
    my_library.findHsmlAndProject.argtypes = [ct.c_long,
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
    my_library.findHsmlAndProject(ct.c_long(NumPart),   # number of particles in map
                                  c_pos,                # positions wrt to centre (NumPart, 3)
                                  c_Hsml,               # SPH kernel
                                  c_QuantityW,          # quantity to be mapped by projection (or weighting for average)
                                  c_QuantityQ,          # quantity to be mapped by averaging
                                  ct.c_float(Xmin),     # left edge of map
                                  ct.c_float(Xmax),     # right edge of map
                                  ct.c_float(Ymin),     # bottom edge of map
                                  ct.c_float(Ymax),     # top edge of map
                                  ct.c_float(Zmin),     # near edge of map
                                  ct.c_float(Zmax),     # far edge of map
                                  ct.c_int(npix_x),     # number of pixels in x direction
                                  ct.c_int(npix_y),     # number of pixels in y direction
                                  ct.c_int(ol.desngb),  # number of neightbours for SPH interpolation
                                  ct.c_int(Axis1),      # horizontal axis (x direction)
                                  ct.c_int(Axis2),      # vertical axis (y direction)
                                  ct.c_int(Axis3),      # projection axis (z direction)
                                  ct.c_float(Hmax),     # maximum allowed smoothing kernel
                                  ct.c_double(BoxSize), # size of simulation box
                                  c_ResultW,            # RESULT: map of projected QuantityW (npix_x, npix_y)
                                  c_ResultQ)            # RESULT: map of projected QuantityQ weighted by QuantityW (npix_x, npix_y)
    
    print '----------'
    
    # check if mapped quantities conserved (but first one also counts some particles outside the map actually)
    print 'Total quantity W in:  %.5e' % (np.sum(c_QuantityW))
    print 'Total quantity W out: %.5e' % (np.sum(ResultW))
    print 'Total quantity Q in:  %.5e' % (np.sum(c_QuantityQ))
    print 'Total quantity Q out: %.5e' % (np.sum(ResultQ))

    return ResultW, ResultQ


