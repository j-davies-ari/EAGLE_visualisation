import numpy as np
import eagle_constants_and_units as c
reload(c)

import ctypes as ct
#import matplotlib
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import mpl_toolkits.axes_grid1 as axgrid
import ion_header as ionh
import string
#import time

import make_maps_opts_locs as ol
reload(ol)
import h5py
#import numbers as num # for instance checking

def findiontables(ion,z):
    # README in dir_iontab:
    # files are hdf5, contain ionisation fraction of a species for rho, T, z
    
    
    #### checks and setup
    
    if not ion in ol.ions:
        print("There will be an error somewhere: %s is not included or misspelled. \n" % ion)
    
    tablefilename = ol.dir_iontab %ion + '.hdf5'      
    tablefile = h5py.File(tablefilename, "r")
    logTK =   np.array(tablefile.get('logt'),dtype=np.float32)  
    lognHcm3 =   np.array(tablefile.get('logd'),dtype=np.float32)
    ziopts = np.array(tablefile.get('redshift'),dtype=np.float32) # monotonically incresing, first entry is zero
    balance_d_t_z = np.array(tablefile.get('ionbal'),dtype=np.float32)
    tablefile.close()

    if z < 0. and z > 1e-4:
        z = 0.0
        zind = 0 
        interp = False
        
    elif z in ziopts:
        # only need one table
        zind = np.argwhere(z == ziopts)
        interp = False
    
    elif z <= ziopts[-1]:
        # linear interpolation between two tables
        zind1 = np.sum(ziopts<z)-1
        zind2 = -sum(ziopts>z)
        interp = True
    else:
        print("Chosen z value requires extrapolation. This has not been implemented. \n") 
        
    
    #### read in the tables; interpolate tables in z if needed and possible
    
    if not interp:
       balance = np.squeeze(balance_d_t_z[:,:,zind]) # for some reason, extra dimensions are tacked on 
    
    if interp: #linear interpolation: 1./(a1-a0) * ( (a1-a)*f0 + (a-a0)*f1 )
        balance1 = balance_d_t_z[:,:,zind1]
        balance2 = balance_d_t_z[:,:,zind2]
        
        print("interpolating 2 emission tables")
        balance = 1./( ziopts[zind2] - ziopts[zind1]) * ( (ziopts[zind2]-z)*balance1 + (z-ziopts[zind1])*balance2 )

    
    return balance, logTK, lognHcm3

def find_ionbal(z,ion,lognH,logT):
    
    # compared to the line emission files, the order of the nH, T indices in the balance tables is switched
    balance, logTK, lognHcm3 = findiontables(ion,z) #(np.array([[0.,0.],[0.,1.],[0.,2.]]), np.array([0.,1.,2.]), np.array([0.,1.]) ) 
    NumPart = len(lognH)
    inbalance = np.zeros(NumPart,dtype=np.float32)
    
    if len(logT) != NumPart:
        print('logrho and logT should have the same length')
        return None

    # need to compile with some extra options to get this to work: make -f make_emission_only
    print("------------------- C interpolation function output --------------------------\n")
    cfile = ol.c_interpfile

    acfile = ct.CDLL(cfile)
    interpfunction = acfile.interpolate_2d # just a linear interpolator; works for non-emission stuff too
    # ion balance tables are density x temperature x redshift 

    interpfunction.argtypes = [np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,)),\
                           ct.c_longlong , \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK)*len(lognHcm3),)), \
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(lognHcm3),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(len(logTK),)), \
                           ct.c_int,\
                           np.ctypeslib.ndpointer(dtype=ct.c_float, shape=(NumPart,))]

 
    res = interpfunction(lognH.astype(np.float32),\
               logT.astype(np.float32),\
               ct.c_longlong(NumPart),\
               np.ndarray.flatten(balance.astype(np.float32)),\
               lognHcm3.astype(np.float32),\
               ct.c_int(len(lognHcm3)),\
               logTK.astype(np.float32),\
               ct.c_int(len(logTK)), \
               inbalance \
              )

    print("-------------- C interpolation function output finished ----------------------\n")
    
    if res != 0:
        print('Something has gone wrong in the C function: output %s. \n',str(res))
        return None
        
    return inbalance

