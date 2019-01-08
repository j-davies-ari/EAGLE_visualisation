# -*- coding: utf-8 -*-

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
import read_eagle as read
from copy import deepcopy
import safe_colours
import numexpr as ne
import matplotlib.offsetbox
from matplotlib.lines import Line2D
safe_colours = safe_colours.initialise()

import metadata
from eagle_imaging_ions import region

import pickle

##################################################################
# Choose which simulation you want
simulation = 'L0100N1504'
run = 'REFERENCE'

mode = 'fixed' # fixed or virial

ion = 'c4'
im_vmin = 11.
im_vmax = 15.
extent = 250. # Radial extent from centre (ie half the frame width) in kpc
depth = 1000. # Depth of the image, sets what is LOADED IN, in kpc

# A list of redshifts that you want to save halos from
redshifts = [0.,]

# Change these directories to suit your system
im_save_path = '/hpcdata7/arijdav1/sph_images/ion_maps/'+simulation+'/'+run+'/'+ion+'/' # the location for the output images

map_save_path = '/hpcdata7/arijdav1/ion_maps/'+simulation+'/'+run+'/'+ion+'/' # the location for the output maps


data_location = '/hpcdata5/simulations/EAGLE/' # point to the EAGLE data - this is the LJMU directory

##################################################################

if not os.path.exists(im_save_path):
    os.makedirs(im_save_path)

if not os.path.exists(map_save_path):
    os.makedirs(map_save_path)

# Loops over redshifts
for k in range(len(redshifts)):
	print 'Running imaging for ',simulation,'-',run
	print 'Redshift: ',redshifts[k]

	# This uses my little 'metadata' module to find the snapshot tag for your chosen redshift
	tag =  metadata.snapnum_search(redshifts[k],returntag=True)
	sim = data_location+simulation+'/'+run+'/data/'

	# FOF quantities
	masslist = np.array(E.readArray("SUBFIND_GROUP",sim,tag,'FOF/Group_M_Crit200'))*1e10
	groupnums = np.arange(len(masslist)) + 1

	gns_for_loop = groupnums[masslist>np.power(10.,11.5)]

	gas = region('ion',ion=ion,sim=simulation)

	for g in tqdm(range(len(gns_for_loop))):

		centre = gas.get_xyz(gns_for_loop[g])
		gas.select(centre,depth/1e3)

		# Generate an image with vmin, vmax for reference
		gas.image(gns_for_loop[g],centre,
			extent=extent,
			save=True,
			show=False,
			theta=0,
			phi=0,
			set_contrast=False,
			vmin=im_vmin,
			vmax=im_vmax,
			showaxes=True,
			path = im_save_path,
			fname='group'+str(gns_for_loop[g])+'.png')


		# Get the true array of column densities, linear units in g cm^-2
		ion_map = gas.image(gns_for_loop[g],centre,
			returnarray=True,
			extent=extent,
			save=False,
			theta=0,
			phi=0,
			set_contrast=False,
			vmin=im_vmin,
			vmax=im_vmax,
			showaxes=True,
			path = im_save_path,
			fname='group'+str(gns_for_loop[g])+'.png')


		outfile = ion+'_map_group%04d.pkl'%(gns_for_loop[g])

		output = open(map_save_path+outfile, 'w')

		pickle.dump(ion_map,output)

		output.close()






