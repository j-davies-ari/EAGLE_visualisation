import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
from sys import exit
from sys import argv
from tqdm import tqdm
import eagle as E
import readmodule # joel's tree reading module

# Hacked from halo isolator
def save(tree_dict, sim, run, vol_dir='/hpcdata7/arijdav1/mergertrees/', overwrite=False):

        vol_file = vol_dir + sim + '_' + run
        exist = os.path.exists(vol_file)

        # This is to prevent overwriting of existing data if you didn't check whether the files already exist
        if exist == True and overwrite == False:
            print 'This tree has already been converted. To overwrite, enter "overwrite" as a command line argument"'
            exit()

        # Generate the output file
        if exist == False or overwrite == True:
            if not os.path.exists(vol_dir):
                os.makedirs(vol_dir)
            f = h5py.File(vol_file + '.hdf5', 'w')

            for k, key in enumerate(tree_dict.keys()):
                f.create_dataset(key,data=tree_dict[key])

            f.close()



sim = 'L0025N0376'
run = 'REFERENCE_ApogeeRun'

tree_location = '/hpcdata5/simulations/EAGLE/L0025N0376/REFERENCE_ApogeeRun/mergerTree/mergerTree.dat'

data = readmodule.read()

print 'Loading...'
data.readfile(tree_location)

print 'Saving...'
save(data,sim,run)

print 'Done.'
