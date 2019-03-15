import sys
import tables
import eagle
#import h5py
import glob as glob
import os
import numpy as np
import scipy.stats as scist
import galmanip.readGalaxy as rg
import galmanip.rotateGalaxy as rotg
import galmanip.writeGalaxy as wg
import galmanip.binRadial as br
import galmanip.binSpherical as bs
import coldens_ben.coldens as coldens

#coords = [[1,1,1],[2,1.5,2],[3,-1.5,3],[1,0,3],[3,1,2]]
coords = ((1,1,1),(2,1.5,2),(3,-1.5,3),(1,0,3),(3,1,2))
x = (4,3,2,1,4)

vradial_hot_3D = bs.sphericalmean(x, coords,0,4,4,4,4)
print vradial_hot_3D
