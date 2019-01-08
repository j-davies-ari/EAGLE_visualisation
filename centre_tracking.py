import numpy as np
import matplotlib.pyplot as plt
import h5py as h5
from scipy import interpolate

class mergertree(object):

    def __init__(self,sim = 'L0025N0376',run = 'REFERENCE_ApogeeRun'):

        tree_location = '/hpcdata7/arijdav1/mergertrees/'+sim+'_'+run+'.hdf5'

        treefile = h5.File(tree_location,'r')

        self.gn = np.array(treefile['GroupNumber'])
        self.sgn = np.array(treefile['SubGroupNumber'])
        self.snapnum = np.array(treefile['snapNum'])
        self.topleafID = np.array(treefile['mainLeafId'])
        self.haloID = np.array(treefile['haloId'])
        COP_x = np.array(treefile['CentreOfPotential_x']) # these are already in 
        COP_y = np.array(treefile['CentreOfPotential_y']) # physical units
        COP_z = np.array(treefile['CentreOfPotential_z'])
        redshift = np.array(treefile['redshift'])
        self.a_exp = 1./(1.+redshift)
        N_tot = len(self.gn)
        self.COP = np.hstack((COP_x.reshape(N_tot,1),COP_y.reshape(N_tot,1),COP_z.reshape(N_tot,1)))
        del COP_x, COP_y, COP_z

        self.tree = treefile

    def printkeys(self):
        print self.tree.keys()

    def find_mainbranch(self,gn_to_track,forwards=True):

        try:
            z0_idx = np.where((self.gn==gn_to_track)&(self.sgn==0)&(self.snapnum==999))[0][0]
        except IndexError:
            raise ValueError('Please enter a valid group number')

        mainbranch = np.arange(self.haloID[z0_idx],self.topleafID[z0_idx]+1)
        self.mainbranch_length = len(mainbranch)

        if forwards: # reverse the order of the indices so the tree advances forwards in time
            self.mainbranch = mainbranch[::-1]
        else:
            self.mainbranch = mainbranch

        # Check that branch is continuous (no gaps)
        if np.any(np.diff(self.snapnum[self.mainbranch])!=1):
            raise ValueError('Gap in main branch!')

    def centre_of_potential(self,comoving=False):

        COP_branch = self.COP[self.mainbranch,:]

        if comoving:
            return COP_branch / self.a_exp[self.mainbranch].reshape(self.mainbranch_length,1)
        else:
            return COP_branch

    def expansion_factors(self):

        return self.a_exp[self.mainbranch]

    def mainbranch_snapnums(self):
        
        return self.snapnum[self.mainbranch]

def interpolate_COP(COP_list,smoothness=10):
    indices = range(len(COP_list))
    indices_sliced = indices[::smoothness] # take every nth point
    indices_sliced.append(indices[-1]) # add final point to ensure correct interpolation

    xs = COP_list[indices_sliced,0]
    ys = COP_list[indices_sliced,1]
    zs = COP_list[indices_sliced,2]

    interp = interpolate.interp1d(indices_sliced,xs)
    x_interp = interp(indices)
    interp = interpolate.interp1d(indices_sliced,ys)
    y_interp = interp(indices)
    interp = interpolate.interp1d(indices_sliced,zs)
    z_interp = interp(indices)

    Nint = len(x_interp)

    return np.hstack((x_interp.reshape(Nint,1),y_interp.reshape(Nint,1),z_interp.reshape(Nint,1)))


        
        

if __name__ == '__main__':

    tree = mergertree()
    tree.find_mainbranch(1)
    COP_comov = tree.centre_of_potential(comoving=True)
    aexp = tree.expansion_factors()
    snapnums = tree.mainbranch_snapnums()

    COP_comov_interpolated = interpolate_COP(COP_comov)

    plt.figure()
    plt.scatter(COP_comov[:,0],COP_comov[:,1],c=aexp,edgecolor='None')
    plt.colorbar()
    plt.show()

    plt.figure()
    plt.scatter(COP_comov_interpolated[:,0],COP_comov_interpolated[:,1],c=aexp,edgecolor='None')
    plt.colorbar()
    plt.show()











