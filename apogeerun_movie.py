import numpy as np
from eagle_imaging import region
import os
import glob
from sys import argv
import eagle as E
from tqdm import tqdm
from scipy.interpolate import interp1d
import centre_tracking

######################################
sim = 'L0025N0376'
run = 'REFERENCE_ApogeeRun'
prop = 'gas'
groupnum = 20
######################################

path = '/hpcdata7/arijdav1/sph_images/group'+str(groupnum)+'_movie/'+prop+'/'

if not os.path.exists(path):
    os.makedirs(path)

resolution = 1024

regionsize = 8. # full box size, in Mpc
extent = 2000. # half the box size, in kpc

root_dir = '/hpcdata5/simulations/EAGLE/L0025N0376/REFERENCE_ApogeeRun/data/'

simpaths = sorted(glob.glob(root_dir+'snapshot*'))
snaps = [f[-12:] for f in simpaths]


'''
# find the COP of the biggest galaxy at z=0
gas = region(prop,sim=sim,run=run,tag=snaps[-1],quiet=True)
centre = gas.get_xyz(groupnum) # centre on the biggest group


# Focus on a particular centre to capture AGN
gas = region(prop,sim=sim,run=run,tag='932_z000p111',quiet=True)
ae = 1./(1.+0.111)
centre = gas.get_xyz(27) / ae

'''

test_snaps = [0,249,499,749,999]

if 'setrange' in argv:
    # Make some test images for setting the dynamic range
    vmax_vals = []
    vmin_vals = []

    print 'Running dynamic range tests... '

    for ts in test_snaps:

        tag = snaps[ts]
        snapnum = int(tag[:3])
        a_exp = E.readAttribute('SNAP', root_dir, tag, "/Header/ExpansionFactor")

        gas = region(prop,sim=sim,run=run,tag=tag,quiet=True)
        gas.select(centre*a_exp,regionsize*a_exp)

        gas.image(groupnum,centre*a_exp,
                    extent=extent*a_exp,
                    resolution=resolution,
                    save=False,
                    show=True)


        while True:
            vmax_temp = np.float32(raw_input('Enter vmax: '))
            vmin_temp = np.float32(raw_input('Enter vmin: '))

            print 'Checking range... '
            gas.image(groupnum,centre*a_exp,
                        extent=extent*a_exp,
                        resolution=resolution,
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

    print vmax_vals
    print vmin_vals

else:
    if prop == 'gas':
        vmax_vals = [1.8, 2.2, 2.5, 1.6, 1.0]
        vmin_vals = [1.0, 0.0, -0.75, -1.3, -1.5]
    elif prop == 'xrays':
        vmax_vals = [30.0, 35.0, 35.0, 35.0, 32.0]
        vmin_vals = [20.0, 22.5, 22.5, 22.5, 22.5]

fmax = interp1d(test_snaps,vmax_vals)
fmin = interp1d(test_snaps,vmin_vals)

snapnumbers = np.arange(1000)
vmax_interpolated = fmax(snapnumbers)
vmin_interpolated = fmin(snapnumbers)


# Use the merger trees to follow the galaxy

tree = centre_tracking.mergertree()
tree.find_mainbranch(groupnum)
COP_phys = tree.centre_of_potential(comoving=False) * 1e3 # in kpc for imaging
#aexp = tree.expansion_factors()
tree_snapnums = tree.mainbranch_snapnums()

# Only do the later bit to catch the AGN
snaps = np.array(snaps)[tree_snapnums]
snapnumbers = np.array(snapnumbers)[tree_snapnums]
vmax_interpolated = np.array(vmax_interpolated)[tree_snapnums]
vmin_interpolated = np.array(vmin_interpolated)[tree_snapnums]

for n in tqdm(range(len(snaps))):
    
    #if snapnumbers[n] != 999:
    #    continue

    tag = snaps[n]

    snapnum = int(tag[:3])

    #print 'Snap ',snapnum

    a_exp = E.readAttribute('SNAP', root_dir, tag, "/Header/ExpansionFactor")

    filepath = path+sim+'_snap%03d'%(snapnum)+'_group'+str(groupnum)+'_'+prop+'_ext'+str(int(extent*a_exp))+'_t000_p000.png'

    if os.path.exists(filepath) and not 'overwrite' in argv:
        print 'already exists'
        continue 

    gas = region(prop,sim=sim,run=run,tag=tag,quiet=True)
    #centre = gas.get_xyz(1) # centre on the biggest group
    gas.select(COP_phys[n,:],regionsize*a_exp)
    gas.image(groupnum,COP_phys[n,:],
                extent=extent*a_exp,
                resolution=resolution,
                vmin = vmin_interpolated[n],
                vmax = vmax_interpolated[n],
                imageonly=False,
                path = path,
                redshift_label=True,
                show=False)







