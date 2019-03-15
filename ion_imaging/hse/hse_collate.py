import sys
import tables
import eagle
#import h5py
import glob as glob
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import scipy.stats as scist
import galmanip.readGalaxy as rg
import galmanip.rotateGalaxy as rotg
import galmanip.writeGalaxy as wg
import galmanip.binRadial as br
import coldens_ben.coldens as coldens
from astropy import constants as const
import matplotlib.pyplot as plt

def init_list_of_objects(size):
    list_of_objects = list()
    for i in range(0,size):
        list_of_objects.append( list() )
    return list_of_objects

def calc_med_and_spread_1d(y):
    med50 = np.percentile(y,50)
    med25 = np.percentile(y,25)
    med75 = np.percentile(y,75)
    return(med50,med75,med25)

def calc_med_and_spread(x, y, nbins,low,hi):
    print "length-x, length-y ", len(x), len(y)
    med50 = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, 50), bins=nbins, range=[(low,hi)])[0]
    med75 = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, 75), bins=nbins, range=[(low,hi)])[0]
    med25 = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, 25), bins=nbins, range=[(low,hi)])[0]
    return(med50,med75,med25)

G_Grav = 6.674e-08
K_Boltz = 1.381e-16
M_P = 1.673e-24
C_S = 2.9989e+10
M_Solar = 1.989e+33
cmpermpc = 3.086e+24

omegaM = 0.307
omegaL = 1-omegaM
hubbleparam = 0.6777

nenh = 1.16
keVtoK = 8.6e-08

lsfile = sys.argv[1]
basedir = "/net/virgo/data5/oppenheimer/"
#basedir = "/net/galaxy/data2/oppenheimer/L100box/"

pversion = int(sys.argv[2])
###pversion = 1

extra_read_column = 0 #1

plotallhalobins = 0 #0

do_angmomaxis = 0 
do_hse = 0

cm = plt.get_cmap('jet') 
mhhi = 15.3
mhlow = 11.7

lMlowbin = 10.75
lMbinsize = 0.5
nhbins = 10
M200bins =  np.linspace(11.0,11.0+0.5*(nhbins-1),nhbins)
#colorbins = ['purple','blueviolet','blue','cyan','lime','gold','orange','red','magenta','brown']
colorbins = ['purple','blue','aqua','lime','orange','red','magenta','pink','tan','brown']
if(pversion):
    colorbins = ['darkblue','blue','aqua','blue','lime','green','violet','violet','brown','brown']
    colorbins = ['darkblue','blue','aqua','blue','darkviolet','green','lime','violet','brown','brown']
    samplename = ['','','','','','','','','','']

#lRfraclow = -2.0
#lRbinsize = 0.05

filename_in = basedir + lsfile
print filename_in 
file_in = open(filename_in,"r").readlines()
nhalos = len(file_in)

if(do_hse):
    fig_hse_all = plt.figure(figsize=(4.5,3.375))
    ax_hse_all = fig_hse_all.add_subplot(111)

fig_velDM_all = plt.figure(figsize=(4.5,3.375))
ax_velDM_all = fig_velDM_all.add_subplot(111)

fig_vradhot_all = plt.figure(figsize=(4.5,3.375))
ax_vradhot_all = fig_vradhot_all.add_subplot(111)

fig_vtanhot_all = plt.figure(figsize=(4.5,3.375))
ax_vtanhot_all = fig_vtanhot_all.add_subplot(111)

fig_sighot_all = plt.figure(figsize=(4.5,3.375))
ax_sighot_all = fig_sighot_all.add_subplot(111)

fig_gasfrac_all = plt.figure(figsize=(4.5,3.375))
ax_gasfrac_all = fig_gasfrac_all.add_subplot(111)

fig_DMfrac_all = plt.figure(figsize=(4.5,3.375))
ax_DMfrac_all = fig_DMfrac_all.add_subplot(111)

fig_DMfrac_all = plt.figure(figsize=(4.5,3.375))
ax_DMfrac_all = fig_DMfrac_all.add_subplot(111)

fig_Jtot_all = plt.figure(figsize=(4.5,3.375))
ax_Jtot_all = fig_Jtot_all.add_subplot(111)

fig_Thot_all = plt.figure(figsize=(4.5,3.375))
ax_Thot_all = fig_Thot_all.add_subplot(111)

fig_nHhot_all = plt.figure(figsize=(4.5,3.375))
ax_nHhot_all = fig_nHhot_all.add_subplot(111)

fig_gasfrac_coll = plt.figure(figsize=(4.5,3.375))
ax_gasfrac_coll = fig_gasfrac_coll.add_subplot(111)

fig_vradhot_coll = plt.figure(figsize=(4.5,3.375))
ax_vradhot_coll = fig_vradhot_coll.add_subplot(111)

fig_vtanhot_coll = plt.figure(figsize=(4.5,3.375))
ax_vtanhot_coll = fig_vtanhot_coll.add_subplot(111)

fig_sighot_coll = plt.figure(figsize=(4.5,3.375))
ax_sighot_coll = fig_sighot_coll.add_subplot(111)

fig_nHhot_coll = plt.figure(figsize=(4.5,3.375))
ax_nHhot_coll = fig_nHhot_coll.add_subplot(111)

fig_Thot_coll = plt.figure(figsize=(4.5,3.375))
ax_Thot_coll = fig_Thot_coll.add_subplot(111)

fig_Zhot_coll = plt.figure(figsize=(4.5,3.375))
ax_Zhot_coll = fig_Zhot_coll.add_subplot(111)

fig_Phot_coll = plt.figure(figsize=(4.5,3.375))
ax_Phot_coll = fig_Phot_coll.add_subplot(111)

fig_Shot_coll = plt.figure(figsize=(4.5,3.375))
ax_Shot_coll = fig_Shot_coll.add_subplot(111)

if(do_hse):
    fig_hse_coll = plt.figure(figsize=(4.5,3.375))
    ax_hse_coll = fig_hse_coll.add_subplot(111)
    
    fig_hsesum_coll = plt.figure(figsize=(4.5,3.375))
    ax_hsesum_coll = fig_hsesum_coll.add_subplot(111)
    
    fig_hsetherm_coll = plt.figure(figsize=(4.5,3.375))
    ax_hsetherm_coll = fig_hsetherm_coll.add_subplot(111)
    
    fig_hsestream_coll = plt.figure(figsize=(4.5,3.375))
    ax_hsestream_coll = fig_hsestream_coll.add_subplot(111)
    
    fig_hserot_coll = plt.figure(figsize=(4.5,3.375))
    ax_hserot_coll = fig_hserot_coll.add_subplot(111)
    
    fig_hseacc_coll = plt.figure(figsize=(4.5,3.375))
    ax_hseacc_coll = fig_hseacc_coll.add_subplot(111)

fig_mass_coll = plt.figure(figsize=(4.5,3.375))
ax_mass_coll = fig_mass_coll.add_subplot(111)

fig_E_coll = plt.figure(figsize=(4.5,3.375))
ax_E_coll = fig_E_coll.add_subplot(111)

fig_lambda_coll = plt.figure(figsize=(4.5,3.375))
ax_lambda_coll = fig_lambda_coll.add_subplot(111)

fig_Jspecific_coll = plt.figure(figsize=(4.5,3.375))
ax_Jspecific_coll = fig_Jspecific_coll.add_subplot(111)


lRfrac_hse_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
lMtot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mtotfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Msumfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mthermfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mrotfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mstreamfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Maccfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

lRfrac_mass_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
DMfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Barfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Gasfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Hotfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Coolfrac_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

vrad_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_DM_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

nH_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
T_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Z_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
P_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
S_hot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]


vrad_cool_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_cool_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_cool_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_cool_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

DM_cum_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fbar_cum_coll =  [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fhot_cum_coll =  [[-99 for i in xrange(1)] for i in xrange(nhbins)]


JDM_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jstars_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jgas_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jcool_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jhot_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

coolfracR200_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
hotfracR200_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
starcenfracR200_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Mtotfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mtotfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mtotfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Msumfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Msumfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Msumfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mthermfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mthermfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mthermfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mrotfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mrotfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mrotfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mstreamfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mstreamfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Mstreamfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Maccfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Maccfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Maccfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

DMfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
DMfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
DMfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Barfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Barfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Barfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Gasfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Gasfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Gasfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Hotfrac_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Hotfrac_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Hotfrac_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

vrad_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vrad_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vrad_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

vtan_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigrad_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigtan_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigtot_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigtot_DM_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_DM_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtot_DM_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

vrad_cool_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vrad_cool_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vrad_cool_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

vtan_cool_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_cool_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
vtan_cool_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigrad_cool_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_cool_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigrad_cool_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

sigtan_cool_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_cool_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
sigtan_cool_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

nH_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
nH_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
nH_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

T_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
T_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
T_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Z_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Z_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Z_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

P_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
P_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
P_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

S_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
S_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
S_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

mass_hot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
mass_hot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
mass_hot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

lMtot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
lMtot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
lMtot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

JDM_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
JDM_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
JDM_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Jstars_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jstars_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jstars_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Jgas_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jgas_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jgas_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Jcool_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jcool_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jcool_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

Jhot_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jhot_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
Jhot_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

DM_cum_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
DM_cum_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
DM_cum_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

fbar_cum_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fbar_cum_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fbar_cum_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]

fhot_cum_50 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fhot_cum_25 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
fhot_cum_75 = [[-99 for i in xrange(1)] for i in xrange(nhbins)]


#EkinR200_coll = [[-99 for i in xrange(1)] for i in xrange(nhbins)]
#EkinR200_coll = [[]]*nhbins
EkinR200_coll = init_list_of_objects(nhbins)
EthermR200_coll = init_list_of_objects(nhbins)
EgravR200_coll = init_list_of_objects(nhbins)
EBHR200_coll = init_list_of_objects(nhbins)
EstarsR200_coll = init_list_of_objects(nhbins)
M200_coll = init_list_of_objects(nhbins)

lambdaDM_coll = init_list_of_objects(nhbins)
lambdastars_coll = init_list_of_objects(nhbins)
lambdastars30_coll = init_list_of_objects(nhbins)
lambdagas_coll = init_list_of_objects(nhbins)
lambdahot_coll = init_list_of_objects(nhbins)
lambdacool_coll = init_list_of_objects(nhbins)

#EkinR200_coll = [None] * nhbins

print "EkinR200_coll= ", EkinR200_coll

f = open("CGM_Props.%s.dat"%lsfile,"w")
fion = open("CIV_Props.%s.dat"%lsfile,"w")


for i in range(nhalos):
    datadir = file_in[i].split()[0]
    galidstr = file_in[i].split()[1]
    substr = galidstr.split("_z")[1] 

    zname = substr.split("_")[0]
    zint = zname.split("p")[0]
    zfrac = zname.split("p")[1]
    z = float(zint)*1.0 + float(zfrac)/1000.
    ID = int(substr.split("_")[2])
    lM200 = float(substr.split("_")[3])
    M200 = 10**lM200
    Mstar = substr.split("_")[4]
    SFR = substr.split("_")[5]

    if(extra_read_column):
        if(pversion == 0):
            extracolumn = file_in[i].split()[2]
        else:
            extracolumn = file_in[i].split()[4]

    omegab = 0.04825
    XH = 0.75    
    proton_mass_cgs = 1.673e-24
    omegaratio = (omegaM+omegaL/(1+z)**3)
    R200 = 1.63e-5*(10**float(lM200)*hubbleparam)**0.333/omegaratio**0.333/(1+z)/hubbleparam
    lR200 = np.log10(R200)
    v200 = np.sqrt(G_Grav*10**lM200*M_Solar/(R200*cmpermpc))/1e+05
    lv200 = np.log10(v200)
    lTvir = 5.69 + 2/3.*(float(lM200)-12.0)+ 1/3.*np.log10(omegaratio) + np.log10(1+z)
    lrhocrit200 = np.log10(XH*1.88e-29*hubbleparam**2*(1+z)**3*omegab/omegaM*200/proton_mass_cgs)
   

    print "HALOINFO= ", substr, z, lM200, R200, v200

    e_filename = basedir + "/" + datadir + "/E_" + galidstr + ".dat"
    R, Ekin, Etherm, Egrav = np.loadtxt(e_filename, usecols=(0,1,2,3), unpack=True)
    k = 0
    while R[k] < lR200-0.0125:
        k = k + 1    

    EkinR200 = Ekin[k]
    EthermR200 = Etherm[k]
    EgravR200 = Egrav[k]

    masses_filename = basedir + "/" + datadir + "/masses_" + galidstr + ".dat"
    #Mstars, Mstarscen, MBH,MBH_max = np.loadtxt(masses_filename, usecols=(7,8,21,22), unpack=True)
    Mstars, Mstarscen, MBH,MBH_max = np.loadtxt(masses_filename, usecols=(7,8,21,21), unpack=True)

    Estars = 1.736e+49/0.55*Mstars #Crain+ 2015
    EBH = 0.015*(3e+10)**2*(MBH-1e+05/0.6777)*1.989e+33

    angmom_filename = basedir + "/" + datadir + "/angmom_" + galidstr + ".dat"
    R, lambda_DM, lambda_gas, lambda_stars, lambda_cool, lambda_hot, J_DM, J_gas, J_stars, J_cool, J_hot, Jtot_DM, Jtot_gas, Jtot_stars, Jtot_cool, Jtot_hot = np.loadtxt(angmom_filename, usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15), unpack=True)
    ###R, lambda_DM, lambda_gas, lambda_stars, lambda_cool, lambda_hot, J_DM, J_gas, J_stars, J_cool, J_hot, Jtot_DM, Jtot_gas, Jtot_stars, Jtot_cool, Jtot_hot = np.loadtxt(angmom_filename, usecols=(0,1,2,3,4,5,11,12,13,14,15,11,12,13,14,15), unpack=True)
    k = 0
    while R[k] < lR200-0.0125:
        k = k + 1    
    k30 = 0
    while R[k30] < np.log10(0.03)-0.0125:
        k30 = k30 + 1
    k50 = 0
    while R[k50] < np.log10(0.05)-0.0125:
        k50 = k50 + 1
    k70 = 0
    while R[k70] < np.log10(0.07)-0.0125:
        k70 = k70 + 1
    k100 = 0
    while R[k100] < np.log10(0.10)-0.0125:
        k100 = k100 + 1
        
    #print "R's: k, k30, k100 = ", k, k30, k100, 10**R[k], 10**R[k30], 10**R[k100], 10**lR200

    rbinhot_filename = basedir + "/" + datadir + "/rbins_hot_" + galidstr + ".dat"
    DM_cum, star_cum, gas_cum, hot_cum, T_hot, nH_hot, P_hot, S_hot, Z_hot = np.loadtxt(rbinhot_filename, usecols=(2,3,4,5,6,7,8,9,14), unpack=True)
    
    if(do_angmomaxis):
        angmomaxis_filename = basedir + "/" + datadir + "/angmomaxis_" + galidstr + ".dat"
        R_kpc, Mstar_kpc, anglestar_kpc, Mcool_kpc, anglecool_kpc, Mhot_kpc, anglehot_kpc, anglehotcool_kpc = np.loadtxt(angmomaxis_filename, usecols=(0,1,2,3,4,5,6,7), unpack=True)

        k30kpc = np.where(R_kpc == 30.0)
        k50kpc = np.where(R_kpc == 50.0)
        k70kpc = np.where(R_kpc == 70.0)
        k100kpc = np.where(R_kpc == 100.0)
        #print "Angmomaxis= ",k30kpc,anglestar_kpc[k30kpc],anglecool_kpc[k50kpc],anglehot_kpc[k50kpc]
    
    lambdaDM = lambda_DM[k]
    lambdastars = lambda_stars[k]
    lambdastars30 = lambda_stars[k30]
    lambdagas = lambda_gas[k]
    lambdacool = lambda_cool[k]
    lambdahot = lambda_hot[k]
    lRfrac = R - lR200

    do_ioncol = 0    
    if(do_ioncol):
        fion.write("%5d %5.2f %5.2f %6.3f"%(ID,float(lM200),float(Mstar),float(SFR)))
        fion.write(" %4.2f %5.3f %5.3f"%(np.log10(MBH_max),(gas_cum[k])/(DM_cum[k]+star_cum[k]+gas_cum[k])*0.307/0.04825,(gas_cum[k100])/(DM_cum[k100]+star_cum[k100]+gas_cum[k100])*0.307/0.04825))

        ions = ['h1','c4','o6']
        for j in range(len(ions)):
            ion_x_filename = basedir + "/" + datadir + "/columnave_lin." + galidstr + "." + ions[j] + ".x.l0.5.dat"
            kpclow,kpchi,ion_x_ave,ion_x_cov1, ion_x_cov2, ion_x_cov3, ion_x_n = np.loadtxt(ion_x_filename, usecols=(0,1,2,3,4,5,6), unpack=True)
        
            ion_y_filename = basedir + "/" + datadir + "/columnave_lin." + galidstr + "." + ions[j] + ".y.l0.5.dat"
            kpclow,kpchi,ion_y_ave,ion_y_cov1, ion_y_cov2, ion_y_cov3, ion_y_n = np.loadtxt(ion_y_filename, usecols=(0,1,2,3,4,5,6), unpack=True)

            ion_z_filename = basedir + "/" + datadir + "/columnave_lin." + galidstr + "." + ions[j] + ".z.l0.5.dat"
            kpclow,kpchi,ion_z_ave,ion_z_cov1, ion_z_cov2, ion_z_cov3, ion_z_n = np.loadtxt(ion_z_filename, usecols=(0,1,2,3,4,5,6), unpack=True)
            
            
            indexes_cov30kpc = np.where(kpchi <= 30.0)
            indexes_cov70kpc = np.where(kpchi <= 70.0)
            indexes_cov100kpc = np.where(kpchi <= 100.0)
            indexes_cov50kpc = np.where(kpchi <= 50.0)
            indexes_cov150kpc = np.where(kpchi <= 150.0)

            ion_mean_30kpc = np.log10(np.sum((ion_x_n[indexes_cov30kpc] * 10**ion_x_ave[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] * 10**ion_y_ave[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc] * 10**ion_z_ave[indexes_cov30kpc]))/np.sum((ion_x_n[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc])))
            ion_cov1_30kpc = np.sum((ion_x_n[indexes_cov30kpc] * ion_x_cov1[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] * ion_y_cov1[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc] * ion_z_cov1[indexes_cov30kpc]))/np.sum((ion_x_n[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc]))
            ion_cov2_30kpc = np.sum((ion_x_n[indexes_cov30kpc] * ion_x_cov2[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] * ion_y_cov2[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc] * ion_z_cov2[indexes_cov30kpc]))/np.sum((ion_x_n[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc]))
            ion_cov3_30kpc = np.sum((ion_x_n[indexes_cov30kpc] * ion_x_cov3[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] * ion_y_cov3[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc] * ion_z_cov3[indexes_cov30kpc]))/np.sum((ion_x_n[indexes_cov30kpc] + ion_y_n[indexes_cov30kpc] +  ion_z_n[indexes_cov30kpc]))

            ion_mean_50kpc = np.log10(np.sum((ion_x_n[indexes_cov50kpc] * 10**ion_x_ave[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] * 10**ion_y_ave[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc] * 10**ion_z_ave[indexes_cov50kpc]))/np.sum((ion_x_n[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc])))
            ion_cov1_50kpc = np.sum((ion_x_n[indexes_cov50kpc] * ion_x_cov1[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] * ion_y_cov1[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc] * ion_z_cov1[indexes_cov50kpc]))/np.sum((ion_x_n[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc]))
            ion_cov2_50kpc = np.sum((ion_x_n[indexes_cov50kpc] * ion_x_cov2[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] * ion_y_cov2[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc] * ion_z_cov2[indexes_cov50kpc]))/np.sum((ion_x_n[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc]))
            ion_cov3_50kpc = np.sum((ion_x_n[indexes_cov50kpc] * ion_x_cov3[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] * ion_y_cov3[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc] * ion_z_cov3[indexes_cov50kpc]))/np.sum((ion_x_n[indexes_cov50kpc] + ion_y_n[indexes_cov50kpc] +  ion_z_n[indexes_cov50kpc]))

            ion_mean_70kpc = np.log10(np.sum((ion_x_n[indexes_cov70kpc] * 10**ion_x_ave[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] * 10**ion_y_ave[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc] * 10**ion_z_ave[indexes_cov70kpc]))/np.sum((ion_x_n[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc])))
            ion_cov1_70kpc = np.sum((ion_x_n[indexes_cov70kpc] * ion_x_cov1[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] * ion_y_cov1[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc] * ion_z_cov1[indexes_cov70kpc]))/np.sum((ion_x_n[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc]))
            ion_cov2_70kpc = np.sum((ion_x_n[indexes_cov70kpc] * ion_x_cov2[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] * ion_y_cov2[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc] * ion_z_cov2[indexes_cov70kpc]))/np.sum((ion_x_n[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc]))
            ion_cov3_70kpc = np.sum((ion_x_n[indexes_cov70kpc] * ion_x_cov3[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] * ion_y_cov3[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc] * ion_z_cov3[indexes_cov70kpc]))/np.sum((ion_x_n[indexes_cov70kpc] + ion_y_n[indexes_cov70kpc] +  ion_z_n[indexes_cov70kpc]))

            ion_mean_100kpc = np.log10(np.sum((ion_x_n[indexes_cov100kpc] * 10**ion_x_ave[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] * 10**ion_y_ave[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc] * 10**ion_z_ave[indexes_cov100kpc]))/np.sum((ion_x_n[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc])))
            ion_cov1_100kpc = np.sum((ion_x_n[indexes_cov100kpc] * ion_x_cov1[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] * ion_y_cov1[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc] * ion_z_cov1[indexes_cov100kpc]))/np.sum((ion_x_n[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc]))
            ion_cov2_100kpc = np.sum((ion_x_n[indexes_cov100kpc] * ion_x_cov2[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] * ion_y_cov2[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc] * ion_z_cov2[indexes_cov100kpc]))/np.sum((ion_x_n[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc]))
            ion_cov3_100kpc = np.sum((ion_x_n[indexes_cov100kpc] * ion_x_cov3[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] * ion_y_cov3[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc] * ion_z_cov3[indexes_cov100kpc]))/np.sum((ion_x_n[indexes_cov100kpc] + ion_y_n[indexes_cov100kpc] +  ion_z_n[indexes_cov100kpc]))

            ion_mean_150kpc = np.log10(np.sum((ion_x_n[indexes_cov150kpc] * 10**ion_x_ave[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] * 10**ion_y_ave[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc] * 10**ion_z_ave[indexes_cov150kpc]))/np.sum((ion_x_n[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc])))
            ion_cov1_150kpc = np.sum((ion_x_n[indexes_cov150kpc] * ion_x_cov1[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] * ion_y_cov1[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc] * ion_z_cov1[indexes_cov150kpc]))/np.sum((ion_x_n[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc]))
            ion_cov2_150kpc = np.sum((ion_x_n[indexes_cov150kpc] * ion_x_cov2[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] * ion_y_cov2[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc] * ion_z_cov2[indexes_cov150kpc]))/np.sum((ion_x_n[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc]))
            ion_cov3_150kpc = np.sum((ion_x_n[indexes_cov150kpc] * ion_x_cov3[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] * ion_y_cov3[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc] * ion_z_cov3[indexes_cov150kpc]))/np.sum((ion_x_n[indexes_cov150kpc] + ion_y_n[indexes_cov150kpc] +  ion_z_n[indexes_cov150kpc]))


#            k25kpc_ion = np.where(kpchi == 30.0)
#            k50kpc_ion = np.where(kpchi == 50.0)
#            k75kpc_ion = np.where(kpchi == 80.0)
#            k100kpc_ion = np.where(kpchi == 100.0)
#            k150kpc_ion = np.where(kpchi == 150.0)
            
#            ion_med_cov1_25kpc = np.median([float(ion_x_cov1[k25kpc_ion]),float(ion_y_cov1[k25kpc_ion]),float(ion_z_cov1[k25kpc_ion])])
#            ion_med_cov1_50kpc = np.median([float(ion_x_cov1[k50kpc_ion]),float(ion_y_cov1[k50kpc_ion]),float(ion_z_cov1[k50kpc_ion])])
#            ion_med_cov1_75kpc = np.median([float(ion_x_cov1[k75kpc_ion]),float(ion_y_cov1[k75kpc_ion]),float(ion_z_cov1[k75kpc_ion])])
#            ion_med_cov1_100kpc = np.median([float(ion_x_cov1[k100kpc_ion]),float(ion_y_cov1[k100kpc_ion]),float(ion_z_cov1[k100kpc_ion])])
#            ion_med_cov1_150kpc = np.median([float(ion_x_cov1[k150kpc_ion]),float(ion_y_cov1[k150kpc_ion]),float(ion_z_cov1[k150kpc_ion])])
            
#            ion_med_cov2_25kpc = np.median([float(ion_x_cov2[k25kpc_ion]),float(ion_y_cov2[k25kpc_ion]),float(ion_z_cov2[k25kpc_ion])])
#            ion_med_cov2_50kpc = np.median([float(ion_x_cov2[k50kpc_ion]),float(ion_y_cov2[k50kpc_ion]),float(ion_z_cov2[k50kpc_ion])])
#            ion_med_cov2_75kpc = np.median([float(ion_x_cov2[k75kpc_ion]),float(ion_y_cov2[k75kpc_ion]),float(ion_z_cov2[k75kpc_ion])])
#            ion_med_cov2_100kpc = np.median([float(ion_x_cov2[k100kpc_ion]),float(ion_y_cov2[k100kpc_ion]),float(ion_z_cov2[k100kpc_ion])])
#            ion_med_cov2_150kpc = np.median([float(ion_x_cov2[k150kpc_ion]),float(ion_y_cov2[k150kpc_ion]),float(ion_z_cov2[k150kpc_ion])])
            
#            ion_med_cov3_25kpc = np.median([float(ion_x_cov3[k25kpc_ion]),float(ion_y_cov3[k25kpc_ion]),float(ion_z_cov3[k25kpc_ion])])
#            ion_med_cov3_50kpc = np.median([float(ion_x_cov3[k50kpc_ion]),float(ion_y_cov3[k50kpc_ion]),float(ion_z_cov3[k50kpc_ion])])
#            ion_med_cov3_75kpc = np.median([float(ion_x_cov3[k75kpc_ion]),float(ion_y_cov3[k75kpc_ion]),float(ion_z_cov3[k75kpc_ion])])
#            ion_med_cov3_100kpc = np.median([float(ion_x_cov3[k100kpc_ion]),float(ion_y_cov3[k100kpc_ion]),float(ion_z_cov3[k100kpc_ion])])
#            ion_med_cov3_150kpc = np.median([float(ion_x_cov3[k150kpc_ion]),float(ion_y_cov3[k150kpc_ion]),float(ion_z_cov3[k150kpc_ion])])

#            ion_med_ave_25kpc = np.median([float(ion_x_ave[k25kpc_ion]),float(ion_y_ave[k25kpc_ion]),float(ion_z_ave[k25kpc_ion])])
#            ion_med_ave_50kpc = np.median([float(ion_x_ave[k50kpc_ion]),float(ion_y_ave[k50kpc_ion]),float(ion_z_ave[k50kpc_ion])])
#            ion_med_ave_75kpc = np.median([float(ion_x_ave[k75kpc_ion]),float(ion_y_ave[k75kpc_ion]),float(ion_z_ave[k75kpc_ion])])
#            ion_med_ave_100kpc = np.median([float(ion_x_ave[k100kpc_ion]),float(ion_y_ave[k100kpc_ion]),float(ion_z_ave[k100kpc_ion])])
#            ion_med_ave_150kpc = np.median([float(ion_x_ave[k150kpc_ion]),float(ion_y_ave[k150kpc_ion]),float(ion_z_ave[k150kpc_ion])])
            
            
            fion.write("  %5.2f %5.2f %5.2f %5.2f %5.2f"%(ion_mean_30kpc,ion_mean_50kpc,ion_mean_70kpc,ion_mean_100kpc,ion_mean_150kpc))
            fion.write("  %5.3f %5.3f %5.3f %5.3f %5.3f"%(ion_cov1_30kpc,ion_cov1_50kpc,ion_cov1_70kpc,ion_cov1_100kpc,ion_cov1_150kpc))
            fion.write("  %5.3f %5.3f %5.3f %5.3f %5.3f"%(ion_cov2_30kpc,ion_cov2_50kpc,ion_cov2_70kpc,ion_cov2_100kpc,ion_cov2_150kpc))
            fion.write("  %5.3f %5.3f %5.3f %5.3f %5.3f"%(ion_cov3_30kpc,ion_cov3_50kpc,ion_cov3_70kpc,ion_cov3_100kpc,ion_cov3_150kpc))        

        fion.write("\n")
        
    if((np.isnan(lambdacool)) | (np.isinf(lambdacool))):
    #if((lambdacool <= 0.0) | (lambda_cool>1.0)):
        lambdacool = 0.0

    if((lM200>13.5) & (lM200<14.6)):
        ax_Jtot_all.plot(lRfrac, np.log10(Jtot_stars), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,alpha=0.5,ls='-')
        ax_Jtot_all.plot(lRfrac, np.log10(Jtot_cool), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=4,alpha=1.0,ls=':')
        ax_Jtot_all.plot(lRfrac, np.log10(Jtot_hot), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=4,alpha=1.0,ls='--')

    print "ENERGIES: ", lM200, EgravR200, EthermR200/EgravR200, EkinR200/EgravR200, Estars/EgravR200, EBH/EgravR200
    print "LAMBDAS: ", lM200, lambdaDM, lambdastars, lambdagas, lambdacool, lambdahot

    if(do_hse):
        hse_filename = basedir + "/" + datadir + "/hse_" + galidstr + ".dat"
        hse_filename = basedir + "/" + datadir + "/hse_sphere_" + galidstr + ".All.dat"

        if(extra_read_column):  hse_filename = basedir + "/" + datadir + "/hse_sphere_" + galidstr + ".All" + extracolumn + ".dat"
        print "hse_filename= ", hse_filename 
        R, Mtot, Mtherm, Mrot, Mstream, Macc, Msum = np.loadtxt(hse_filename, usecols=(0,1,2,3,4,5,6), unpack=True)
        #Macc = Macc - Mrot
        #Msum = Msum - Mrot
        #R, Mtot, Mtherm, Mrot, Mstream, Macc, Msum, fcov = np.loadtxt(hse_filename, usecols=(0,1,2,3,4,5,6,7), unpack=True)
        

        #Mrot = Mtandisp # REDEFINITION of Mrot to be Mtandisp
        lRfrac_hse = R - lR200
        lMtot = np.log10(Mtot)
        lMtherm = np.log10(Mtherm)
        lMrot = np.log10(Mrot)
        lMstream = np.log10(Mstream)
        lMacc = np.log10(Macc)    
        lMsum = np.log10(Msum)    

        ax_hse_all.plot(lRfrac_hse, lMtot-lM200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,alpha=0.5,ls='--')
        ax_hse_all.plot(lRfrac_hse, lMrot-lM200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2, ls='-')
        ax_hse_all.plot(lRfrac_hse, lMtherm-lM200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2, ls=':')
        ax_hse_all.plot(lRfrac_hse, lMacc-lM200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2, ls='-.')

    velDM_filename = basedir + "/" + datadir + "/vel_DM_" + galidstr + ".dat"
    R, sigrad_DM, vrad_DM, vtan_DM, sigtan_DM = np.loadtxt(velDM_filename, usecols=(0,1,2,3,4), unpack=True)
    lsigtot_DM = np.log10(np.sqrt(sigrad_DM**2 + sigtan_DM**2))
    lRfrac_velDM = R - lR200    
    lsigtot_DMfrac = lsigtot_DM-lv200
    
    ax_velDM_all.plot(lRfrac_velDM, sigtan_DM/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-', zorder=(lM200*10))

    velhot_filename = basedir + "/" + datadir + "/vel_hot_" + galidstr + ".dat"
    R, sigrad_hot, vrad_hot, vtan_hot, sigtan_hot = np.loadtxt(velhot_filename, usecols=(0,1,2,3,4), unpack=True)
    lsigtot_hot = np.log10(np.sqrt(sigrad_hot**2 + sigtan_hot**2))
    lRfrac_velhot = R - lR200    

    if(pversion==0):
        velcool_filename = basedir + "/" + datadir + "/vel_cool_" + galidstr + ".dat"
        R, sigrad_cool, vrad_cool, vtan_cool, sigtan_cool = np.loadtxt(velcool_filename, usecols=(0,1,2,3,4), unpack=True)
    else: # HACK- we'll not show cool for pversion==1.  
        R, sigrad_cool, vrad_cool, vtan_cool, sigtan_cool = np.loadtxt(velhot_filename, usecols=(0,1,2,3,4), unpack=True)
        
    #lsigtotfrac = lsigtot-lv200
    
    ax_vradhot_all.plot(lRfrac_velhot, vrad_hot/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-')
    ax_vradhot_all.plot(lRfrac_velhot, sigrad_hot/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls=':')
    ax_vtanhot_all.plot(lRfrac_velhot, vtan_hot/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-')
    ax_vtanhot_all.plot(lRfrac_velhot, sigtan_hot/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls=':')
    ax_sighot_all.plot(lRfrac_velhot, np.sqrt(sigrad_hot**2+sigtan_hot**2)/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-')
    ax_sighot_all.plot(lRfrac_velhot, np.sqrt(sigrad_DM**2+sigtan_DM**2)/v200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls=':')
    
    mass_filename = basedir + "/" + datadir + "/mass_" + galidstr + ".dat"
    if(do_angmomaxis):
        R, DM_mass, gas_mass, star_mass, cool_mass, hot_mass, coolnoISM_mass = np.loadtxt(mass_filename, usecols=(0,1,2,3,4,5,6), unpack=True)
    else:
        R, DM_mass, gas_mass, star_mass, cool_mass, hot_mass = np.loadtxt(mass_filename, usecols=(0,1,2,3,4,5), unpack=True)
        
    lRfrac_mass = R - lR200    
    lMDM = np.log10(DM_mass)
    lMstars = np.log10(star_mass)
    lMgas = np.log10(gas_mass)
    lMcool = np.log10(cool_mass)
    lMhot = np.log10(hot_mass)
    if(do_hse==0):
        lMtot = np.log10(DM_mass+star_mass+gas_mass)


    rot_over_tot_specmom_tan_cool = np.zeros(len(R))
    rot_over_tot_specmom_tan_hot = np.zeros(len(R))

    totv_coolsum = 0
    rotv_coolsum = 0
    totv_hotsum = 0
    rotv_hotsum = 0
    for j in range(len(R)):
        if(j>0):
            totv_coolsum += (cool_mass[j]-cool_mass[j-1])*sigtan_cool[j]*(10**R[j]+10**R[j-1])*1e+03/2.
            totv_hotsum += (hot_mass[j]-hot_mass[j-1])*sigtan_hot[j]*(10**R[j]+10**R[j-1])*1e+03/2.
            rotv_coolsum += (cool_mass[j]-cool_mass[j-1])*vtan_cool[j]*(10**R[j]+10**R[j-1])*1e+03/2.
            rotv_hotsum += (hot_mass[j]-hot_mass[j-1])*vtan_hot[j]*(10**R[j]+10**R[j-1])*1e+03/2.
            rot_over_tot_specmom_tan_cool[j] = rotv_coolsum/totv_coolsum
            rot_over_tot_specmom_tan_hot[j] = rotv_hotsum/totv_hotsum
            #totspecmom_tan_cool[j] = coolsum/cool_mass[j]
            #totspecmom_tan_hot[j] = hotsum/hot_mass[j]
            #print R[j], rot_over_tot_specmom_tan_cool[j], rot_over_tot_specmom_tan_hot[j]
            
            
    if(do_angmomaxis):
        f.write("%5d %5.2f %5.2f %6.3f"%(ID,float(lM200),float(Mstar),float(SFR)))
        f.write(" %5.1f"%anglecool_kpc[k50kpc])
        f.write(" %6.2f %5.3f  %6.1f %5.3f  %5.3f %5.3f  %5.1f"%(coolnoISM_mass[k30]/1e+08,hot_mass[k30]/(coolnoISM_mass[k30]+hot_mass[k30]),J_cool[k30],J_hot[k30]/J_cool[k30],rot_over_tot_specmom_tan_cool[k30],rot_over_tot_specmom_tan_hot[k30], anglehotcool_kpc[k30kpc]))
        f.write(" %6.2f %5.3f  %6.1f %5.3f  %5.3f %5.3f  %5.1f"%(coolnoISM_mass[k50]/1e+08,hot_mass[k50]/(coolnoISM_mass[k50]+hot_mass[k50]),J_cool[k50],J_hot[k50]/J_cool[k50],rot_over_tot_specmom_tan_cool[k50],rot_over_tot_specmom_tan_hot[k50], anglehotcool_kpc[k50kpc]))
        f.write(" %6.2f %5.3f  %6.1f %5.3f  %5.3f %5.3f  %5.1f"%(coolnoISM_mass[k70]/1e+08,hot_mass[k70]/(coolnoISM_mass[k70]+hot_mass[k70]),J_cool[k70],J_hot[k70]/J_cool[k70],rot_over_tot_specmom_tan_cool[k70],rot_over_tot_specmom_tan_hot[k70], anglehotcool_kpc[k70kpc]))
        f.write(" %6.2f %5.3f  %6.1f %5.3f  %5.3f %5.3f  %5.1f"%(coolnoISM_mass[k100]/1e+08,hot_mass[k100]/(coolnoISM_mass[k100]+hot_mass[k100]),J_cool[k100],J_hot[k100]/J_cool[k100],rot_over_tot_specmom_tan_cool[k100],rot_over_tot_specmom_tan_hot[k100], anglehotcool_kpc[k100kpc]))

        f.write("\n")
    

    #R_kpc, Mstar_kpc, anglestar_kpc, Mcool_kpc, anglecool_kpc, Mhot_kpc, anglehot_kpc, anglehotcool_kpc = np.loadtxt(angmomaxis_filename, usecols=(0,1,2,3,4,5,6,7), unpack=True)

    k = 0
    while R[k] < lR200-0.0125:
        k = k + 1    
    barfracR200 = (star_mass[k]+gas_mass[k])/M200*0.307/0.04825

    coolfracR200 = cool_mass[k]/gas_mass[k]
    hotfracR200 = hot_mass[k]/gas_mass[k]
    starcenfracR200 = Mstarscen/Mstars
    print "FRACTIONS: ", coolfracR200, hotfracR200, starcenfracR200

    #if(lM200>12.7): 
    ax_gasfrac_all.plot(lRfrac_mass, np.log10(gas_mass/(M200)*0.307/0.04825), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-')
    ax_gasfrac_all.plot(lRfrac_mass, np.log10(hot_mass/(M200)*0.307/0.04825), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls=':')

    #if(lM200>12.7): 
    ax_DMfrac_all.plot(lRfrac_mass, np.log10(DM_mass/(M200)*0.307/(0.307-0.04825)), color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=2,ls='-')

    ax_Thot_all.plot(lRfrac_velhot, np.log10(T_hot)-lTvir, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=1,ls='-',zorder=lM200*100)
    ax_nHhot_all.plot(lRfrac_velhot, np.log10(nH_hot)-lrhocrit200, color=cm((float(lM200)-mhlow)/(mhhi-mhlow)), lw=1,ls='-',zorder=lM200*100)


    if(plotallhalobins):
        hbin = int((lM200-lMlowbin)/lMbinsize)
    else:
        if((lM200>11.7) & (lM200<12.4)): hbin = 2
        if((lM200>12.7) & (lM200<=13.25)): hbin = 4
        if((lM200>13.25) & (lM200<13.85)): hbin = 6
        if((lM200>14.7) & (lM200<15.6)): hbin = 8

    if(pversion):
        hbin = int(file_in[i].split()[2])
        #print 'hbin = ', hbin
        samplename[hbin]  = file_in[i].split()[3]
        #samplename = '$M5.3$'
        #samplename[i] = 
        print 'samplename = ', samplename

    #print "lM200 = ", lM200, " hbin = ", hbin 
    accept = 1
    if(extra_read_column==0):
        if (np.abs(lM200 - lMDM[np.where(lRfrac_mass>-0.05)][0])>0.2):
            accept = 0
            print "Reject ", lM200,  ", wrong integrated mass = ", lM200, lMDM[np.where(lRfrac_mass>-0.05)][0], hbin
            
        if (lv200 - lsigtot_DM[np.where(lRfrac_mass>-0.30)][0] > 0.15):
            accept = 0
            print "Reject ", lM200,  ", wrong v200 at ~R500 = ",lv200,  lsigtot_DM[np.where(lRfrac_mass>-0.30)][0]

        if (lv200 - lsigtot_DM[np.where(lRfrac_mass>-1.05)][0] > 0.15):
            accept = 0
            print "Reject ", lM200,  ", wrong v200 at 0.1 R200 = ",lv200,  lsigtot_DM[np.where(lRfrac_mass>-1.05)][0]

    ###accept = 1
    #if((lv200 - lsigtot_DM[np.where(lRfrac_mass>-2.05)][0] > 0.15) & (lM200 < 14.7)): ###BDO taken out 8/12/17
        #accept = 0
        #print "Reject ", lM200,  ", wrong v200 at 0.01 R200 = ",lv200,  lsigtot_DM[np.where(lRfrac_mass>-2.05)][0]
       
    ###lsigtotfrac_coll[hbin].extend(lsigtotfrac)
    #if(lsigtotfrac[
    #print "hbin = ", hbin, lM200
    #rbin = (lRfrac-lRfraclow)/lRbinsize
    #for j in range(len(rbin)):
    if(accept):
        if(do_hse):
            lRfrac_hse_coll[hbin].extend(lRfrac_hse[np.where(np.isnan(Mrot)==0)])
            lMtot_coll[hbin].extend(lMtot[np.where(np.isnan(Mrot)==0)])
            Mtotfrac_coll[hbin].extend(Mtot[np.where(np.isnan(Mrot)==0)]/M200)
            Mthermfrac_coll[hbin].extend(Mtherm[np.where(np.isnan(Mrot)==0)]/Mtot[np.where(np.isnan(Mrot)==0)])
            Msumfrac_coll[hbin].extend(Msum[np.where(np.isnan(Mrot)==0)]/Mtot[np.where(np.isnan(Mrot)==0)])
            Mrotfrac_coll[hbin].extend(Mrot[np.where(np.isnan(Mrot)==0)]/Mtot[np.where(np.isnan(Mrot)==0)])
            Maccfrac_coll[hbin].extend(Macc[np.where(np.isnan(Mrot)==0)]/Mtot[np.where(np.isnan(Mrot)==0)])
            Mstreamfrac_coll[hbin].extend(Mstream[np.where(np.isnan(Mrot)==0)]/Mtot[np.where(np.isnan(Mrot)==0)])
        else:
            lMtot_coll[hbin].extend(lMtot[np.where(np.isnan(lMtot)==0)])
 

        vrad_hot_coll[hbin].extend(vrad_hot/v200)
        sigrad_hot_coll[hbin].extend(sigrad_hot/v200)
        vtan_hot_coll[hbin].extend(vtan_hot/v200)
        sigtan_hot_coll[hbin].extend(sigtan_hot/v200)
        sigtot_hot_coll[hbin].extend(np.sqrt(sigtan_hot**2+sigrad_hot**2)/v200)

        vrad_cool_coll[hbin].extend(vrad_cool/v200)
        sigrad_cool_coll[hbin].extend(sigrad_cool/v200)
        vtan_cool_coll[hbin].extend(vtan_cool/v200)
        sigtan_cool_coll[hbin].extend(sigtan_cool/v200)

        DM_cum_coll[hbin].extend(np.log10(DM_cum))
        fbar_cum_coll[hbin].extend((gas_cum+star_cum)/(DM_cum+star_cum+gas_cum))
        fhot_cum_coll[hbin].extend((hot_cum)/(DM_cum+star_cum+gas_cum))
        
        nH_hot_coll[hbin].extend(np.log10(nH_hot)-lrhocrit200)
        T_hot_coll[hbin].extend(np.log10(T_hot)-lTvir)
        Z_hot_coll[hbin].extend(np.log10(Z_hot))
        P_hot_coll[hbin].extend(np.log10(P_hot)-lTvir-lrhocrit200)
        S_hot_coll[hbin].extend(np.log10(S_hot)-lTvir+2./3*lrhocrit200-np.log10(keVtoK)+2/3.*np.log10(nenh))
        
        sigtot_DM_coll[hbin].extend(np.sqrt(sigtan_DM**2+sigrad_DM**2)/v200)

        lRfrac_mass_coll[hbin].extend(lRfrac_mass)
        DMfrac_coll[hbin].extend((DM_mass+star_mass+gas_mass)/10**lM200*(0.307)/0.04825)
        Barfrac_coll[hbin].extend((star_mass+gas_mass)/10**lM200*(0.307)/0.04825)
        Gasfrac_coll[hbin].extend(gas_mass/10**lM200*(0.307)/0.04825)
        Hotfrac_coll[hbin].extend(hot_mass/10**lM200*(0.307)/0.04825)
        Coolfrac_coll[hbin].extend(cool_mass/gas_mass)
        #print "gas_mass = ", gas_mass, 10**lM200

        JDM_coll[hbin].extend(J_DM/(R200*1e+03*v200))
        Jstars_coll[hbin].extend(J_stars/(R200*1e+03*v200))
        Jgas_coll[hbin].extend(J_gas/(R200*1e+03*v200))
        Jcool_coll[hbin].extend(J_cool/(R200*1e+03*v200))
        Jhot_coll[hbin].extend(J_hot/(R200*1e+03*v200))

        EkinR200_coll[hbin].append(EkinR200/EgravR200)
        EthermR200_coll[hbin].append(EthermR200/EgravR200)
        EBHR200_coll[hbin].append(EBH*barfracR200/EgravR200)
        EstarsR200_coll[hbin].append(Estars*barfracR200/EgravR200)
        EgravR200_coll[hbin].append(EgravR200)
        M200_coll[hbin].append(M200)

        lambdaDM_coll[hbin].append(lambdaDM)
        lambdastars_coll[hbin].append(lambdastars)
        lambdastars30_coll[hbin].append(lambdastars30)
        lambdagas_coll[hbin].append(lambdagas)
        lambdacool_coll[hbin].append(lambdacool)
        lambdahot_coll[hbin].append(lambdahot)

        coolfracR200_coll[hbin].append(coolfracR200)
        hotfracR200_coll[hbin].append(hotfracR200)
        starcenfracR200_coll[hbin].append(starcenfracR200)

rbinlow = -2.499
rbinhi = 0.501
nRbins = 31

#print "EkinR200_coll= ", EkinR200_coll

lRbins = np.linspace(rbinlow, rbinhi, num=nRbins)

for h in range(nhbins):

    if (len(lRfrac_mass_coll[h])>1):

        EgravR200_50, EgravR200_75, EgravR200_25 = calc_med_and_spread_1d(EgravR200_coll[h])
        EthermR200_50, EthermR200_75, EthermR200_25 = calc_med_and_spread_1d(EthermR200_coll[h])
        EkinR200_50, EkinR200_75, EkinR200_25 = calc_med_and_spread_1d(EkinR200_coll[h])
        EstarsR200_50, EstarsR200_75, EstarsR200_25 = calc_med_and_spread_1d(EstarsR200_coll[h])
        EBHR200_50, EBHR200_75, EBHR200_25 = calc_med_and_spread_1d(EBHR200_coll[h])
        M200_50, M200_75, M200_25 = calc_med_and_spread_1d(np.log10(M200_coll[h]))

        R200_50 = 1.63e-5*(10**(M200_50)*hubbleparam)**0.333/omegaratio**0.333/(1+z)/hubbleparam
        v200_50 = np.sqrt(G_Grav*10**M200_50*M_Solar/(R200_50*cmpermpc))/1e+05

        lambdaDM_50, lambdaDM_75, lambdaDM_25 = calc_med_and_spread_1d(lambdaDM_coll[h])
        lambdastars_50, lambdastars_75, lambdastars_25 = calc_med_and_spread_1d(lambdastars_coll[h])
        lambdastars30_50, lambdastars30_75, lambdastars30_25 = calc_med_and_spread_1d(lambdastars30_coll[h])
        lambdagas_50, lambdagas_75, lambdagas_25 = calc_med_and_spread_1d(lambdagas_coll[h])
        lambdacool_50, lambdacool_75, lambdacool_25 = calc_med_and_spread_1d(lambdacool_coll[h])
        lambdahot_50, lambdahot_75, lambdahot_25 = calc_med_and_spread_1d(lambdahot_coll[h])

        coolfracR200_50, coolfracR200_75, coolfracR200_25 = calc_med_and_spread_1d(coolfracR200_coll[h])
        hotfracR200_50, hotfracR200_75, hotfracR200_25 = calc_med_and_spread_1d(hotfracR200_coll[h])
        starcenfracR200_50, starcenfracR200_75, starcenfracR200_25 = calc_med_and_spread_1d(starcenfracR200_coll[h])

        #print "EkinR200_spread= ", M200_25, M200_50, M200_75, EkinR200_50, EkinR200_25, EkinR200_75
        #print "lambdacool(50,75,25)= ",lambdacool_50, lambdacool_75, lambdacool_25

        if(do_hse):
            Mtotfrac_50[h], Mtotfrac_75[h], Mtotfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Mtotfrac_coll[h],nRbins,rbinlow,rbinhi)
            Mrotfrac_50[h], Mrotfrac_75[h], Mrotfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Mrotfrac_coll[h],nRbins,rbinlow,rbinhi)
            Mthermfrac_50[h], Mthermfrac_75[h], Mthermfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Mthermfrac_coll[h],nRbins,rbinlow,rbinhi)
            Mstreamfrac_50[h], Mstreamfrac_75[h], Mstreamfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Mstreamfrac_coll[h],nRbins,rbinlow,rbinhi)
            Maccfrac_50[h], Maccfrac_75[h], Maccfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Maccfrac_coll[h],nRbins,rbinlow,rbinhi)
            Msumfrac_50[h], Msumfrac_75[h], Msumfrac_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], Msumfrac_coll[h],nRbins,rbinlow,rbinhi)
            lMtot_50[h], lMtot_75[h], lMtot_25[h] = calc_med_and_spread(lRfrac_hse_coll[h], lMtot_coll[h],nRbins,rbinlow,rbinhi)
        else:
            lMtot_50[h], lMtot_75[h], lMtot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], lMtot_coll[h],nRbins,rbinlow,rbinhi)
            
        DMfrac_50[h], DMfrac_75[h], DMfrac_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], DMfrac_coll[h],nRbins,rbinlow,rbinhi)
        Barfrac_50[h], Barfrac_75[h], Barfrac_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Barfrac_coll[h],nRbins,rbinlow,rbinhi)
        Gasfrac_50[h], Gasfrac_75[h], Gasfrac_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Gasfrac_coll[h],nRbins,rbinlow,rbinhi)
        Hotfrac_50[h], Hotfrac_75[h], Hotfrac_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Hotfrac_coll[h],nRbins,rbinlow,rbinhi)
        vrad_hot_50[h], vrad_hot_75[h], vrad_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], vrad_hot_coll[h],nRbins,rbinlow,rbinhi)
        vtan_hot_50[h], vtan_hot_75[h], vtan_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], vtan_hot_coll[h],nRbins,rbinlow,rbinhi)
        sigrad_hot_50[h], sigrad_hot_75[h], sigrad_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigrad_hot_coll[h],nRbins,rbinlow,rbinhi)
        sigtan_hot_50[h], sigtan_hot_75[h], sigtan_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigtan_hot_coll[h],nRbins,rbinlow,rbinhi)
        sigtot_hot_50[h], sigtot_hot_75[h], sigtot_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigtot_hot_coll[h],nRbins,rbinlow,rbinhi)
        sigtot_DM_50[h], sigtot_DM_75[h], sigtot_DM_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigtot_DM_coll[h],nRbins,rbinlow,rbinhi)

        vrad_cool_50[h], vrad_cool_75[h], vrad_cool_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], vrad_cool_coll[h],nRbins,rbinlow,rbinhi)
        vtan_cool_50[h], vtan_cool_75[h], vtan_cool_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], vtan_cool_coll[h],nRbins,rbinlow,rbinhi)
        sigrad_cool_50[h], sigrad_cool_75[h], sigrad_cool_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigrad_cool_coll[h],nRbins,rbinlow,rbinhi)
        sigtan_cool_50[h], sigtan_cool_75[h], sigtan_cool_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], sigtan_cool_coll[h],nRbins,rbinlow,rbinhi)

        DM_cum_50[h], DM_cum_75[h], DM_cum_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], DM_cum_coll[h],nRbins,rbinlow,rbinhi)
        fhot_cum_50[h], fhot_cum_75[h], fhot_cum_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], fhot_cum_coll[h],nRbins,rbinlow,rbinhi)
        fbar_cum_50[h], fbar_cum_75[h], fbar_cum_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], fbar_cum_coll[h],nRbins,rbinlow,rbinhi)

        nH_hot_50[h], nH_hot_75[h], nH_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], nH_hot_coll[h],nRbins,rbinlow,rbinhi)
        T_hot_50[h], T_hot_75[h], T_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], T_hot_coll[h],nRbins,rbinlow,rbinhi)
        Z_hot_50[h], Z_hot_75[h], Z_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Z_hot_coll[h],nRbins,rbinlow,rbinhi)
        P_hot_50[h], P_hot_75[h], P_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], P_hot_coll[h],nRbins,rbinlow,rbinhi)
        S_hot_50[h], S_hot_75[h], S_hot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], S_hot_coll[h],nRbins,rbinlow,rbinhi)
        
        JDM_50[h], JDM_75[h], JDM_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], JDM_coll[h],nRbins,rbinlow,rbinhi)#*(R200_50*1e+03*v200_50)
        Jstars_50[h], Jstars_75[h], Jstars_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Jstars_coll[h],nRbins,rbinlow,rbinhi)#*(R200_50*1e+03*v200_50)
        Jgas_50[h], Jgas_75[h], Jgas_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Jgas_coll[h],nRbins,rbinlow,rbinhi)#*(R200_50*1e+03*v200_50)
        Jcool_50[h], Jcool_75[h], Jcool_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Jcool_coll[h],nRbins,rbinlow,rbinhi)#*(R200_50*1e+03*v200_50)
        Jhot_50[h], Jhot_75[h], Jhot_25[h] = calc_med_and_spread(lRfrac_mass_coll[h], Jhot_coll[h],nRbins,rbinlow,rbinhi)#*(R200_50*1e+03*v200_50)
        JDM_50[h] *= (R200_50*1e+03*v200_50)
        JDM_25[h] *= (R200_50*1e+03*v200_50)
        JDM_75[h] *= (R200_50*1e+03*v200_50)
        Jstars_50[h] *= (R200_50*1e+03*v200_50)
        Jstars_25[h] *= (R200_50*1e+03*v200_50)
        Jstars_75[h] *= (R200_50*1e+03*v200_50)
        Jgas_50[h] *= (R200_50*1e+03*v200_50)
        Jgas_25[h] *= (R200_50*1e+03*v200_50)
        Jgas_75[h] *= (R200_50*1e+03*v200_50)
        Jcool_50[h] *= (R200_50*1e+03*v200_50)
        Jcool_25[h] *= (R200_50*1e+03*v200_50)
        Jcool_75[h] *= (R200_50*1e+03*v200_50)
        Jhot_50[h] *= (R200_50*1e+03*v200_50)
        Jhot_25[h] *= (R200_50*1e+03*v200_50)
        Jhot_75[h] *= (R200_50*1e+03*v200_50)

        fradial = open("radial.%s.%4.1f.dat"%(lsfile,M200_50),"w")

        for i in range(len(nH_hot_50[h])):
            if(np.isfinite(nH_hot_50[h][i])):
                #print "lRbins= ", h, lRbins
                #print "lMtot_50= ", lMtot_50[h]
                #print "nH_hot_50= ", nH_hot_50[h]                
                fradial.write("%6.3f  %5.3f % 5.3f % 5.3f % 5.3f % 5.3f % 5.3f  %6.3f %5.3f %5.3f\n"%(lRbins[i],lMtot_50[h][i],nH_hot_50[h][i],T_hot_50[h][i],Z_hot_50[h][i],P_hot_50[h][i],S_hot_50[h][i],DM_cum_50[h][i],fbar_cum_50[h][i],fhot_cum_50[h][i]))

        fradial.close()
        
        #print "Gas_frac[h] = ", h, Gasfrac_50[h], lMtot_50[h][60]
        #print "Mthermfrac[h] = ", h, Mthermfrac_50[h], lMtot_50[h][30]

        if(h%1==0):
            labelwrite = "log[$M_{\mathrm{200}}$]$=%6.3f$"%M200_50   #M200bins[h]
            if(pversion):
                labelwrite = samplename[h]

            #print 'labelwrite = ', labelwrite, ' colorbins = ', colorbins[h]
            lhotfrac = ax_gasfrac_coll.plot(10**lRbins, np.log10(Hotfrac_50[h]), color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            lgasfrac = ax_gasfrac_coll.plot(10**lRbins, np.log10(Gasfrac_50[h]), color='black', lw=4,alpha=1.0,ls=':',zorder=20-h)
            lbarfrac = ax_gasfrac_coll.plot(10**lRbins, np.log10(Barfrac_50[h]), color='black', lw=2,alpha=1.0,ls='--',zorder=20-h)
            ax_gasfrac_coll.plot(10**lRbins, np.log10(Hotfrac_50[h]), color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            #print "lRbins[30] = ", lRbins[30], " lhotgrac= ", Hotfrac_50[h][30], " lgasfrac = ",Gasfrac_50[h][30]
            if(pversion==0): ax_gasfrac_coll.fill_between(10**lRbins,np.log10(Hotfrac_25[h]),np.log10(Hotfrac_75[h]),facecolor=colorbins[h],alpha=0.2)
            ax_gasfrac_coll.plot(10**lRbins, np.log10(Gasfrac_50[h]), color=colorbins[h], lw=4,alpha=1.0,ls=':',zorder=20-h)
            ax_gasfrac_coll.plot(10**lRbins, np.log10(Barfrac_50[h]), color=colorbins[h], lw=2,alpha=1.0,ls='--',zorder=20-h)
            #ax_gasfrac_coll.plot(10**lRbins, np.log10(DMfrac_50[h]), color=colorbins[h], lw=2,alpha=0.2,ls='-.',zorder=20-h)

            #ax_gasfrac_coll.fill_between(10**lRbins,np.log10(Hotfrac_25[h]),np.log10(Hotfrac_75[h]),facecolor=colorbins[h],alpha=0.2)

            lvradhot = ax_vradhot_coll.plot(10**lRbins, vrad_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_vradhot_coll.plot(10**lRbins, vrad_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_vradhot_coll.fill_between(10**lRbins,vrad_hot_25[h],vrad_hot_75[h],facecolor=colorbins[h],alpha=0.2)
            if((h==2) & (pversion==0) & (plotallhalobins==0)):lvradcoolLstar = ax_vradhot_coll.plot(10**lRbins, vrad_cool_50[h], color=colorbins[h-1], lw=1,alpha=1.0,ls='-',zorder=20-h)


            lsigradhot = ax_vradhot_coll.plot(10**lRbins, sigrad_hot_50[h], color='black', lw=4,alpha=1.0,ls=':',zorder=20-h)
            ax_vradhot_coll.plot(10**lRbins, sigrad_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls=':',zorder=20-h)
            if(pversion==0): ax_vradhot_coll.fill_between(10**lRbins,sigrad_hot_25[h],sigrad_hot_75[h],facecolor=colorbins[h],alpha=0.2)
            if((h==2) & (pversion==0) & (plotallhalobins==0)):lsigradcoolLstar = ax_vradhot_coll.plot(10**lRbins, sigrad_cool_50[h], color=colorbins[h-1], lw=1,alpha=1.0,ls=':',zorder=20-h)

            lvtanhot = ax_vtanhot_coll.plot(10**lRbins, vtan_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_vtanhot_coll.plot(10**lRbins, vtan_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_vtanhot_coll.fill_between(10**lRbins,vtan_hot_25[h],vtan_hot_75[h],facecolor=colorbins[h],alpha=0.2)
            if((h==2) & (pversion==0) & (plotallhalobins==0)):lvtancoolLstar = ax_vtanhot_coll.plot(10**lRbins, vtan_cool_50[h], color=colorbins[h-1], lw=1,alpha=1.0,ls='-',zorder=20-h)

            lsigtanhot = ax_vtanhot_coll.plot(10**lRbins, sigtan_hot_50[h], color='black', lw=4,alpha=1.0,ls=':',zorder=20-h)
            ax_vtanhot_coll.plot(10**lRbins, sigtan_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls=':',zorder=20-h)
            if(pversion==0): ax_vtanhot_coll.fill_between(10**lRbins,sigtan_hot_25[h],sigtan_hot_75[h],facecolor=colorbins[h],alpha=0.2)
            if((h==2) & (pversion==0) & (plotallhalobins==0)):lsigtancoolLstar = ax_vtanhot_coll.plot(10**lRbins, sigtan_cool_50[h], color=colorbins[h-1], lw=1,alpha=1.0,ls=':',zorder=20-h)

            lsighot = ax_sighot_coll.plot(10**lRbins, sigtot_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_sighot_coll.plot(10**lRbins, sigtot_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_sighot_coll.fill_between(10**lRbins,sigtot_hot_25[h],sigtot_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            lsigDM = ax_sighot_coll.plot(10**lRbins, sigtot_DM_50[h], color='black', lw=4,alpha=1.0,ls='--',zorder=20-h)
            ax_sighot_coll.plot(10**lRbins, sigtot_DM_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='--',zorder=20-h)
            if(pversion==0): ax_sighot_coll.fill_between(10**lRbins,sigtot_DM_25[h],sigtot_DM_75[h],facecolor=colorbins[h],alpha=0.2)

            nHhot = ax_nHhot_coll.plot(10**lRbins, nH_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_nHhot_coll.plot(10**lRbins, nH_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_nHhot_coll.fill_between(10**lRbins,nH_hot_25[h],nH_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            Thot = ax_Thot_coll.plot(10**lRbins, T_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_Thot_coll.plot(10**lRbins, T_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_Thot_coll.fill_between(10**lRbins,T_hot_25[h],T_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            Zhot = ax_Zhot_coll.plot(10**lRbins, Z_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_Zhot_coll.plot(10**lRbins, Z_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_Zhot_coll.fill_between(10**lRbins,Z_hot_25[h],Z_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            Phot = ax_Phot_coll.plot(10**lRbins, P_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_Phot_coll.plot(10**lRbins, P_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_Phot_coll.fill_between(10**lRbins,P_hot_25[h],P_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            Shot = ax_Shot_coll.plot(10**lRbins, S_hot_50[h], color='black', lw=4,alpha=1.0,ls='-',zorder=20-h)
            ax_Shot_coll.plot(10**lRbins, S_hot_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
            if(pversion==0): ax_Shot_coll.fill_between(10**lRbins,S_hot_25[h],S_hot_75[h],facecolor=colorbins[h],alpha=0.2)

            if(do_hse):
                ax_hse_coll.plot(10**lRbins, Mthermfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hse_coll.fill_between(10**lRbins,Mthermfrac_25[h],Mthermfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hse_coll.plot(10**lRbins, Mrotfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='--',zorder=20-h)#,label=labelwrite)
                if(pversion==0): ax_hse_coll.fill_between(10**lRbins,Mrotfrac_25[h],Mrotfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hse_coll.plot(10**lRbins, Mstreamfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls=':',zorder=20-h)#,label=labelwrite)
                if(pversion==0): ax_hse_coll.fill_between(10**lRbins,Mstreamfrac_25[h],Mstreamfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hse_coll.plot(10**lRbins, Maccfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-.',zorder=20-h)#,label=labelwrite)
                if(pversion==0): ax_hse_coll.fill_between(10**lRbins,Maccfrac_25[h],Maccfrac_75[h],facecolor=colorbins[h],alpha=0.2)

                ax_hsetherm_coll.plot(10**lRbins, Mthermfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hsetherm_coll.fill_between(10**lRbins,Mthermfrac_25[h],Mthermfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hserot_coll.plot(10**lRbins, Mrotfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hserot_coll.fill_between(10**lRbins,Mrotfrac_25[h],Mrotfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hsestream_coll.plot(10**lRbins, Mstreamfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hsestream_coll.fill_between(10**lRbins,Mstreamfrac_25[h],Mstreamfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hseacc_coll.plot(10**lRbins, Maccfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hseacc_coll.fill_between(10**lRbins,Maccfrac_25[h],Maccfrac_75[h],facecolor=colorbins[h],alpha=0.2)
                ax_hsesum_coll.plot(10**lRbins, Msumfrac_50[h], color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h,label=labelwrite)
                if(pversion==0): ax_hsesum_coll.fill_between(10**lRbins,Msumfrac_25[h],Msumfrac_75[h],facecolor=colorbins[h],alpha=0.2)
            
            if(pversion):
                M200_50 = h 

            sE_tot = ax_E_coll.scatter(M200_50,np.log10(EthermR200_50+EkinR200_50),s=120,facecolors='black', edgecolors='black')
            sE_therm = ax_E_coll.scatter(M200_50,np.log10(EthermR200_50),color='black',s=120,facecolors='black', marker='^')
            sE_kin = ax_E_coll.scatter(M200_50,np.log10(EkinR200_50),color='black',s=120,facecolors='black', marker='s')
            sE_feed = ax_E_coll.scatter(M200_50,np.log10(EstarsR200_50+EBHR200_50),s=100,marker='D',facecolors='None', edgecolors='black')
            sE_stars = ax_E_coll.scatter(M200_50,np.log10(EstarsR200_50),facecolors='None',marker='*',s=120, edgecolors='black')
            sE_BH = ax_E_coll.scatter(M200_50,np.log10(EBHR200_50),edgecolors='black',s=40,marker='o', facecolors='None')

            ax_E_coll.scatter(M200_50,np.log10(EthermR200_50+EkinR200_50),s=120,facecolors=colorbins[h], edgecolors=colorbins[h])
            ax_E_coll.scatter(M200_50,np.log10(EthermR200_50),color=colorbins[h],s=120,facecolors=colorbins[h], marker='^')
            ax_E_coll.scatter(M200_50,np.log10(EkinR200_50),color=colorbins[h],s=120,facecolors=colorbins[h], marker='s')
            ax_E_coll.scatter(M200_50,np.log10(EstarsR200_50+EBHR200_50),s=100,marker='D',facecolors='None', edgecolors=colorbins[h])
            ax_E_coll.scatter(M200_50,np.log10(EstarsR200_50),facecolors='None',marker='*',s=120, edgecolors=colorbins[h])
            ax_E_coll.scatter(M200_50,np.log10(EBHR200_50),edgecolors=colorbins[h],s=40,marker='o', facecolors='None')

            if(h==2):
                slambda_DM =ax_lambda_coll.scatter(M200_50,lambdaDM_50,s=120,facecolors='black', edgecolors='black')
                slambda_stars30 =ax_lambda_coll.scatter(M200_50,lambdastars30_50,facecolors='black',marker='*',s=120)
                slambda_stars =ax_lambda_coll.scatter(M200_50,lambdastars_50,facecolors='None',marker='*',s=120, edgecolors='black')#,alpha=float(1-starcenfracR200_50))
                slambda_gas =ax_lambda_coll.scatter(M200_50,lambdagas_50,color='black',s=120,facecolors='None', edgecolors='black')
                slambda_hot =ax_lambda_coll.scatter(M200_50,lambdahot_50,edgecolors='black',s=100,marker='D', facecolors='None',alpha=float(hotfracR200_50))
                slambda_cool =ax_lambda_coll.scatter(M200_50,lambdacool_50,edgecolors='black',s=120,marker='s', facecolors='None',alpha=float(np.sqrt(coolfracR200_50)))
            else:
                ax_lambda_coll.scatter(M200_50,lambdaDM_50,s=120,facecolors='black', edgecolors='black')
                ax_lambda_coll.scatter(M200_50,lambdastars30_50,facecolors='black',marker='*',s=120)
                ax_lambda_coll.scatter(M200_50,lambdastars_50,facecolors='None',marker='*',s=120, edgecolors='black')#,alpha=float(1-starcenfracR200_50))
                ax_lambda_coll.scatter(M200_50,lambdagas_50,color='black',s=120,facecolors='None', edgecolors='black')
                ax_lambda_coll.scatter(M200_50,lambdahot_50,edgecolors='black',s=100,marker='D', facecolors='None',alpha=float(hotfracR200_50))
                ax_lambda_coll.scatter(M200_50,lambdacool_50,edgecolors='black',s=120,marker='s', facecolors='None',alpha=float(np.sqrt(coolfracR200_50)))

            ax_lambda_coll.scatter(M200_50,lambdaDM_50,color=colorbins[h],s=120,facecolors=colorbins[h], edgecolors=colorbins[h])
            ax_lambda_coll.scatter(M200_50,lambdastars30_50,facecolors=colorbins[h],marker='*',s=120)
            ax_lambda_coll.scatter(M200_50,lambdastars_50,facecolors='None',marker='*',s=120, edgecolors=colorbins[h])#,alpha=float(1-starcenfracR200_50))
            ax_lambda_coll.scatter(M200_50,lambdagas_50,color=colorbins[h],s=120,facecolors='None', edgecolors=colorbins[h])
            ax_lambda_coll.scatter(M200_50,lambdahot_50,edgecolors=colorbins[h],s=100,marker='D', facecolors='None',alpha=float(hotfracR200_50))
            ax_lambda_coll.scatter(M200_50,lambdacool_50,edgecolors=colorbins[h],s=120,marker='s', facecolors='None', alpha=float(np.sqrt(coolfracR200_50)))
            if(pversion): ax_lambda_coll.text(h,-0.01,labelwrite,fontsize=14,verticalalignment='center', horizontalalignment='center')

            ax_Jspecific_coll.plot(10**lRbins, np.log10(Jstars_50[h]), color=colorbins[h], lw=2,alpha=1.0,ls='--',zorder=20-h)
            #ax_Jspecific_coll.fill_between(10**lRbins,Jstars_25[h],Jstars_75[h],facecolor=colorbins[h],alpha=0.2)
            ax_Jspecific_coll.plot(10**lRbins, np.log10(Jhot_50[h]), color=colorbins[h], lw=4,alpha=1.0,ls='-',zorder=20-h)
            if(pversion==0): ax_Jspecific_coll.fill_between(10**lRbins,np.log10(Jhot_25[h]),np.log10(Jhot_75[h]),facecolor=colorbins[h],alpha=0.2)
            ax_Jspecific_coll.plot(10**lRbins, np.log10(Jcool_50[h]), color=colorbins[h], lw=4,alpha=1.0,ls=':',zorder=20-h)
            #ax_Jspecific_coll.fill_between(10**lRbins,Jcool_25[h],Jcool_75[h],facecolor=colorbins[h],alpha=0.2)

        #out = scist.binned_statistic(tmp, tmp2,statistic='median', bins=nRbins, range=[(rbinlow,rbinhi)])[0]
        #print h,  lMtot_50[h][60], len(lRfrac_hse_coll[h]) #-1)/71


#out #lMttotfrac_coll[h],Rbins
if(do_ioncol):
    fion.close()

f.close()

Rlow = -1.5
if(plotallhalobins == 1):
    Rlow = -1.2

if(do_hse):
    ax_hse_coll.set_xscale('log')
    ax_hse_coll.set_xlim(10**Rlow,10**0.33)
    ax_hse_coll.set_ylim(-0.5,2.0)
    ax_hse_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hse_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hse_coll.set_yticks([-0.5,0.0,0.5,1.0,1.5,2.0])
    ax_hse_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hse_coll.set_ylabel('$M_{term}/M_{\mathrm{tot}}$',fontsize=16)
    ax_hse_coll.set_ylabel('$\mathcal{S}_{term}$',fontsize=16)
    ax_hse_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hse_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hse_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
    ax_hse_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hse_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hse_coll.savefig('Hse_coll.' + lsfile + '.png')
    
    ax_hsesum_coll.set_xscale('log')
    ax_hsesum_coll.set_xlim(10**Rlow,10**0.33)
    ax_hsesum_coll.set_ylim(-0.5,2.0)
    ax_hsesum_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hsesum_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hsesum_coll.set_yticks([-0.5,0.0,0.5,1.0,1.5,2.0])
    ax_hsesum_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hsesum_coll.set_ylabel('$M_{\mathrm{sum}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hsesum_coll.set_ylabel('$\mathcal{S}_{\mathrm{sum}}$',fontsize=16)
    ax_hsesum_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hsesum_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hsesum_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
    ax_hsesum_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hsesum_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hsesum_coll.savefig('Hsesum_coll.' + lsfile + '.png')
    
    ax_hsetherm_coll.set_xscale('log')
    ax_hsetherm_coll.set_xlim(10**Rlow,10**0.33)
    ax_hsetherm_coll.set_ylim(-0.5,2.0)
    ax_hsetherm_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hsetherm_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hsetherm_coll.set_yticks([-0.5,0.0,0.5,1.0,1.5,2.0])
    ax_hsetherm_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hsetherm_coll.set_ylabel('$M_{\mathrm{therm}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hsetherm_coll.set_ylabel('$\mathcal{S}_{\mathrm{therm}}$',fontsize=16)
    ax_hsetherm_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hsetherm_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hsetherm_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
    ax_hsetherm_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hsetherm_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hsetherm_coll.savefig('Hsetherm_coll.' + lsfile + '.png')
    
    ax_hsesum_coll.set_xscale('log')
    ax_hsesum_coll.set_xlim(10**Rlow,10**0.33)
    ax_hsesum_coll.set_ylim(0.0,2.0)
    ax_hsesum_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hsesum_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hsesum_coll.set_yticks([-0.5,0.0,0.5,1.0,1.5,2.0])
    ax_hsesum_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hsesum_coll.set_ylabel('$M_{\mathrm{sum}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hsesum_coll.set_ylabel('$\mathcal{S}_{\mathrm{sum}}$',fontsize=16)
    ax_hsesum_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hsesum_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hsesum_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
    ax_hsesum_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hsesum_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hsesum_coll.savefig('Hsesum_coll.' + lsfile + '.png')
    
    ax_hserot_coll.set_xscale('log')
    ax_hserot_coll.set_xlim(10**Rlow,10**0.33)
    ax_hserot_coll.set_ylim(-0.1,1.5)
    ax_hserot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hserot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hserot_coll.set_yticks([0.0,0.5,1.0,1.5])
    ax_hserot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hserot_coll.set_ylabel('$M_{\mathrm{rot}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hserot_coll.set_ylabel('$\mathcal{S}_{\mathrm{rot}}$',fontsize=16)
    ax_hserot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hserot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hserot_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hserot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hserot_coll.savefig('Hserot_coll.' + lsfile + '.png')
    
    ax_hsestream_coll.set_xscale('log')
    ax_hsestream_coll.set_xlim(10**Rlow,10**0.33)
    ax_hsestream_coll.set_ylim(-0.2,0.2)
    ax_hsestream_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hsestream_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hsestream_coll.set_yticks([-0.2,-0.1,0.0,0.1,0.2])
    ax_hsestream_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hsestream_coll.set_ylabel('$M_{\mathrm{stream}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hsestream_coll.set_ylabel('$\mathcal{S}_{\mathrm{stream}}$',fontsize=16)
    ax_hsestream_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hsestream_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hsestream_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hsestream_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hsestream_coll.savefig('Hsestream_coll.' + lsfile + '.png')
    
    ax_hseacc_coll.set_xscale('log')
    ax_hseacc_coll.set_xlim(10**Rlow,10**0.33)
    ax_hseacc_coll.set_ylim(-0.6,0.6)
    ax_hseacc_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
    ax_hseacc_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_hseacc_coll.set_yticks([-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6])
    ax_hseacc_coll.set_xlabel('$R/R_{200}$',fontsize=16)
    #ax_hseacc_coll.set_ylabel('$M_{\mathrm{acc}}$/$M_{\mathrm{tot}}$',fontsize=16)
    ax_hseacc_coll.set_ylabel('$\mathcal{S}_{\mathrm{accel}}$',fontsize=16)
    ax_hseacc_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
    ax_hseacc_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
    ax_hseacc_coll.legend(loc="upper right", ncol=2,fontsize=10)
    fig_hseacc_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
    fig_hseacc_coll.savefig('Hseacc_coll.' + lsfile + '.png')

ax_gasfrac_coll.set_xscale('log')
ax_gasfrac_coll.set_xlim(10**Rlow,10**0.33)
ax_gasfrac_coll.set_ylim(-3.0,0.5)
ax_gasfrac_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_gasfrac_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_gasfrac_coll.set_yticks([-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5])
ax_gasfrac_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_gasfrac_coll.set_ylabel(r'log $M($<$R)/M_{200} \times (\Omega_M/\Omega_b)$',fontsize=16)
ax_gasfrac_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_gasfrac_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
if(pversion==0):first_legend = ax_gasfrac_coll.legend(loc="upper left", ncol=1,fontsize=10)
else:first_legend = ax_gasfrac_coll.legend(loc="upper left", ncol=2,fontsize=10)
fig_gasfrac_coll.gca().add_artist(first_legend)
ax_gasfrac_coll.legend((lhotfrac[0],lgasfrac[0],lbarfrac[0]),('$M_{\mathrm{hot}}$','$M_{\mathrm{gas}}$','$M_{\mathrm{bar}}$'),loc='lower right', ncol=1, fontsize=14)
fig_gasfrac_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_gasfrac_coll.savefig('Gasfrac_coll.' + lsfile + '.png')

ax_vradhot_coll.set_xscale('log')
ax_vradhot_coll.set_xlim(10**Rlow,10**0.33)
ax_vradhot_coll.set_ylim(-0.4,1.6)
ax_vradhot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_vradhot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_vradhot_coll.set_yticks([0.0,0.5,1.0,1.5])
ax_vradhot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_vradhot_coll.set_ylabel('$v/v_{200}$',fontsize=16)
ax_vradhot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_vradhot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
first_legend = ax_vradhot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_vradhot_coll.gca().add_artist(first_legend)
#second_legend = ax_vradhot_coll.legend((lvradcoolLstar[0],lsigradcoolLstar[0]),(r'$L^*$ $\langle v_{\mathrm{rad,cool}} \rangle$',r'$L^*$ $\langle v_{\mathrm{rad,cool}}^2 \rangle^{1/2}$'),loc="center right", ncol=1,fontsize=10)
if('lvradcoolLstar' in vars()):
    second_legend = ax_vradhot_coll.legend((lvradcoolLstar[0],lsigradcoolLstar[0]),(r'$L^*$ $v_{\mathrm{rad,cool}}$',r'$L^*$ $\sigma_{\mathrm{rad,cool}}$'),loc="center right", ncol=1,fontsize=10)
    fig_vradhot_coll.gca().add_artist(second_legend)
ax_vradhot_coll.legend((lvradhot[0],lsigradhot[0]),(r'$v_{\mathrm{rad,hot}}$',r'$\sigma_{\mathrm{rad,hot}}$'),loc="upper left", ncol=1,fontsize=14)
fig_vradhot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_vradhot_coll.savefig('Vrad_hot_coll.' + lsfile + '.png')

ax_vtanhot_coll.set_xscale('log')
ax_vtanhot_coll.set_xlim(10**Rlow,10**0.33)
ax_vtanhot_coll.set_ylim(0,2.0)
ax_vtanhot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_vtanhot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_vtanhot_coll.set_yticks([0.0,0.5,1.0,1.5,2.0])
ax_vtanhot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_vtanhot_coll.set_ylabel('$v/v_{200}$',fontsize=16)
ax_vtanhot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_vtanhot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
first_legend = ax_vtanhot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_vtanhot_coll.gca().add_artist(first_legend)
#second_legend = ax_vtanhot_coll.legend((lvtancoolLstar[0],lsigtancoolLstar[0]),(r'$L^*$ $\langle v_{\mathrm{tan,cool}} \rangle$',r'$L^*$ $\langle v_{\mathrm{tan,cool}}^2 \rangle^{1/2}$'),loc="center right", ncol=1,fontsize=10)
if('lvtancoolLstar' in vars()):
    second_legend = ax_vtanhot_coll.legend((lvtancoolLstar[0],lsigtancoolLstar[0]),(r'$L^*$ $v_{\mathrm{tan,cool}}$',r'$L^*$ $\sigma_{\mathrm{tan,cool}}$'),loc="center right", ncol=1,fontsize=10)
    fig_vtanhot_coll.gca().add_artist(second_legend)
#ax_vtanhot_coll.legend((lvtanhot[0],lsigtanhot[0]),(r'$\langle v_{\mathrm{tan,hot}} \rangle$',r'$\langle v_{\mathrm{tan,hot}}^2 \rangle^{1/2}$'),loc="upper left", ncol=1,fontsize=14)
ax_vtanhot_coll.legend((lvtanhot[0],lsigtanhot[0]),(r'$v_{\mathrm{tan,hot}}$',r'$\sigma_{\mathrm{tan,hot}}$'),loc="upper left", ncol=1,fontsize=14)
fig_vtanhot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_vtanhot_coll.savefig('Vtan_hot_coll.' + lsfile + '.png')

ax_sighot_coll.set_xscale('log')
ax_sighot_coll.set_xlim(10**Rlow,10**0.33)
ax_sighot_coll.set_ylim(0,2.0)
ax_sighot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_sighot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_sighot_coll.set_yticks([0.0,0.5,1.0,1.5,2.0])
ax_sighot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_sighot_coll.set_ylabel('$v/v_{200}$',fontsize=16)
ax_sighot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_sighot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
first_legend = ax_sighot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_sighot_coll.gca().add_artist(first_legend)
#ax_sighot_coll.legend((lsighot[0],lsigDM[0]),(r'$\langle v_{\mathrm{hot}}^2 \rangle^{1/2}$',r'$\langle v_{\mathrm{DM}}^2 \rangle^{1/2}$'),loc="upper left", ncol=1,fontsize=14)
ax_sighot_coll.legend((lsighot[0],lsigDM[0]),(r'$\sigma_{\mathrm{hot}}$',r'$\sigma_{\mathrm{DM}}$'),loc="upper left", ncol=1,fontsize=14)
fig_sighot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_sighot_coll.savefig('Sig_coll.' + lsfile + '.png')

ax_Jspecific_coll.set_xlim(Rlow,0.3)
ax_Jspecific_coll.set_ylim(11.0,16.0)
ax_Jspecific_coll.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_Jspecific_coll.set_ylabel('$J$ [$M_{\odot}$ kpc km s$^{-1}$]',fontsize=16)
ax_Jspecific_coll.plot([0,0],[-10,10],color='k',lw=1,ls=':')
#ax_Jspecific_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Jspecific_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_Jspecific_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_Jspecific_coll.savefig('Jspecific_coll.' + lsfile + '.png')

ax_nHhot_coll.set_xscale('log')
ax_nHhot_coll.set_xlim(10**Rlow,10**0.33)
ax_nHhot_coll.set_ylim(-1.0,3.0)
ax_nHhot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_nHhot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_nHhot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_nHhot_coll.set_ylabel('log $n_{\mathrm{H}}/n_{\mathrm{H},200}$',fontsize=16)
ax_nHhot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_nHhot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_nHhot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_nHhot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_nHhot_coll.savefig('nH_hot_coll.' + lsfile + '.png')

ax_Thot_coll.set_xscale('log')
ax_Thot_coll.set_xlim(10**Rlow,10**0.33)
ax_Thot_coll.set_ylim(-1.0,1.0)
ax_Thot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Thot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Thot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Thot_coll.set_ylabel('log $T/T_{200}$',fontsize=16)
ax_Thot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Thot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Thot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_Thot_coll.subplots_adjust(left=0.19, bottom=0.17,top=0.98,right=0.98)
fig_Thot_coll.savefig('T_hot_coll.' + lsfile + '.png')

ax_Zhot_coll.set_xscale('log')
ax_Zhot_coll.set_xlim(10**Rlow,10**0.33)
ax_Zhot_coll.set_ylim(-4.0,-1.0)
ax_Zhot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Zhot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Zhot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Zhot_coll.set_ylabel('log $Z$',fontsize=16)
ax_Zhot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Zhot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Zhot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_Zhot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_Zhot_coll.savefig('Z_hot_coll.' + lsfile + '.png')

ax_Phot_coll.set_xscale('log')
ax_Phot_coll.set_xlim(10**Rlow,10**0.33)
ax_Phot_coll.set_ylim(-3.0,3.0)
ax_Phot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Phot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Phot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Phot_coll.set_ylabel('log $P/P_{200}$',fontsize=16)
ax_Phot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Phot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Phot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_Phot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_Phot_coll.savefig('P_hot_coll.' + lsfile + '.png')

ax_Shot_coll.set_xscale('log')
ax_Shot_coll.set_xlim(10**Rlow,10**0.33)
ax_Shot_coll.set_ylim(-3.0,3.0)
ax_Shot_coll.set_xticks([0.03,0.1,0.3,1.0,2.0])
ax_Shot_coll.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
ax_Shot_coll.set_xlabel('$R/R_{200}$',fontsize=16)
ax_Shot_coll.set_ylabel('log $S/S_{200}$',fontsize=16)
ax_Shot_coll.plot([1,1],[-10,10],color='k',lw=1,ls=':')
ax_Shot_coll.plot([-10,10],[0,0],color='k',lw=1,ls=':')
ax_Shot_coll.legend(loc="upper right", ncol=1,fontsize=10)
fig_Shot_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.98,right=0.98)
fig_Shot_coll.savefig('S_hot_coll.' + lsfile + '.png')

if (do_hse):
    ax_hse_all.set_xlim(Rlow,0.3)
    ax_hse_all.set_ylim(-4.0,0.5)
    ax_hse_all.set_xlabel('log $R/R_{200}$',fontsize=16)
    ax_hse_all.set_ylabel('log $M/M_{200}$',fontsize=16)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
    sm._A = []
    cbar = fig_hse_all.colorbar(sm)
    cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
    cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
    fig_hse_all.subplots_adjust(left=0.17, bottom=0.17)
    fig_hse_all.savefig('Hse_all.' + lsfile + '.png')
    
ax_velDM_all.set_xlim(Rlow,0.3)
ax_velDM_all.set_ylim(0.0,2.0)
ax_velDM_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_velDM_all.set_ylabel('log $\sigma_{\mathrm{tot}}/v_{200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_velDM_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_velDM_all.subplots_adjust(left=0.17, bottom=0.17)
fig_velDM_all.savefig('SigDM_all.' + lsfile + '.png')

ax_vradhot_all.set_xlim(Rlow,0.3)
ax_vradhot_all.set_ylim(-1.0,2.0)
ax_vradhot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_vradhot_all.set_ylabel('$v_{\mathrm{rad}}/v_{200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_vradhot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_vradhot_all.subplots_adjust(left=0.17, bottom=0.17)
fig_vradhot_all.savefig('Vradhot_all.' + lsfile + '.png')

ax_vtanhot_all.set_xlim(Rlow,0.3)
ax_vtanhot_all.set_ylim(-1.0,2.0)
ax_vtanhot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_vtanhot_all.set_ylabel('$v_{\mathrm{tan}}/v_{200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_vtanhot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_vtanhot_all.subplots_adjust(left=0.17, bottom=0.17)
fig_vtanhot_all.savefig('Vtanhot_all.' + lsfile + '.png')

ax_sighot_all.set_xlim(Rlow,0.3)
ax_sighot_all.set_ylim(-1.0,2.0)
ax_sighot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_sighot_all.set_ylabel('$\sigma/v_{200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_sighot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_sighot_all.subplots_adjust(left=0.17, bottom=0.17)
fig_sighot_all.savefig('Sig_all.' + lsfile + '.png')

ax_gasfrac_all.set_xlim(Rlow,0.3)
ax_gasfrac_all.set_ylim(-3.0,0.5)
ax_gasfrac_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_gasfrac_all.set_ylabel('log $f_{gas}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_gasfrac_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_gasfrac_all.subplots_adjust(left=0.17, bottom=0.17)
fig_gasfrac_all.savefig('Gasfrac_all.' + lsfile + '.png')

ax_DMfrac_all.set_xlim(Rlow,0.3)
ax_DMfrac_all.set_ylim(-3.0,0.5)
ax_DMfrac_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_DMfrac_all.set_ylabel('log $f_{DM}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_DMfrac_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_DMfrac_all.subplots_adjust(left=0.17, bottom=0.17)
fig_DMfrac_all.savefig('DMfrac_all.' + lsfile + '.png')

ax_Jtot_all.set_xlim(Rlow,0.3)
ax_Jtot_all.set_ylim(10, 15)
ax_Jtot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_Jtot_all.set_ylabel('log $J_{\mathrm{tot}}$ [kpc km s$^{-1}$ M$_\odot$]',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_Jtot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_Jtot_all.subplots_adjust(left=0.17, bottom=0.17)
fig_Jtot_all.savefig('Jtot_all.' + lsfile + '.png')

ax_Thot_all.set_xlim(Rlow,0.3)
ax_Thot_all.set_ylim(-1.0, 1.0) #5.0,8.0)
ax_Thot_all.set_yticks([-1.0,-0.5,0.0,0.5,1.0])
ax_Thot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_Thot_all.set_ylabel('log $T/T_{200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_Thot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_Thot_all.subplots_adjust(left=0.19, bottom=0.17)
fig_Thot_all.savefig('Thot_all.' + lsfile + '.png')


ax_nHhot_all.set_xlim(Rlow,0.3)
ax_nHhot_all.set_ylim(-1.0, 3.0) #5.0,8.0)
ax_nHhot_all.set_xlabel('log $R/R_{200}$',fontsize=16)
ax_nHhot_all.set_ylabel('$n_{\mathrm{H}}/n_{\mathrm{H},200}$',fontsize=16)
sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=mhlow, vmax=mhhi))
sm._A = []
cbar = fig_nHhot_all.colorbar(sm)
cbar.set_label('log $M_{\mathrm{200}}$ [$M_{\odot}$]', fontsize=10)
cbar.set_ticks([mhlow,(3*mhlow+mhhi)/4.,(mhlow+mhhi)/2.,(mhlow+3*mhhi)/4.,mhhi])
fig_nHhot_all.subplots_adjust(left=0.17, bottom=0.17)
fig_nHhot_all.savefig('nHhot_all.' + lsfile + '.png')

if(plotallhalobins):
    ax_E_coll.set_xlim(10.9,16.5)
    ax_E_coll.set_ylim(-1.3,2.5)
else:
    ax_E_coll.set_xlim(11.8,15.99)
    ax_E_coll.set_ylim(-1.3,1.5)
if(pversion):
    ax_E_coll.set_xlim(1.5,10.5)
    ax_E_coll.xaxis.set_major_formatter(plt.NullFormatter())
    ax_E_coll.xaxis.set_ticks_position('none') 
ax_E_coll.set_xlabel('log $M_{200}$ [M$_\odot$]',fontsize=16)
ax_E_coll.set_ylabel('log $E/E_{\mathrm{bind}}$',fontsize=16)
ax_E_coll.plot([-10,20],[-0.301,-0.301],color='k',lw=1,ls=':')
#ax_E_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
#ax_E_coll.legend(loc="upper right", ncol=2,fontsize=10)
ax_E_coll.legend([sE_therm,sE_kin,sE_tot,sE_stars,sE_BH,sE_feed],['$E_{\mathrm{therm}}$','$E_{\mathrm{kin}}$','$E_{\mathrm{halo}}$','$E_{*}$','$E_{\mathrm{BH}}$','$E_{\mathrm{feed}}$'],loc='upper right',ncol=1,fontsize=14)
fig_E_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.97,right=0.97)
fig_E_coll.savefig('E_coll.' + lsfile + '.png')

if(plotallhalobins):
    ax_lambda_coll.set_xlim(10.9,16.5)
    ax_lambda_coll.set_ylim(0,0.16)
else:
    ax_lambda_coll.set_xlim(11.8,15.99)
    ax_lambda_coll.set_ylim(0,0.15)
if(pversion):
    ax_lambda_coll.set_xlim(1.5,10.5)
    ax_lambda_coll.xaxis.set_major_formatter(plt.NullFormatter())
    ax_lambda_coll.xaxis.set_ticks_position('none') 
if(pversion==0):
    ax_lambda_coll.set_xlabel('log $M_{200}$ [M$_\odot$]',fontsize=16)
ax_lambda_coll.set_ylabel('$\lambda$',fontsize=16)
ax_lambda_coll.plot([-10,20],[-0.301,-0.301],color='k',lw=1,ls=':')
#ax_lambda_coll.plot([-10,10],[1,1],color='k',lw=1,ls=':')
ax_lambda_coll.legend([slambda_DM,slambda_stars,slambda_stars30,slambda_gas,slambda_hot,slambda_cool],['$\lambda_{\mathrm{DM}}$','$\lambda_{\mathrm{*}}$','$\lambda_{\mathrm{*,c}}$','$\lambda_{\mathrm{gas}}$','$\lambda_{\mathrm{hot}}$','$\lambda_{\mathrm{cool}}$'],loc='upper right', ncol=1, fontsize=14)
fig_lambda_coll.subplots_adjust(left=0.17, bottom=0.17,top=0.97,right=0.97)
fig_lambda_coll.savefig('lambda_coll.' + lsfile + '.png')

