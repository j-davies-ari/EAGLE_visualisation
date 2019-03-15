import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

coldensmap_file = sys.argv[1]
rlow = 0.0
rhi = float(sys.argv[2])
colorbar = sys.argv[3]

nh, no, nh1, nc2, nc3, nc4, no1, no6, nmg2, nsi2, nsi3, r = np.loadtxt(coldensmap_file,usecols=(2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)

rbins = np.linspace(rhi*0.05,rhi*0.95,10)
nbins = 10

indexes = np.where((r>=rlow) & (r<rhi))
indexes150_170 = np.where((nh1>=15.0) & (nh1<17.0) & (r>=rlow) & (r<rhi))
indexes160_185 = np.where((nh1>=16.0) & (nh1<18.5) & (r>=rlow) & (r<rhi))
indexes150 = np.where((nh1>=15.0) & (r>=rlow) & (r<rhi))
indexes160 = np.where((nh1>=16.0) & (r>=rlow) & (r<rhi))
indexes170 = np.where((nh1>=17.0) & (r>=rlow) & (r<rhi))
indexesDLA = np.where((nh1>=20.3) & (r>=rlow) & (r<rhi))

covfrac150 = len(nh1[indexes150])*1./len(nh1[indexes])
covfrac160 = len(nh1[indexes160])*1./len(nh1[indexes])

n_bins = np.histogram(r, nbins, range=[rlow,rhi])[0]
n_bins150 = np.histogram(r[indexes150], nbins, range=[rlow,rhi])[0]
n_bins160 = np.histogram(r[indexes160], nbins, range=[rlow,rhi])[0]
n_bins170 = np.histogram(r[indexes170], nbins, range=[rlow,rhi])[0]
n_binsDLA = np.histogram(r[indexesDLA], nbins, range=[rlow,rhi])[0]
nh1_bins = np.histogram(r, nbins, range=[rlow,rhi], weights=nh1)[0]

for i in range(len(rbins)):
    print rbins[i], n_bins150[i]*1./n_bins[i], n_bins160[i]*1./n_bins[i], n_bins170[i]*1./n_bins[i], n_binsDLA[i]*1./n_bins[i], nh1_bins[i]/n_bins[i]

nZbins=30
no1_h1 = no1 - nh1
no1_h1_low = -8.0
no1_h1_hi = -2.0
no_h = no - nh
no_h_low = -8.0
no_h_hi = -2.0
nmg2_h1 = nmg2 - nh1
nmg2_h1_low = -8.0
nmg2_h1_hi = -2.0

Zbinplot = np.zeros(nZbins)
no_h_160_185, Zbins =  np.histogram(no_h[indexes160_185], nZbins, range=[no_h_low, no_h_hi])
no1_h1_160_185 = np.histogram(no1_h1[indexes160_185], nZbins, range=[no1_h1_low, no1_h1_hi])[0]
nmg2_h1_160_185 = np.histogram(nmg2_h1[indexes160_185], nZbins, range=[nmg2_h1_low, nmg2_h1_hi])[0]
indexes160_185_mg2_125 = np.where((nh1>=16.0) & (nh1<18.5) & (nmg2>12.5) & (r>=rlow) & (r<rhi))
nmg2_h1_160_185_mg2_125 = np.histogram(nmg2_h1[indexes160_185_mg2_125], nZbins, range=[nmg2_h1_low, nmg2_h1_hi])[0]

for i in range(nZbins):
    Zbinplot[i] = (Zbins[i]+Zbins[i+1])/2.
    print (Zbins[i]+Zbins[i+1])/2., no_h_160_185[i]*1./np.sum(no_h_160_185), no1_h1_160_185[i]*1./np.sum(no1_h1_160_185), nmg2_h1_160_185[i]*1./np.sum(nmg2_h1_160_185), nmg2_h1_160_185_mg2_125[i]*1./np.sum(nmg2_h1_160_185_mg2_125)

f_160_185_no_nh = plt.figure(figsize=(4.0,3.5))
ax160_185_no_nh = f_160_185_no_nh.add_subplot(111)

ax160_185_no_nh.hist(Zbinplot,30,weights=no_h_160_185*1./np.sum(no_h_160_185),label="",color=colorbar,alpha=0.5)

ax160_185_no_nh.set_xlim(-6.0,-2.0)
#ax160_185_no_nh.set_ylim(0,1e+06)
ax160_185_no_nh.set_xlabel('log $Z$',fontsize=10)
ax160_185_no_nh.set_ylabel('$f$',fontsize=10)
ax160_185_no_nh.legend(loc='upper left',fontsize=10)
ax160_185_no_nh.yaxis.set_ticklabels([])
f_160_185_no_nh.savefig('HI_160_185_no_nh_hists.' + coldensmap_file + '.png')


print covfrac150
