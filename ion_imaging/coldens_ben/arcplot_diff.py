import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter
from matplotlib import gridspec



sim= '.'
diff_file_base = sys.argv[1]
direction = sys.argv[2]
length = sys.argv[3]
theta_lo = float(sys.argv[4])
theta_hi = float(sys.argv[5])

npix = 400

r_scale = float(length)*1000./npix

r = np.arange(200)*r_scale

#colion_map_diff.snap_noneq_034_z000p149_shalo_1_12.24_9.81_1.260.o7.z.l1.0.dat

cen_x = npix/2
cen_y = npix/2

rlow = 0.0
rhi = 600  

fh1_name = '%s.h1.%s.l%s.dat'%(diff_file_base,direction,length)
fc4_name = '%s.c4.%s.l%s.dat'%(diff_file_base,direction,length)
fn5_name = '%s.n5.%s.l%s.dat'%(diff_file_base,direction,length)
fo6_name = '%s.o6.%s.l%s.dat'%(diff_file_base,direction,length)
fo7_name = '%s.o7.%s.l%s.dat'%(diff_file_base,direction,length)
fo8_name = '%s.o8.%s.l%s.dat'%(diff_file_base,direction,length)
fo9_name = '%s.o9.%s.l%s.dat'%(diff_file_base,direction,length)
fne8_name = '%s.ne8.%s.l%s.dat'%(diff_file_base,direction,length)
fmg10_name = '%s.mg10.%s.l%s.dat'%(diff_file_base,direction,length)
fsi12_name = '%s.si12.%s.l%s.dat'%(diff_file_base,direction,length)
ffe17_name = '%s.fe17.%s.l%s.dat'%(diff_file_base,direction,length)

print fh1_name

f = plt.figure(figsize = (6,8))

gs = gridspec.GridSpec(5,1,width_ratios=[1], height_ratios=[1,1,1,2,2], hspace=0)
axnh = plt.subplot(gs[0])
axT = plt.subplot(gs[1])
axv = plt.subplot(gs[2])
axcolne = plt.subplot(gs[3])
axcoldiff = plt.subplot(gs[4])

#axnh = f.add_subplot(511)
#axT = f.add_subplot(512)
#axv = f.add_subplot(513)
#axcolne = f.add_subplot(514)
#axcoldiff = f.add_subplot(515)


x, y, h1_ne, h1_eq, colh, nh, T, v = np.loadtxt(fc4_name,usecols=(0,1,2,3,4,5,6,9),unpack=True)
c4_ne, c4_eq = np.loadtxt(fc4_name,usecols=(2,3),unpack=True)
n5_ne, n5_eq = np.loadtxt(fn5_name,usecols=(2,3),unpack=True)
o6_ne, o6_eq = np.loadtxt(fo6_name,usecols=(2,3),unpack=True)
o7_ne, o7_eq = np.loadtxt(fo7_name,usecols=(2,3),unpack=True)
o8_ne, o8_eq = np.loadtxt(fo8_name,usecols=(2,3),unpack=True)
o9_ne, o9_eq = np.loadtxt(fo9_name,usecols=(2,3),unpack=True)
ne8_ne, ne8_eq = np.loadtxt(fne8_name,usecols=(2,3),unpack=True)
mg10_ne, mg10_eq = np.loadtxt(fmg10_name,usecols=(2,3),unpack=True)
si12_ne, si12_eq = np.loadtxt(fsi12_name,usecols=(2,3),unpack=True)
fe17_ne, fe17_eq = np.loadtxt(ffe17_name,usecols=(2,3),unpack=True)

dx = x-cen_x
dy = y-cen_y
dist = np.sqrt(dx**2+dy**2)

if((theta_lo>90) & (theta_lo< 270)):
    theta_lo -= 180.
    theta_hi -= 180.
    dx = -dx
    dy = -dy

indexes = np.where((np.arccos(dx/(dist))*180/np.pi<theta_hi) & (np.arccos(dx/(dist))*180/np.pi>theta_lo) & (np.arcsin(dy/(dist))*180/np.pi<theta_hi) & (np.arcsin(dy/(dist))*180/np.pi>theta_lo))

print dx
print dy
print dist[indexes]
#print indexes
print x[indexes]
print y[indexes]
print len(x), len(x[indexes])

n = np.histogram(dist[indexes],200,range=[0,200])[0]
T_ave = np.histogram(dist[indexes],200,range=[0,200],weights=T[indexes])[0]
T_ave = T_ave/n

nh_ave = np.histogram(dist[indexes],200,range=[0,200],weights=nh[indexes])[0]
nh_ave = nh_ave/n

v_ave = np.histogram(dist[indexes],200,range=[0,200],weights=v[indexes])[0]
v_ave = v_ave/n

h1_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=h1_ne[indexes])[0]
h1_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=h1_eq[indexes])[0]
h1_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=h1_ne[indexes]-h1_eq[indexes])[0]
h1_ne_ave = h1_ne_ave/n
h1_eq_ave = h1_eq_ave/n
h1_diff_ave = h1_diff_ave/n

c4_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=c4_ne[indexes])[0]
c4_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=c4_eq[indexes])[0]
c4_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=c4_ne[indexes]-c4_eq[indexes])[0]
c4_ne_ave = c4_ne_ave/n
c4_eq_ave = c4_eq_ave/n
c4_diff_ave = c4_diff_ave/n

n5_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=n5_ne[indexes])[0]
n5_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=n5_eq[indexes])[0]
n5_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=n5_ne[indexes]-n5_eq[indexes])[0]
n5_ne_ave = n5_ne_ave/n
n5_eq_ave = n5_eq_ave/n
n5_diff_ave = n5_diff_ave/n

o6_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o6_ne[indexes])[0]
o6_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o6_eq[indexes])[0]
o6_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o6_ne[indexes]-o6_eq[indexes])[0]
o6_ne_ave = o6_ne_ave/n
o6_eq_ave = o6_eq_ave/n
o6_diff_ave = o6_diff_ave/n

o7_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o7_ne[indexes])[0]
o7_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o7_eq[indexes])[0]
o7_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o7_ne[indexes]-o7_eq[indexes])[0]
o7_ne_ave = o7_ne_ave/n
o7_eq_ave = o7_eq_ave/n
o7_diff_ave = o7_diff_ave/n

o8_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o8_ne[indexes])[0]
o8_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o8_eq[indexes])[0]
o8_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o8_ne[indexes]-o8_eq[indexes])[0]
o8_ne_ave = o8_ne_ave/n
o8_eq_ave = o8_eq_ave/n
o8_diff_ave = o8_diff_ave/n

o9_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o9_ne[indexes])[0]
o9_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o9_eq[indexes])[0]
o9_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=o9_ne[indexes]-o9_eq[indexes])[0]
o9_ne_ave = o9_ne_ave/n
o9_eq_ave = o9_eq_ave/n
o9_diff_ave = o9_diff_ave/n

ne8_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=ne8_ne[indexes])[0]
ne8_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=ne8_eq[indexes])[0]
ne8_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=ne8_ne[indexes]-ne8_eq[indexes])[0]
ne8_ne_ave = ne8_ne_ave/n
ne8_eq_ave = ne8_eq_ave/n
ne8_diff_ave = ne8_diff_ave/n

mg10_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=mg10_ne[indexes])[0]
mg10_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=mg10_eq[indexes])[0]
mg10_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=mg10_ne[indexes]-mg10_eq[indexes])[0]
mg10_ne_ave = mg10_ne_ave/n
mg10_eq_ave = mg10_eq_ave/n
mg10_diff_ave = mg10_diff_ave/n

si12_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=si12_ne[indexes])[0]
si12_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=si12_eq[indexes])[0]
si12_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=si12_ne[indexes]-si12_eq[indexes])[0]
si12_ne_ave = si12_ne_ave/n
si12_eq_ave = si12_eq_ave/n
si12_diff_ave = si12_diff_ave/n

fe17_ne_ave = np.histogram(dist[indexes],200,range=[0,200],weights=fe17_ne[indexes])[0]
fe17_eq_ave = np.histogram(dist[indexes],200,range=[0,200],weights=fe17_eq[indexes])[0]
fe17_diff_ave = np.histogram(dist[indexes],200,range=[0,200],weights=fe17_ne[indexes]-fe17_eq[indexes])[0]
fe17_ne_ave = fe17_ne_ave/n
fe17_eq_ave = fe17_eq_ave/n
fe17_diff_ave = fe17_diff_ave/n

r 

print T_ave
print o8_ne_ave
print o8_eq_ave
print o8_diff_ave

axT.plot(r, T_ave, color="red", lw=2)
axnh.plot(r, nh_ave, color="blue", lw=2)
axv.plot(r, v_ave, color="lime", lw=2)

#axcolne.plot(r,h1_ne_ave,color="black", lw=1)
#axcolne.plot(r,n5_ne_ave,color="blue", lw=1)
axcolne.plot(r,o6_ne_ave,color="blue", lw=2)
axcolne.plot(r,o7_ne_ave,color="cyan", lw=2)
axcolne.plot(r,o8_ne_ave,color="lime", lw=2)
axcolne.plot(r,o9_ne_ave,color="gold", lw=2)
axcolne.plot(r,c4_ne_ave,color="purple", lw=2)
axcolne.plot(r,ne8_ne_ave,color="darkorange", lw=2)
#axcolne.plot(r,mg10_ne_ave,color="darkorange", lw=2)
axcolne.plot(r,si12_ne_ave,color="red", lw=2)
axcolne.plot(r,fe17_ne_ave,color="magenta", lw=2)

#axcolne.plot(r,h1_eq_ave,color="black", lw=1)
#axcolne.plot(r,n5_eq_ave,color="blue", lw=1)
axcolne.plot(r,o6_eq_ave,color="blue", lw=2, linestyle=":")
axcolne.plot(r,o7_eq_ave,color="cyan", lw=2, linestyle=":")
axcolne.plot(r,o8_eq_ave,color="lime", lw=2, linestyle=":")
axcolne.plot(r,o9_eq_ave,color="gold", lw=2, linestyle=":")
axcolne.plot(r,c4_eq_ave,color="purple", lw=2, linestyle=":")
axcolne.plot(r,ne8_eq_ave,color="darkorange", lw=2, linestyle=":")
#axcolne.plot(r,mg10_eq_ave,color="darkorange", lw=2, linestyle=":")
axcolne.plot(r,si12_eq_ave,color="red", lw=2, linestyle=":")
axcolne.plot(r,fe17_eq_ave,color="magenta", lw=2, linestyle=":")

#axcoldiff.plot(r,h1_diff_ave,color="black", lw=2, label='HI')
#axcoldiff.plot(r,n5_diff_ave,color="blue", lw=2, label='NV')
axcoldiff.plot(r,o6_diff_ave,color="blue", lw=2, label='OVI', zorder=60)
axcoldiff.plot(r,o7_diff_ave,color="cyan", lw=2, label='OVII', zorder=60)
axcoldiff.plot(r,o8_diff_ave,color="lime", lw=2, label='OVIII', zorder=60)
axcoldiff.plot(r,o9_diff_ave,color="gold", lw=2, label='OIX', zorder=60)
axcoldiff.plot(r,c4_diff_ave,color="purple", lw=2, label='CIV')
axcoldiff.plot(r,ne8_diff_ave,color="darkorange", lw=2, label='NeVIII')
#axcoldiff.plot(r,mg10_diff_ave,color="darkorange", lw=2, label='MgX')
axcoldiff.plot(r,si12_diff_ave,color="red", lw=2, label='SiXII')
axcoldiff.plot(r,fe17_diff_ave,color="magenta", lw=2, label='FeXVII')

axnh.set_xlim(rlow,rhi)
axnh.set_ylim(-5.49,-2.0)
#axnh.set_xlabel('$r$ [kpc]', fontsize=14)
axnh.set_ylabel('log $n_{\mathrm{H}}$ [cm$^{-3}$]', fontsize=14)
axnh.xaxis.set_major_formatter( NullFormatter() )

axT.set_xlim(rlow,rhi)
axT.set_ylim(5.0,7.5)
#axT.set_xlabel('$r$ [kpc]', fontsize=14)
axT.set_ylabel('log $T$ [K]', fontsize=14)
axT.xaxis.set_major_formatter( NullFormatter() )

axv.set_xlim(rlow,rhi)
axv.set_ylim(-600,199)
#axv.set_xlabel('$r$ [kpc]', fontsize=14)
axv.set_ylabel('$v$ [km/s]', fontsize=14)
axv.xaxis.set_major_formatter( NullFormatter() )

axcolne.set_xlim(rlow,rhi)
axcolne.set_ylim(11.01,16.5)
#axcolne.set_xlabel('$r$ [kpc]', fontsize=14)
axcolne.set_ylabel('log $N$ [cm$^{-2}$]', fontsize=14)
axcolne.xaxis.set_major_formatter( NullFormatter() )

axcoldiff.set_xlim(rlow,rhi)
axcoldiff.set_ylim(-2.0,1.0)
axcoldiff.set_xlabel('$r$ [kpc]', fontsize=14)
axcoldiff.set_ylabel('$\Delta$log $N$ [cm$^{-2}$]', fontsize=14)
axcoldiff.legend(loc='lower left', fontsize=12, ncol=2)

#f.savefig('nh_T_colne_colde' + sys.argv[1] + '.png')
f.subplots_adjust(left=0.15)
f.savefig('nh_T_v_colne_coldiff.%s.%s.l%s.theta%s_%s.png'%(diff_file_base,direction,length,sys.argv[4],sys.argv[5]))


#print n

#for i in x:
#    dx = x[i]-cen_x
#    dy = y[i]-cen_y
#    if((np.arccos(dx/(dx**2+dy**2))*180/np.pi<theta_hi) & (np.arccos(dx/(dx**2+dy**2))*180/np.pi>theta_lo) & (np.arcsin(dy/(dx**2+dy**2))*180/np.pi<theta_hi) & (np.arcsin(dy/(dx**2+dy**2))*180/np.pi>theta_lo))

