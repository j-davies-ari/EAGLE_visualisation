#import eagle
import sys
#import coldens
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axgrid
from pylab import get_cmap
from matplotlib.ticker import NullFormatter,MultipleLocator, FormatStrFormatter

#def plot_o6_c4_h1_plot

#o6col =  np.loadtxt(o6file, usecols=(2), unpack=True)
#c4col =  np.loadtxt(c4file, usecols=(2), unpack=True)
#h1col =  np.loadtxt(h1file, usecols=(2), unpack=True) 
#ax_.plot


def plot_ion_ratio(xval,yval,xlo,xhi,ylo,yhi,xlab,ylab,plotname,fraction=None):

    f = plt.figure(figsize = (6,6))

    ax = f.add_subplot(111)

    hex1 = ax.hexbin(xval,yval,gridsize=100, C=xval*0.0+1.0,vmin=0,vmax=npixels_radius*0.01,reduce_C_function=np.nansum,zorder=10,cmap=get_cmap("gist_ncar_r"))#get_cmap("YlGnBu"))
    ax.axis([xlo,xhi,ylo,yhi])
    ax.set_xlabel(r'%s'%xlab,fontsize=16)
    ax.set_ylabel(r'%s'%ylab,fontsize=16)
    ax.set_title(r'%s'%title,fontsize=16)

    div1 = axgrid.make_axes_locatable(ax)
    cax1 = div1.append_axes("right", size="5%", pad=-0.2, zorder=20)
    cb1 = f.colorbar(hex1,cax=cax1)
    #cb1.set_label('n$_\mathrm{LOS}$')
    cb1.ax.set_xlabel('n$_\mathrm{LOS}$',labelpad=12)

    if fraction != None:
        ax.text((xhi-xlo)*0.86+xlo,(yhi-ylo)*0.92+ylo,r'f=%5.3f'%(float(fraction)),fontsize=14,horizontalalignment='center', verticalalignment='center',zorder=100)

    ax2 = ax.twinx()
    ax2.yaxis.set_major_formatter( NullFormatter() )
    n,bins,patches = ax2.hist(xval,20, normed=1,range=[xlo,xhi], histtype='step',zorder=20)
    ax2.axis([xlo,xhi,0,2])
    ax2.set_ylim(0,5.5)
    
    ax3 = ax.twiny()
    ax3.xaxis.set_major_formatter( NullFormatter() )
    n,bins,patches = ax3.hist(yval, 20, normed=1,range=[ylo,yhi], histtype='step',orientation='horizontal',zorder=20)
    ax3.axis([0,2,ylo,yhi])
    ax3.set_xlim(0,5.5)

    f.subplots_adjust(left=0.135, bottom=0.135)
    f.savefig(coldens_basename + '.' + plotname + '.' + dir + '.l' + size + '.png') 





coldens_basename = sys.argv[1]
title = sys.argv[2] 

dir = 'z'
size = '0.6'
radius = 200.

o6filename = coldens_basename + '.o6.' + dir + '.l' + size + '.dat' 
c4filename = coldens_basename + '.n5.' + dir + '.l' + size + '.dat' 
h1filename = coldens_basename + '.h1.' + dir + '.l' + size + '.dat' 
n5filename = coldens_basename + '.n5.' + dir + '.l' + size + '.dat' 

print o6filename
print c4filename
 
h1lo = 13.0
h1hi = 18.0
o6lo = 12.0
o6hi = 16.0
n5lo = 11.0
n5hi = 15.0
c4lo = 12.0
c4hi = 16.0
o6c4lo = -1.5
o6c4hi = 1.5
o6n5lo = -1.0
o6n5hi = 2.0
o6h1lo = -3.5
o6h1hi = 1.5
n5h1lo = -4.5
n5h1hi = 0.5
c4h1lo = -3.5
c4h1hi = 1.5
n5po6lo = 26.0
n5po6hi = 29.0

o6collim = 13.5
n5collim = 12.5
c4collim = 13.0

o6col, dist =  np.loadtxt(o6filename, usecols=(2,3), unpack=True)
c4col, dist =  np.loadtxt(c4filename, usecols=(2,3), unpack=True)
h1col, dist =  np.loadtxt(h1filename, usecols=(2,3), unpack=True)
n5col, dist =  np.loadtxt(n5filename, usecols=(2,3), unpack=True)

indexes_radius = np.where(dist<radius)
npixels_radius = len(h1col[indexes_radius])

indexes_c4_o6 = np.where((o6col>o6collim) & (c4col>c4collim) & (dist<radius))
indexes_n5_o6 = np.where((o6col>o6collim) & (n5col>n5collim) & (dist<radius))
indexes_h1_o6 = np.where((dist<radius))
indexes_h1_n5 = np.where((dist<radius))
indexes_h1_c4 = np.where((dist<radius))
indexes_h1 = np.where((dist<radius))


print "length_radius= ", npixels_radius
print "length_n5_o6= ", len(o6col[indexes_n5_o6]), len(o6col)
print "length_c4_o6= ", len(o6col[indexes_c4_o6]), len(o6col)

plot_ion_ratio(h1col[indexes_c4_o6],o6col[indexes_c4_o6]-c4col[indexes_c4_o6],h1lo,h1hi,o6c4lo,o6c4hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{CIV}}$]','h1_o6c4',fraction=len(h1col[indexes_c4_o6])*1./npixels_radius)

plot_ion_ratio(h1col[indexes_n5_o6],o6col[indexes_n5_o6]-n5col[indexes_n5_o6],h1lo,h1hi,o6n5lo,o6n5hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{NV}}$]','h1_o6n5',fraction=len(h1col[indexes_n5_o6])*1./npixels_radius)

plot_ion_ratio(h1col[indexes_h1_o6],o6col[indexes_h1_o6]-h1col[indexes_h1_o6],h1lo,h1hi,o6h1lo,o6h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','h1_o6h1',fraction=len(h1col[indexes_h1_o6])*1./npixels_radius)

plot_ion_ratio(h1col[indexes_h1_n5],n5col[indexes_h1_n5]-h1col[indexes_h1_n5],h1lo,h1hi,n5h1lo,n5h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','h1_n5h1',fraction=len(h1col[indexes_h1_n5])*1./npixels_radius)

plot_ion_ratio(h1col[indexes_h1_c4],c4col[indexes_h1_c4]-h1col[indexes_h1_c4],h1lo,h1hi,c4h1lo,c4h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{CIV}}/N_{\mathrm{HI}}$]','h1_c4h1',fraction=len(h1col[indexes_h1_c4])*1./npixels_radius)

plot_ion_ratio(o6col[indexes_n5_o6]-h1col[indexes_n5_o6],n5col[indexes_n5_o6]-h1col[indexes_n5_o6],o6h1lo,o6h1hi,n5h1lo,n5h1hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','o6h1_n5h1',fraction=len(h1col[indexes_n5_o6])*1./npixels_radius)

plot_ion_ratio(h1col[indexes_h1],n5col[indexes_h1]+o6col[indexes_h1],h1lo,h1hi,n5po6lo,n5po6hi,'log[$N_{\mathrm{HI}}$ [cm^-2]]','log[$N_{\mathrm{NV}}+N_{\mathrm{OVI}}$ [cm$^{-2}$+cm$^{-2}$]','h1_n5po6',fraction=len(h1col[indexes_h1])*1./npixels_radius)

plot_ion_ratio(o6col[indexes_h1],n5col[indexes_h1],o6lo,o6hi,n5lo,n5hi,'log[$N_{\mathrm{OVI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}$ (cm$^{-2}$)]','o6_n5',fraction=len(h1col[indexes_h1])*1./npixels_radius)
 
plot_ion_ratio(o6col[indexes_h1],c4col[indexes_h1],o6lo,o6hi,c4lo,c4hi,'log[$N_{\mathrm{OVI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}$ (cm$^{-2}$)]','o6_c4',fraction=len(h1col[indexes_h1])*1./npixels_radius)


#1. N_HI vs N_OVI /  N_HI
#2. N_HI vs N_NV / N_HI
#3. N_NV / N_HI vs N_OVI / N_HI
#4. N_HI vs N_NV + N_OVI




#f_h1_o6n5 = plt.figure(figsize=(6.0,5.0))
#ax_h1_o6n5 = f_h1_o6n5.add_subplot(111)

#ax_h1_o6n5.hexbin(h1col[indexes_n5_o6],o6col[indexes_n5_o6]-n5col[indexes_n5_o6],gridsize=100, C=h1col[indexes_n5_o6]*0.0+1.0,reduce_C_function=np.nansum,zorder=10)
#ax_h1_o6n5.axis([h1lo,h1hi,o6n5lo,o6n5hi])
#ax_h1_o6n5.set_xlabel(r'log $N_{\mathrm{HI}}$',fontsize=16)
#ax_h1_o6n5.set_ylabel(r'log $N_{\mathrm{OVI}}/N_{\mathrm{NV}}$',fontsize=16)

#f_h1_o6n5.subplots_adjust(left=0.135, bottom=0.135)
#f_h1_o6n5.savefig(coldens_basename + '.h1_o6n5.' + dir + '.l' + size + '.png') 

