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

def plot_ion_ratio_contours(xval,yval,xlo,xhi,ylo,yhi,xlab,ylab,plotname,fraction=None):

    f = plt.figure(figsize = (6,6))

    ax = f.add_subplot(111)


    #H, xedges, yedges = np.histogram2d(aa, bb, range=[[293.,1454.0], [464.,1896.0]], bins=(50, 50))
    #extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    #subplots_adjust(bottom=0.15, left=0.15)
    #levels = (1.0e4, 1.0e3, 1.0e2, 2.0e1)
    #cset = contour(H, levels, origin=lower,colors=[black,green,blue,red],linewidths=(1.9, 1.6, 1.5, 1.4),extent=extent)
    #plt.clabel(cset, inline=1, fontsize=10, fmt=%1.0i)
    #for c in cset.collections:
    #    c.set_linestyle(solid)


    dens, xedges, yedges = np.histogram2d(xval, yval, range=[[xlo,xhi],[ylo,yhi]], bins=[50,50])
    cont1 = ax.contour(dens,extent=[xlo,xhi,ylo,yhi],cmap=get_cmap("Blues"))

    ax.axis([xlo,xhi,ylo,yhi])
    ax.set_xlabel(r'%s'%xlab,fontsize=16)
    ax.set_ylabel(r'%s'%ylab,fontsize=16)
    ax.set_title(r'%s'%title1,fontsize=16)

    #div1 = axgrid.make_axes_locatable(ax)
    #cax1 = div1.append_axes("right", size="5%", pad=-0.2, zorder=20)
    #cb1 = f.colorbar(hex1,cax=cax1)
    #cb1.ax.set_xlabel('n$_\mathrm{LOS}$',labelpad=12)

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
    f.savefig(coldens_basename1 + '.' + plotname + '.' + dir + '.l' + size + '.png') 


def plot_ion_ratio_compare_contours(xval1,yval1,xval2,yval2,xlo,xhi,ylo,yhi,xlab,ylab,plotname,fraction=None):

    f = plt.figure(figsize = (6,6))

    ax = f.add_subplot(111)

    dens1, xedges, yedges = np.histogram2d(yval1, xval1, range=[[ylo,yhi],[xlo,xhi]], bins=[50,50])
    #print 'xedges= ', xedges
    #print 'yedges= ', yedges
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    cont1 = ax.contour(dens1,origin='lower',extent=[xlo,xhi,ylo,yhi],cmap=get_cmap("Blues"),linestyles='dashed')
    dens2, xedges, yedges = np.histogram2d(yval2, xval2, range=[[ylo,yhi],[xlo,xhi]], bins=[50,50])
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    cont2 = ax.contour(dens2,origin='lower',extent=[xlo,xhi,ylo,yhi],cmap=get_cmap("Reds"))

    ax.axis([xlo,xhi,ylo,yhi])
    ax.set_xlabel(r'%s'%xlab,fontsize=16)
    ax.set_ylabel(r'%s'%ylab,fontsize=16)

    if fraction != None:
        ax.text((xhi-xlo)*0.86+xlo,(yhi-ylo)*0.92+ylo,r'f=%5.3f'%(float(fraction)),fontsize=14,horizontalalignment='center', verticalalignment='center',zorder=100)

    ax2 = ax.twinx()
    ax2.yaxis.set_major_formatter( NullFormatter() )
    n1,bins1,patches1 = ax2.hist(xval1,20, normed=1,range=[xlo,xhi], histtype='step',zorder=20,color="blue",linestyle='dashed')
    n2,bins2,patches2 = ax2.hist(xval2,20, normed=1,range=[xlo,xhi], histtype='step',zorder=20,color="red")
    ax2.axis([xlo,xhi,0,2])
    ax2.set_ylim(0,5.5)
    
    ax3 = ax.twiny()
    ax3.xaxis.set_major_formatter( NullFormatter() )
    n1,bins1,patches1 = ax3.hist(yval1, 20, normed=1,range=[ylo,yhi], histtype='step',orientation='horizontal',zorder=20,color="blue",label=title1,linestyle='dashed')
    n2,bins2,patches2 = ax3.hist(yval2, 20, normed=1,range=[ylo,yhi], histtype='step',orientation='horizontal',zorder=20,color="red",label=title2)
    ax3.axis([0,2,ylo,yhi])
    ax3.set_xlim(0,5.5)
    ax3.legend(loc="upper left")

    f.subplots_adjust(left=0.135, bottom=0.135)
    f.savefig(outname + '_comp.' + plotname + '.' + dir + '.l' + size + '.png') 
    #f.savefig(coldens_basename1 + '_' +  coldens_basename2 + '.' + plotname + '.' + dir + '.l' + size + '.png') 



def plot_ion_ratio(xval,yval,xlo,xhi,ylo,yhi,xlab,ylab,plotname,fraction=None):

    f = plt.figure(figsize = (6,6))

    ax = f.add_subplot(111)

    hex1 = ax.hexbin(xval,yval,gridsize=100, C=xval*0.0+1.0,vmin=npixels_radius1*0.0,vmax=npixels_radius1*0.01,reduce_C_function=np.nansum,zorder=10,cmap=get_cmap("gist_ncar_r"))#get_cmap("YlGnBu"))
    ax.axis([xlo,xhi,ylo,yhi])
    ax.set_xlabel(r'%s'%xlab,fontsize=16)
    ax.set_ylabel(r'%s'%ylab,fontsize=16)
    ax.set_title(r'%s'%title1,fontsize=16)

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
    f.savefig(coldens_basename1 + '.' + plotname + '.' + dir + '.l' + size + '.png') 

def plot_ion_ratio_compare(xval1,yval1,xval2,yval2,xlo,xhi,ylo,yhi,xlab,ylab,plotname,fraction=None):

    f = plt.figure(figsize = (6,6))

    ax = f.add_subplot(111)

    hex1 = ax.hexbin(xval1,yval1,gridsize=100, C=xval1*0.0+1.0,vmin=npixels_radius1*0.001,vmax=npixels_radius1*0.01,reduce_C_function=np.nansum,zorder=10,cmap=get_cmap("Blues"),alpha=0.5)
    hex2 = ax.hexbin(xval2,yval2,gridsize=100, C=xval1*0.0+1.0,vmin=npixels_radius2*0.001,vmax=npixels_radius2*0.01,reduce_C_function=np.nansum,zorder=10,cmap=get_cmap("Reds"),alpha=0.5)

    ax.axis([xlo,xhi,ylo,yhi])
    ax.set_xlabel(r'%s'%xlab,fontsize=16)
    ax.set_ylabel(r'%s'%ylab,fontsize=16)
    #ax.set_title(r'%s and %s'%(title1,title2),fontsize=16)


    #div1 = axgrid.make_axes_locatable(ax)
    #cax1 = div1.append_axes("right", size="5%", pad=-0.2, zorder=20)
    #cb1 = f.colorbar(hex1,cax=cax1)
    #cb1.set_label('n$_\mathrm{LOS}$')
    #cb1.ax.set_xlabel('n$_\mathrm{LOS}$',labelpad=12)

    if fraction != None:
        ax.text((xhi-xlo)*0.86+xlo,(yhi-ylo)*0.92+ylo,r'f=%5.3f'%(float(fraction)),fontsize=14,horizontalalignment='center', verticalalignment='center',zorder=100)

    ax2 = ax.twinx()
    ax2.yaxis.set_major_formatter( NullFormatter() )
    n1,bins1,patches1 = ax2.hist(xval1,20, normed=1,range=[xlo,xhi], histtype='step',zorder=20,color="blue")
    n2,bins2,patches2 = ax2.hist(xval2,20, normed=1,range=[xlo,xhi], histtype='step',zorder=20,color="red")
    ax2.axis([xlo,xhi,0,2])
    ax2.set_ylim(0,5.5)
    
    ax3 = ax.twiny()
    ax3.xaxis.set_major_formatter( NullFormatter() )
    n1,bins1,patches1 = ax3.hist(yval1, 20, normed=1,range=[ylo,yhi], histtype='step',orientation='horizontal',zorder=20,color="blue",label=title1)
    n2,bins2,patches2 = ax3.hist(yval2, 20, normed=1,range=[ylo,yhi], histtype='step',orientation='horizontal',zorder=20,color="red",label=title2)
    ax3.axis([0,2,ylo,yhi])
    ax3.set_xlim(0,5.5)
    ax3.legend(loc="upper left")

    f.subplots_adjust(left=0.135, bottom=0.135)
    f.savefig(outname + '_comp.' + plotname + '.' + dir + '.l' + size + '.png') 
    #f.savefig(coldens_basename1 + '_' +  coldens_basename2 + '.' + plotname + '.' + dir + '.l' + size + '.png') 



coldens_basename1 = sys.argv[1]
title1 = sys.argv[2] 
coldens_basename2 = sys.argv[3]
title2 = sys.argv[4] 
outname = sys.argv[5]

dir = 'z'
size = '0.6'
radius = 200.

o6filename1 = coldens_basename1 + '.o6.' + dir + '.l' + size + '.dat' 
si3filename1 = coldens_basename1 + '.si3.' + dir + '.l' + size + '.dat' 
h1filename1 = coldens_basename1 + '.h1.' + dir + '.l' + size + '.dat' 
n5filename1 = coldens_basename1 + '.n5.' + dir + '.l' + size + '.dat' 
c4filename1 = coldens_basename1 + '.c4.' + dir + '.l' + size + '.dat' 
si4filename1 = coldens_basename1 + '.si4.' + dir + '.l' + size + '.dat' 

o6filename2 = coldens_basename2 + '.o6.' + dir + '.l' + size + '.dat' 
si3filename2 = coldens_basename2 + '.si3.' + dir + '.l' + size + '.dat' 
h1filename2 = coldens_basename2 + '.h1.' + dir + '.l' + size + '.dat' 
n5filename2 = coldens_basename2 + '.n5.' + dir + '.l' + size + '.dat' 
c4filename2 = coldens_basename2 + '.c4.' + dir + '.l' + size + '.dat' 
si4filename2 = coldens_basename2 + '.si4.' + dir + '.l' + size + '.dat' 


print o6filename1
print o6filename2
 
h1lo = 13.0
h1hi = 18.0
o6lo = 12.0
o6hi = 16.0
n5lo = 11.0
n5hi = 15.0
si3lo = 11.0
si3hi = 15.0
o6si3lo = -1.5
o6si3hi = 1.5
o6n5lo = -1.0
o6n5hi = 2.0
o6h1lo = -3.5
o6h1hi = 1.5
n5h1lo = -4.5
n5h1hi = 0.5
#si3h1lo = -3.5
#si3h1hi = 1.5
n5po6lo = 26.0
n5po6hi = 29.0
n5si3lo = -2
n5si3hi = 5
si3h1lo = -4
si3h1hi = -1
c4si3lo = -1
c4si3hi = 5
si4si3lo = -1.5
si4si3hi = 2.5
c4si4lo = 0.0
c4si4hi = 5.0


o6collim = 13.5
n5collim = 12.5
si3collim = 12.0

o6col1, dist1 =  np.loadtxt(o6filename1, usecols=(2,3), unpack=True)
si3col1, dist1 =  np.loadtxt(si3filename1, usecols=(2,3), unpack=True)
h1col1, dist1 =  np.loadtxt(h1filename1, usecols=(2,3), unpack=True)
n5col1, dist1 =  np.loadtxt(n5filename1, usecols=(2,3), unpack=True)
c4col1, dist1 =  np.loadtxt(c4filename1, usecols=(2,3), unpack=True)
si4col1, dist1 =  np.loadtxt(si4filename1, usecols=(2,3), unpack=True)

o6col2, dist2 =  np.loadtxt(o6filename2, usecols=(2,3), unpack=True)
si3col2, dist2 =  np.loadtxt(si3filename2, usecols=(2,3), unpack=True)
h1col2, dist2 =  np.loadtxt(h1filename2, usecols=(2,3), unpack=True)
n5col2, dist2 =  np.loadtxt(n5filename2, usecols=(2,3), unpack=True)
c4col2, dist2 =  np.loadtxt(c4filename2, usecols=(2,3), unpack=True)
si4col2, dist2 =  np.loadtxt(si4filename2, usecols=(2,3), unpack=True)

indexes_radius1 = np.where(dist1<radius)
npixels_radius1 = len(h1col1[indexes_radius1])

indexes_radius2 = np.where(dist2<radius)
npixels_radius2 = len(h1col1[indexes_radius2])

indexes_si3_o6 = np.where((o6col1>o6collim) & (si3col1>si3collim) & (dist1<radius))
indexes_n5_o6 = np.where((o6col1>o6collim) & (n5col1>n5collim) & (dist1<radius))
indexes_h1_o6 = np.where((dist1<radius))
indexes_h1_n5 = np.where((dist1<radius))
indexes_h1_si3 = np.where((dist1<radius))
indexes_h1 = np.where((dist1<radius))


print "length_radius= ", npixels_radius1
print "length_n5_o6= ", len(o6col1[indexes_n5_o6]), len(o6col1)
print "length_si3_o6= ", len(o6col1[indexes_si3_o6]), len(o6col1)

plot_ion_ratio_compare_contours(o6col1[indexes_radius1]-h1col1[indexes_radius1],n5col1[indexes_radius1]-h1col1[indexes_radius1],o6col2[indexes_radius2]-h1col2[indexes_radius2],n5col2[indexes_radius2]-h1col2[indexes_radius2],o6h1lo,o6h1hi,n5h1lo,n5h1hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','o6h1_n5h1')

plot_ion_ratio_compare_contours(o6col1[indexes_radius1]-h1col1[indexes_radius1],n5col1[indexes_radius1]-si3col1[indexes_radius1],o6col2[indexes_radius2]-h1col2[indexes_radius2],n5col2[indexes_radius2]-si3col2[indexes_radius2],o6h1lo,o6h1hi,n5si3lo,n5si3hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}/N_{\mathrm{SiIII}}$]','o6h1_n5si3')

plot_ion_ratio_compare_contours(n5col1[indexes_radius1]-h1col1[indexes_radius1],si4col1[indexes_radius1]-si3col1[indexes_radius1],n5col2[indexes_radius2]-h1col2[indexes_radius2],si4col2[indexes_radius2]-si3col2[indexes_radius2],n5h1lo,n5h1hi,si4si3lo,si4si3hi,'log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{SiIV}}/N_{\mathrm{SiIII}}$]','n5h1_si4si3')

plot_ion_ratio_compare_contours(n5col1[indexes_radius1]-h1col1[indexes_radius1],c4col1[indexes_radius1]-si3col1[indexes_radius1],n5col2[indexes_radius2]-h1col2[indexes_radius2],c4col2[indexes_radius2]-si3col2[indexes_radius2],n5h1lo,n5h1hi,c4si3lo,c4si3hi,'log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{CIV}}/N_{\mathrm{SiIII}}$]','n5h1_c4si3')


plot_ion_ratio_compare_contours(n5col1[indexes_radius1]-h1col1[indexes_radius1],c4col1[indexes_radius1]-si4col1[indexes_radius1],n5col2[indexes_radius2]-h1col2[indexes_radius2],c4col2[indexes_radius2]-si4col2[indexes_radius2],n5h1lo,n5h1hi,c4si4lo,c4si4hi,'log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{CIV}}/N_{\mathrm{SiIV}}$]','n5h1_c4si4')

plot_ion_ratio_compare_contours(o6col1[indexes_radius1]-n5col1[indexes_radius1],si3col1[indexes_radius1]-h1col1[indexes_radius1],o6col2[indexes_radius2]-n5col2[indexes_radius2],si3col2[indexes_radius2]-h1col2[indexes_radius2],o6n5lo,o6n5hi,si3h1lo,si3h1hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{NV}}$]','log[$N_{\mathrm{SiIII}}/N_{\mathrm{HI}}$]','o6n5_si3h1')

plot_ion_ratio_compare_contours(o6col1[indexes_radius1]-h1col1[indexes_radius1],o6col1[indexes_radius1]-n5col1[indexes_radius1],o6col2[indexes_radius2]-h1col2[indexes_radius2],o6col2[indexes_radius2]-n5col2[indexes_radius2],o6h1lo,o6h1hi,o6n5lo,o6n5hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{OVI}}/N_{\mathrm{NV}}$]','o6h1_o6n5')

plot_ion_ratio_compare_contours(h1col1[indexes_radius1],o6col1[indexes_radius1]-n5col1[indexes_radius1],h1col2[indexes_radius2],o6col2[indexes_radius2]-n5col2[indexes_radius2],h1lo,h1hi,o6n5lo,o6n5hi,'log[$N_{\mathrm{HI}}$]','log[$N_{\mathrm{OVI}}/N_{\mathrm{NV}}$]','h1_o6n5')

plot_ion_ratio_compare_contours(h1col1[indexes_radius1],n5col1[indexes_radius1],h1col2[indexes_radius2],n5col2[indexes_radius2],h1lo,h1hi,n5lo,n5hi,'log[$N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}$]','h1_n5')


plot_ion_ratio_compare_contours(o6col1[indexes_radius1],n5col1[indexes_radius1],o6col2[indexes_radius2],n5col2[indexes_radius2],o6lo,o6hi,n5lo,n5hi,'log[$N_{\mathrm{OVI}}$]','log[$N_{\mathrm{NV}}$]','o6_n5')

plot_ion_ratio_compare_contours(si3col1[indexes_radius1],n5col1[indexes_radius1],si3col2[indexes_radius2],n5col2[indexes_radius2],si3lo,si3hi,n5lo,n5hi,'log[$N_{\mathrm{SiIII}}$]','log[$N_{\mathrm{NV}}$]','si3_n5')

plot_ion_ratio_contours(h1col1[indexes_si3_o6],o6col1[indexes_si3_o6]-si3col1[indexes_si3_o6],h1lo,h1hi,o6si3lo,o6si3hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{SiIII}}$]','h1_o6si3',fraction=len(h1col1[indexes_si3_o6])*1./npixels_radius1)

plot_ion_ratio_contours(o6col1[indexes_radius1]-h1col1[indexes_radius1],n5col1[indexes_radius1]-si3col1[indexes_radius1],o6h1lo,o6h1hi,n5si3lo,n5si3hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}/N_{\mathrm{SiIII}}$]','o6h1_n5si3')


plot_ion_ratio_contours(h1col1[indexes_n5_o6],o6col1[indexes_n5_o6]-n5col1[indexes_n5_o6],h1lo,h1hi,o6n5lo,o6n5hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{NV}}$]','h1_o6n5',fraction=len(h1col1[indexes_n5_o6])*1./npixels_radius1)

plot_ion_ratio_contours(h1col1[indexes_h1_o6],o6col1[indexes_h1_o6]-h1col1[indexes_h1_o6],h1lo,h1hi,o6h1lo,o6h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','h1_o6h1',fraction=len(h1col1[indexes_h1_o6])*1./npixels_radius1)

plot_ion_ratio_contours(h1col1[indexes_h1_n5],n5col1[indexes_h1_n5]-h1col1[indexes_h1_n5],h1lo,h1hi,n5h1lo,n5h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','h1_n5h1',fraction=len(h1col1[indexes_h1_n5])*1./npixels_radius1)

plot_ion_ratio_contours(h1col1[indexes_h1_si3],si3col1[indexes_h1_si3]-h1col1[indexes_h1_si3],h1lo,h1hi,si3h1lo,si3h1hi,'log[$N_{\mathrm{HI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{SiIII}}/N_{\mathrm{HI}}$]','h1_si3h1',fraction=len(h1col1[indexes_h1_si3])*1./npixels_radius1)

plot_ion_ratio_contours(o6col1[indexes_n5_o6]-h1col1[indexes_n5_o6],n5col1[indexes_n5_o6]-h1col1[indexes_n5_o6],o6h1lo,o6h1hi,n5h1lo,n5h1hi,'log[$N_{\mathrm{OVI}}/N_{\mathrm{HI}}$]','log[$N_{\mathrm{NV}}/N_{\mathrm{HI}}$]','o6h1_n5h1',fraction=len(h1col1[indexes_n5_o6])*1./npixels_radius1)

plot_ion_ratio_contours(h1col1[indexes_h1],n5col1[indexes_h1]+o6col1[indexes_h1],h1lo,h1hi,n5po6lo,n5po6hi,'log[$N_{\mathrm{HI}}$ [cm^-2]]','log[$N_{\mathrm{NV}}+N_{\mathrm{OVI}}$ [cm$^{-2}$+cm$^{-2}$]','h1_n5po6',fraction=len(h1col1[indexes_h1])*1./npixels_radius1)

plot_ion_ratio_contours(o6col1[indexes_h1],n5col1[indexes_h1],o6lo,o6hi,n5lo,n5hi,'log[$N_{\mathrm{OVI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}$ (cm$^{-2}$)]','o6_n5',fraction=len(h1col1[indexes_h1])*1./npixels_radius1)
 
plot_ion_ratio_contours(o6col1[indexes_h1],si3col1[indexes_h1],o6lo,o6hi,si3lo,si3hi,'log[$N_{\mathrm{OVI}}$ (cm$^{-2}$)]','log[$N_{\mathrm{NV}}$ (cm$^{-2}$)]','o6_si3',fraction=len(h1col1[indexes_h1])*1./npixels_radius1)

