#import eagle
import sys
#import coldens
import numpy as np
import matplotlib.pyplot as plt

#def plot_o6_c4_h1_plot

#o6col =  np.loadtxt(o6file, usecols=(2), unpack=True)
#c4col =  np.loadtxt(c4file, usecols=(2), unpack=True)
#h1col =  np.loadtxt(h1file, usecols=(2), unpack=True) 
#ax_.plot





coldens_basename = sys.argv[1]

dir = 'z'
size = '0.6'

o6filename = coldens_basename + '.o6.' + dir + '.l' + size + '.dat' 
c4filename = coldens_basename + '.n5.' + dir + '.l' + size + '.dat' 
h1filename = coldens_basename + '.h1.' + dir + '.l' + size + '.dat' 
n5filename = coldens_basename + '.n5.' + dir + '.l' + size + '.dat' 

print o6filename
print c4filename

h1lo = 13.0
h1hi = 17.0
o6c4lo = -1.5
o6c4hi = 1.5
o6n5lo = -1.5
o6n5hi = 1.5

o6col, dist =  np.loadtxt(o6filename, usecols=(2,3), unpack=True)
c4col, dist =  np.loadtxt(c4filename, usecols=(2,3), unpack=True)
h1col, dist =  np.loadtxt(h1filename, usecols=(2,3), unpack=True)
n5col, dist =  np.loadtxt(n5filename, usecols=(2,3), unpack=True)


indexes_c4_o6 = np.where((o6col>14.0) & (c4col>13.5) & (dist<150.0))

print "length_c4_o6= ", len(o6col[indexes_c4_o6]), len(o6col)

f_h1_o6c4 = plt.figure(figsize=(6.0,5.0))
ax_h1_o6c4 = f_h1_o6c4.add_subplot(111)

ax_h1_o6c4.hexbin(h1col[indexes_c4_o6],o6col[indexes_c4_o6]-c4col[indexes_c4_o6],gridsize=100, C=h1col[indexes_c4_o6]*0.0+1.0,reduce_C_function=np.nansum,zorder=10)
ax_h1_o6c4.axis([h1lo,h1hi,o6c4lo,o6c4hi])
ax_h1_o6c4.set_xlabel(r'log $N_{\mathrm{HI}}$',fontsize=16)
ax_h1_o6c4.set_ylabel(r'log $N_{\mathrm{OVI}}/N_{\mathrm{CIV}}$',fontsize=16)

f_h1_o6c4.subplots_adjust(left=0.135, bottom=0.135)
f_h1_o6c4.savefig(coldens_basename + '.h1_o6c4.' + dir + '.l' + size + '.png') 



indexes_n5_o6 = np.where((o6col>14.0) & (n5col>13.0) & (dist<150.0))

print "length_n5_o6= ", len(o6col[indexes_n5_o6]), len(o6col)

f_h1_o6n5 = plt.figure(figsize=(6.0,5.0))
ax_h1_o6n5 = f_h1_o6n5.add_subplot(111)

ax_h1_o6n5.hexbin(h1col[indexes_n5_o6],o6col[indexes_n5_o6]-n5col[indexes_n5_o6],gridsize=100, C=h1col[indexes_n5_o6]*0.0+1.0,reduce_C_function=np.nansum,zorder=10)
ax_h1_o6n5.axis([h1lo,h1hi,o6n5lo,o6n5hi])
ax_h1_o6n5.set_xlabel(r'log $N_{\mathrm{HI}}$',fontsize=16)
ax_h1_o6n5.set_ylabel(r'log $N_{\mathrm{OVI}}/N_{\mathrm{NV}}$',fontsize=16)

f_h1_o6n5.subplots_adjust(left=0.135, bottom=0.135)
f_h1_o6n5.savefig(coldens_basename + '.h1_o6n5.' + dir + '.l' + size + '.png') 
