import numpy as np
import sys
import matplotlib.pyplot as plt


particlefilelist = sys.argv[1]

fparts = file(particlefilelist, 'r')

IDpart, nhpart, Tpart, Zpart, SFRpart, dpart = np.loadtxt(fparts,usecols=(0,1,2,3,4,5),unpack=True)

nproc = 72

for i in range(nproc):
    eaglegasvarname = 'eaglegasvars.' + str(i) +'.txt20'
    #eaglegasvarname = 'test.' + str(i) +'.txt'
    print eaglegasvarname
    
    feagle = file(eaglegasvarname, 'r')
    aex, ID, nh, T, Z, dt, uold, dL, dLeq, ne, nh1, nc4, no1, no2, no3, no4, no5, no6, no7, no8, no9 = np.loadtxt(feagle,usecols=(1,2,3,4,5,6,7,8,9,10,15,24,37,38,39,40,41,42,43,44,45),unpack=True)

    feagle.close
    for j in range(len(IDpart)):
        if(dpart[j] < 1.0):
            print int(IDpart[j])
            indexes = np.where(ID == IDpart[j])
            aexout = aex[indexes]
            nhout = nh[indexes]
            Tout = T[indexes]
            Zout = Z[indexes]
            dtout = dt[indexes]
            uoldout = uold[indexes]
            dLout = dL[indexes]
            dLeqout = dLeq[indexes]
            neout = ne[indexes]
            nh1out = nh1[indexes]
            nc4out = nc4[indexes]
            noout = no1[indexes]+no2[indexes]+no3[indexes]+no4[indexes]+no5[indexes]+no6[indexes]+no7[indexes]+no8[indexes]+no9[indexes]
            no6out = no6[indexes]
            no7out = no7[indexes]
            no8out = no8[indexes]
            no9out = no9[indexes]
            #print aexout
            #print nhout
            #print Tout
            if(i==0):
                fIDpart = file(str(int(IDpart[j])) + '.dat','w')
            else:
                fIDpart = file(str(int(IDpart[j])) + '.dat','a')
            #np.savetxt(IDpart + '.dat', (aexout, nhout, Tout))
            #ndarray.tofile(fIDpart, format='%d % 5.3e %5.3e\n'%(aexout, nhout, Tout))
            for k in range(len(aexout)):
                fIDpart.write('%9.7f % 5.3f %5.3f %6.4f %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e %5.3e\n'%(aexout[k], nhout[k], Tout[k], Zout[k], dtout[k], uoldout[k], dLout[k], dLeqout[k], neout[k], nh1out[k], nc4out[k], noout[k], no6out[k], no7out[k], no8out[k]))
            fIDpart.close
            #np.savetxt('test.out', x, fmt='%1.4e')
            #ndarray.tofile(fid, sep="", format="%s")
