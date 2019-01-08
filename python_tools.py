import numpy as np
import pickle
import matplotlib.pyplot as plt

def mapload(fname,location): # Loads the pickle file for a particular map, given a filename
    infile = location+'%s.pkl'%(fname)
    pkl_file = open(infile, 'rb')
    myData = pickle.load(pkl_file)
    pkl_file.close()
    
    return myData

def get_binedges(z_store): # Creates an array of bin edges for variable-width histogram plotting, given the bin centres
    bins = []
    bins.append(z_store[0]-(z_store[1]-z_store[0])/2)
    for j in range(len(z_store)):
        if j+1 == len(z_store):
            break
        else:
            bins.append(z_store[j]+(z_store[j+1]-z_store[j])/2)
    bins.append(z_store[-1]+(z_store[-1]-z_store[-2])/2)
    
    return bins
    
def get_bincentres(binedges): # Finds the centre points of a set of bin edges
    bincentres = []
    for i in range(len(binedges)):
        if i+1 == len(binedges):
            break
        else:
            bincentres.append((binedges[i+1]+binedges[i])/2)
    return bincentres
    
def get_binsizes(binedges):
    binsizes = []
    for i in range(len(binedges)):
        if i+1 == len(binedges):
            break
        else:
            binsizes.append(binedges[i+1]-binedges[i])
    return binsizes

def progress_bar(step,tot_steps):
    s = np.float(step)
    t = np.float(tot_steps)
    if s/t < 0.1:
        print '{0}\r'.format('Progress: [          ]'),
    elif 0.1 <= s/t < 0.2:
        print '{0}\r'.format('Progress: [=         ]'),
    elif 0.2 <= s/t < 0.3:
        print '{0}\r'.format('Progress: [==        ]'),
    elif 0.3 <= s/t < 0.4:
        print '{0}\r'.format('Progress: [===       ]'),
    elif 0.4 <= s/t < 0.5:
        print '{0}\r'.format('Progress: [====      ]'),
    elif 0.5 <= s/t < 0.6:
        print '{0}\r'.format('Progress: [=====     ]'),
    elif 0.6 <= s/t < 0.7:
        print '{0}\r'.format('Progress: [======    ]'),
    elif 0.7 <= s/t < 0.8:
        print '{0}\r'.format('Progress: [=======   ]'),
    elif 0.8 <= s/t < 0.9:
        print '{0}\r'.format('Progress: [========  ]'),
    elif 0.9 <= s/t < 1.:
        print '{0}\r'.format('Progress: [========= ]'),

    if step == tot_steps-1:
        print '{0}\r'.format('Progress: [==========]'),
        print
