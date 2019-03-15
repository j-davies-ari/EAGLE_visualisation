# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

u = np.linspace(0.,1.)

ben_out = np.zeros(len(u))
jon_out = np.zeros(len(u))

ben_out[u<0.5] = (u[u<0.5])**2 * (u[u<0.5] - 1.)
ben_out[u>0.5] = (1. - u[u>0.5])**3

jon_out = (1.-u)**4 * (1.+4.*u)

plt.figure()

plt.plot(u,ben_out,c='k')
plt.plot(u,jon_out,c='r')
plt.show()