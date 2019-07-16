import matplotlib
import matplotlib.pyplot as plt
import numpy as np

x= np.loadtxt("angoli.dat", usecols=(0), unpack='true')
n_bins = 150
n, bins, patches = plt.hist(x, n_bins, range=(0,2*np.pi))

plt.grid(True)

plt.show()
