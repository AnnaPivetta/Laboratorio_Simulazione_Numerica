import matplotlib
import matplotlib.pyplot as plt
import numpy as np


x, f, error = np.loadtxt("randomwalk.dat", usecols=(0,1,2), delimiter=',', unpack='true')
plt.errorbar(x,f,yerr=error)
g=x**0.5
plt.plot(x,g)
plt.xlabel('N')
plt.ylabel('r')
plt.title('random walk 3D discreto')

plt.show()


file.close();
