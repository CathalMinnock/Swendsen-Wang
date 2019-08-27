import matplotlib.pyplot as plt
import numpy as np

n=4
q=2
beta = ".5000"
size=8
directory_name="SW_q_" + str(q) + "_size_" + str(size) + "_n_" + str(n)

filename=directory_name + "/" + directory_name + "_beta_" + beta + "_autocorrelation.txt"
x, y  = np.loadtxt(filename, unpack = True)
plt.plot(x, y, 'o', label = "beta = 0.1")

beta = ".2000"
filename=directory_name + "/" + directory_name + "_beta_" + beta + "_autocorrelation.txt"
x, y  = np.loadtxt(filename, unpack = True)
plt.plot(x, y, 'o', label = "beta = 0.2")

beta = ".3000"
filename=directory_name + "/" + directory_name + "_beta_" + beta + "_autocorrelation.txt"
x, y  = np.loadtxt(filename, unpack = True)
plt.plot(x, y, 'o', label = "beta = 0.3")

beta = ".4000"
filename=directory_name + "/" + directory_name + "_beta_" + beta + "_autocorrelation.txt"
x, y  = np.loadtxt(filename, unpack = True)
plt.plot(x, y, 'o', label = "beta = 0.4")

beta = ".5000"
filename=directory_name + "/" + directory_name + "_beta_" + beta + "_autocorrelation.txt"
x, y  = np.loadtxt(filename, unpack = True)
plt.plot(x, y, 'o', label = "beta = 0.5")

plt.xlabel("Number of Swendsen-Wang iterations")
plt.ylabel("Autocorrelation Function")
plt.title("16^3 lattice, Q=2, various inverse temperatures (beta)")
plt.grid(True)
plt.xlim(0,10)
plt.legend(loc = "upper right")
plt.show()
