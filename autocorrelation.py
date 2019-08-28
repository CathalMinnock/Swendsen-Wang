import numpy as np
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]
tau = int(sys.argv[3])

data = np.loadtxt(input_filename, unpack = True)
N = len(data)
mean = np.mean(data)

top = 0
bottom = 0
for i in range(N-tau):
	top += (data[i] - mean) * (data[i+tau] - mean)
	bottom += (data[i] - mean)**2
for i in range(N-tau, N):
		bottom += (data[i] - mean)**2
autocorrelation = top / bottom

f=open(output_filename, "a+")
f.write(str(tau) + " " + str(autocorrelation))
f.close() 
