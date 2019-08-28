import numpy as np
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]
beta = sys.argv[3]

f = open(input_filename, "r")
f1 = f.readlines()
s4 = 0
s2 = 0
samples_count = 0

for x in f1:
	s4 += float(x)**4
	s2 += float(x)**2
	samples_count += 1
f.close()

s4 /= samples_count
s2 /= samples_count
binder_parameter = 1 - (s4 / (3 * s2**2))

f=open(output_filename, "a+")
f.write(str(beta) + " " + str(binder_parameter) + "\n")
f.close() 
