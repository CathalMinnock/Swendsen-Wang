import numpy as np
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]
beta = sys.argv[3]

f = open(input_filename, "r")
f1 = f.readlines()
total = 0
samples_count = 0
for x in f1:
	total += float(x)
	samples_count += 1
f.close()

f = open(input_filename, "r")
f1 = f.readlines()
estimates = np.zeros(samples_count)
i=0
for x in f1:
	estimates[i] = (total - float(x)) / (samples_count - 1)
	i += 1
estimate = np.mean(estimates)
f.close()

variance = 0
for x in range(samples_count):
	variance += (estimates[x] - estimate)**2
variance *= (samples_count - 1.0) / samples_count 
error = np.sqrt(variance)

f=open(output_filename, "a+")
f.write(str(beta) + " " + str(estimate) + " " + str(error) + "\n")
f.close()
