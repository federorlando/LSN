import numpy as np
import math

def find_min(vect):
	N = len(vect)
	runner = vect[0]
	index = 0

	for i in range(N):
		if(vect[i] < runner):
			runner = vect[i]
			index = i 
  
	return runner, index

N=441

sigma = np.zeros(N)
mi = np.zeros(N)
fc_value = np.zeros(N)

file = open("./E0_sigma_mi")
line = file.readlines()

for i in range(N):
	data=line[i].split()
	sigma[i] = data[0] 
	mi[i] = data[1]
	fc_value[i] = data[2]

fc_minimum, fc_minimum_index = find_min(fc_value)
sigma_min = sigma[fc_minimum_index]
mi_min = mi[fc_minimum_index]

print('Minimum: ',fc_minimum, ' at line ', fc_minimum_index)
print('Corresponds to sigma = ', sigma_min, ' and mi = ', mi_min)


