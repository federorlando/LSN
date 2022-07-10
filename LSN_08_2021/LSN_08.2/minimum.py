import numpy as np
import math

def find_max(matrix):
    N = len(matrix)  #rows
    M = len(matrix[0])  #cols
    runner = matrix[0,0]
    index_i = 0
    index_j = 0

    for i in range(N):
        for j in range(M):
            if(matrix[i,j] > runner):
                runner = matrix[i,j]
                index_i = i
                index_j = j
  
    return runner, index_i, index_j

M=200
L = 0.1
Nbin = 100
binsize = L/Nbin
lowlim_sigma = 0.55
uplim_sigma = 0.65
lowlim_mi = 0.75
uplim_mi = 0.85
hist = np.zeros((Nbin, Nbin))
Sigma = np.zeros(M)
Mi = np.zeros(M)

file2 = open("./points_iter100000")
line = file2.readlines()
for i in range(M):
    data = line[i].split("   ")
    Sigma[i] = data[0]
    Mi[i] = data[1]
file2.close()

for i in range(M):
    for j in range(M):
        value_sigma = Sigma[i]
        value_mi = Mi[j]
        for k in range(Nbin):
            for l in range(Nbin):
                if ((value_sigma >= (lowlim_sigma+k*binsize)) and (value_sigma < (lowlim_sigma+(k+1)*binsize)) and(value_mi >= (lowlim_mi+l*binsize)) and(value_mi < (lowlim_mi+(l+1)*binsize))):
                    hist[k,l]=hist[k,l]+1

max_value, max_sigma, max_mi = find_max(hist)
max_sigma = lowlim_sigma+max_sigma*binsize
max_mi = lowlim_mi+max_mi*binsize
#print(hist)
#print(max_value)
print('sigma = ', max_sigma, ' and mi = ', max_mi)


