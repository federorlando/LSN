import numpy as np

N = 50
sum_prog = np.zeros(N)
err_prog = np.zeros(N)

file = open("./American_capitals.dat")
line=file.readlines()
for i in range(N):
    data=line[i].split(" ")
    sum_prog[i]=data[0]
    err_prog[i]=data[1]
    print(sum_prog[i], err_prog[i])
