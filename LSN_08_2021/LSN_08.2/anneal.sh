#!/bin/bash

#rm E0_sigma_mi

for ((S=20000;((S<=100000));S=S+10000))

do

#sigma=$(dc -e "3k $S 0.001*p")
#mi=$(dc -e "3k $M 0.001*p")

cat > input.dat << EOF
10
$S
100
0.001
100
10000
1
0.625
0.8

  ReadInput >> temp;
  ReadInput >> niter;
  ReadInput >> vol;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart_mode;
  ReadInput >> sigma;
  ReadInput >> mi;
EOF

make
./Monte_Carlo_NVT.exe

done
