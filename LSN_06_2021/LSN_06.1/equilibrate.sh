#!/bin/bash

for((T=10;((T<=20));T=T+10))
do 

i=$(dc -e "3k $T 0.1*p")

for((h=20;((h<=20));h=h+2))
do

k=$(dc -e "3k $h 0.01*p")

for((alg=1;((alg<=1));alg=alg+1))
do
		
cat > input.dat << EOF
$i
50
1.0
$k
$alg
50
20
0

  ReadInput >> temp;
  ReadInput >> nspin;
  ReadInput >> J;
  ReadInput >> h;
  ReadInput >> metro;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart_mode;
EOF

make
./Monte_Carlo_ISING_1D.exe

done

done

done



