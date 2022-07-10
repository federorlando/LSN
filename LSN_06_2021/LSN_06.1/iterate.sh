#!/bin/bash

for((T=11;((T<=20));T=T+1))
do 

i=$(dc -e "3k $T 0.1*p")

for((h=0;((h<=2));h=h+2))
do

k=$(dc -e "3k $h 0.01*p")

for((alg=0;((alg<=1));alg=alg+1))
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

echo Running with T = $i, h = $k, algorithm = $alg
echo Preliminary run to equilibrate: length of run is 1000 steps
echo

make
./Monte_Carlo_ISING_1D.exe

if ((alg==0)) && ((h==0)); then rm ./T_$i/output.chi.0 ./T_$i/output.cv.0 ./T_$i/output.ene.0 ./T_$i/output.mag.0; fi

if ((alg==0)) && ((h==2)); then rm ./T_$i/output.chi.h_ext.0 ./T_$i/output.cv.h_ext.0 ./T_$i/output.ene.h_ext.0 ./T_$i/output.mag.h_ext.0; fi

if ((alg==1)) && ((h==0)); then	rm ./T_$i/output.chi.1 ./T_$i/output.cv.1 ./T_$i/output.ene.1 ./T_$i/output.mag.1; fi

if ((alg==1)) && ((h==2)); then	rm ./T_$i/output.chi.h_ext.1 ./T_$i/output.cv.h_ext.1 ./T_$i/output.ene.h_ext.1 ./T_$i/output.mag.h_ext.1; fi

cat > input.dat << EOF
$i
50
1.0
$k
$alg
50
100000
1

  ReadInput >> temp;
  ReadInput >> nspin;
  ReadInput >> J;
  ReadInput >> h;
  ReadInput >> metro;
  ReadInput >> nblk;
  ReadInput >> nstep;
  ReadInput >> restart_mode;
EOF

echo Now long run starts, with restart mode 1, i.e. from previous config

make
./Monte_Carlo_ISING_1D.exe

echo Simulation completed
echo -------------------------------------------------------------
echo -------------------------------------------------------------
echo -------------------------------------------------------------

done

done

done



