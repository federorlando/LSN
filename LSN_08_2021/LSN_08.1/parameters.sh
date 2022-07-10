#!/bin/bash

#rm E0_sigma_mi

for ((S=500;((S<=1000));S=S+25))

do

for ((M=500;((M<=1000));M=M+25))

do 

sigma=$(dc -e "3k $S 0.001*p")
mi=$(dc -e "3k $M 0.001*p")

cat > input.dat << EOF
$sigma
$mi

  ReadInput >> sigma;
  ReadInput >> mi;
EOF

make
./main.exe

done

done
