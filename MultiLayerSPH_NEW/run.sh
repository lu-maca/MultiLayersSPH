#!/bin/bash

rm *.dat
cd output
rm *.dat
cd ..
cd src
rm *.o
rm *.mod
cd ..

cd src
OMP_NUM_THREADS=4;
export OMP_NUM_THREADS
f95 -c  mod_global.F03 mod_integration.F03 mod_sph.F03 mod_boundary.F03 mod_eulerian.F03 main.F03
f95 -o MAIN mod_global.o mod_integration.o mod_sph.o mod_boundary.o mod_eulerian.o main.o -fopenmp
mv MAIN ..
cd ..


echo compiled

./MAIN

mv *.dat output
echo simulazione finita
