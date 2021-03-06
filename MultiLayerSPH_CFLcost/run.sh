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
f95 -c -fcheck=all -Wno-tabs -g mod_global.f03 mod_integration.f03 mod_sph.f03 mod_boundary.f03 mod_eulerian.f03 main.f03
f95 -o MAIN mod_global.o mod_integration.o mod_sph.o mod_boundary.o mod_eulerian.o main.o
mv MAIN ..
cd ..


echo compiled

./MAIN

mv *.dat output
echo simulazione finita
