#!/bin/bash

#Compile, run, plot
gfortran SodShockTube_1D.f90
./a.out 
python3 plot.py

#Move result to the out folder
mkdir -p out
mkdir -p out/data/
mkdir -p out/picture/
mv Lax.csv L_W.csv Mac.csv FVS.csv Roe.csv out/data/
mv Lax.png L_W.png Mac.png FVS.png Roe.png out/picture/

#clean
rm a.out
rm sodshock.mod

echo 'Finish! Result is in the out folder'