# SodShockTube

## Problem:

The Sod shock tube problem, named after Gary A. Sod, is a common test for the accuracy of computational fluid codes, like Riemann solvers, and was heavily investigated by Sod in 1978. The test consists of a one-dimensional Riemann problem with the following parameters, for left and right states of an ideal gas.

![problem](problem.png 'the Sod shock tube problem')

The time evolution of this problem can be described by solving the Euler equations, which leads to three characteristics, describing the propagation speed of the various regions of the system. Namely the rarefaction wave, the contact discontinuity and the shock discontinuity. If this is solved numerically, one can test against the analytical solution, and get information how well a code captures and resolves shocks and contact discontinuities and reproduce the correct density profile of the rarefaction wave.

Reference: https://en.wikipedia.org/wiki/Sod_shock_tube

## Run: 
    chmod +x run.sh
    ./run.sh
## Result:
![Lax](out/picture/Lax.png 'Lax')
![L-W](out/picture/L_W.png 'L-W')
![Mac](out/picture/Mac.png 'Mac')
![FVS](out/picture/FVS.png 'FVS')
![Roe](out/picture/Roe.png 'Roe')

