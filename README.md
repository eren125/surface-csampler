# surface-csampler
C++ library for material surface energy sampling

# Includes
Gemmi
https://github.com/project-gemmi/gemmi.git

ThomsonProblemSolver
https://github.com/inlmouse/ThomsonProblemSolver.git

# Compilation
There are four versions of the code: 
1. main\_initial.cpp contains the initial implementation with a sampling at $R_{min}$ and no rejection condition
2. main\_radius.cpp contains the implementation without rejection but with a possibility to change the sampling sphere with the $\lambda$ parameter
3. main\_rejection.cpp contains the implementation with a sampling at $R_{min}$ coupled with a rejection condition
4. main\_final.cpp contains the final implementation with the possibility to change the sampling sphere with the $\lambda$ parameter and the rejection conditithreshold with the $\mu$ parameter

To compile:
c++ -std=c++11 -I./include -O2 src/main_**option**.cpp -o surface_**option**.out

# Usage of the binary
Depending on the option you chose:
1. ./surface\_initial.out   CIF\_FILE FF\_DEF TEMPERATURE CUTOFF SAMPLING\_POINTS\_NUM ADSORBENT
2. ./surface\_radius.out    CIF\_FILE FF\_DEF TEMPERATURE CUTOFF SAMPLING\_POINTS\_NUM ADSORBENT LAMBDA
3. ./surface\_rejection.out CIF\_FILE FF\_DEF TEMPERATURE CUTOFF SAMPLING\_POINTS\_NUM ADSORBENT MU
4. ./surface\_final.out     CIF\_FILE FF\_DEF TEMPERATURE CUTOFF SAMPLING\_POINTS\_NUM ADSORBENT MU LAMBDA

CIF\_FILE: path to the cif file
FF\_DEF: forcefield definition as defined in Raspa2 (see example in test/)
TEMPERATURE: the temperature at which the Boltzmann average will be done
CUTOFF: the cutoff distance of the LJ potential (usually 12 angstrom)
SAMPLING\_POINTS\_NUM: The number of points around a sampling sphere in the simulation
ADSORBENT: for now only atoms (e.g. Xe or Kr). Need an update to take account of multi-atomic adsorbent
MU: the $\mu$ parameter defined in the paper (rejection condition)
LAMBDA: the $\lambda$ parameter as defined in the paper (sampling sphere)

LINK TO THE PAPER
