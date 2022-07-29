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

