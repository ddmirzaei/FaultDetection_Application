# FaultDetection_Application

    By: Davoud Mirzaei and Navid Soodbakhsh, December 2022

------------

This file contains the MATLAB code of section 8  of paper
 
"D. Mirzaei, N. Soodbakhsh, A fault detection method based on partition of unity and kernel approximation, Numerical Algorithms, 2022."

for an application for solving conservation laws by combining the WENO reconstruction FVM and the given fault detection algorithm.

------------
The code works for Burger's equation 

                   u_t + 0.5(u^2)_x = 0 

with periodic boundary conditions on rectangular domain [xl,xr]x[yl,yr].

The user should download all .m files and run 'Start_Run.m' to see the result.  

------------
