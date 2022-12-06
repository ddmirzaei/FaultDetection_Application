
% Solving Burgers' equation u_t +u*u_x + u*u_y=0 on [-0.5,0,5]^2 by 
% combination of WENO Finite Volume method with the fault detection algorithm
% See section 8 of the following paper: 
%   D. Mirzaei, N. Soodbakhsh, A fault detection method based on partition
%        of unity and kernel approximation, Numerical Algorithm, 2022. 

% Initial setting
clc
close all
clear all

h = 1/64;                                               % Mesh size.
xl = -0.5; xr = 0.5; yl = -0.5; yr = 0.5;   
Rect.xl = xl; Rect.xr = xr; Rect.yl = yl; Rect.yr = yr; % Rectangular domain of the problem [xl,xr]*[yl,yr]
t0 = 0;                                                 % Initial time.
tfinal = 1;                                             % Final time.
dt = 0.0025;                                            % dt: time discretization step.
GhostLen = 3;                                           % level of ghost cells added on each side of the domain
NumGaussP = 3;                                          % Number of Gaussian integration points on each edge
eps1 = 10^-6;                                           % used in WENO reconstruction indicator
rho = 2;                                                % used WENO reconstruction indicator
ScalePar = h;                                           % scalling parameter for computing the polyharmonic spline rbf
% 
% Hybrid method 
StencilTypeWENO = 1;      % 7 stencils(1 center, 3 forward, 3 backward) of size 7.
StencilTypeCnt = 2;       % A central stencil of size 'size_centered stencil'.
SizeCentSt = 7;                       % Size of central stencil for non fault triangles.
[U_hybrid,MeanFaultTri,CPUtime_hybrid] = ...
    SolveEqHybrid(h,Rect,t0,tfinal,dt,GhostLen,eps1,rho,StencilTypeWENO,StencilTypeCnt,SizeCentSt,NumGaussP,ScalePar);
% 
% WENO method
StencilType = 1;           % 7 stencils(1 center, 3 forward, 3 backward) of size 7.
[U_weno,CPUtime_weno] = SolveEqWENO(h,Rect,t0,tfinal,dt,GhostLen,eps1,rho,StencilTypeWENO,SizeCentSt,NumGaussP,ScalePar);
%
% Results
fprintf('\n Average number of fault triangles: %2f .\n',MeanFaultTri);
fprintf('\n CPU time for the hybrid WENO method: %2f .\n',CPUtime_hybrid);
fprintf('\n CPU time for the WENO method: %2f .\n',CPUtime_weno);
RMSE = norm(U_hybrid-U_weno,2)/sqrt(size(U_hybrid,1));
fprintf('\n RMS error: %d .\n',RMSE);
