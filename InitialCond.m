function u = InitialCond(X)
% Initial condition of the PDE

% Inputs: 
%  X = [x, y]: x and y coordinates
% Output:
%  u: vector of initial values 
%  
c0 = [-0.2,-0.2];
r0 = 0.15;
r = DistMat(X,c0);
func = @(r,r0) exp(r.^2./(r.^2-r0^2));
u = func(r,r0);
u(r >= r0) = 0;

