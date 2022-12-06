function Flux = FluxFunct(u_in,u_out,n)
% This functiom computes the Lax-Friedrich flux at a Gaussian point for F(u) = [u^2 u^2]

% Inputs: 
%  u_in: inside reconstruction at a Gaussian point
%  u_out: outside reconstruction at a Gaussian point
%  n: normal vecotr.

% Outputs:
%  Flux: flux value
%
sigma = max(abs([sum(n)*u_in,sum(n)*u_out]));
if round(u_in,10) == round(u_out,10)
    Flux = 0.5*((dot(([0.5*(u_in.^2);0.5*(u_in.^2)])+([0.5*(u_out.^2);0.5*(u_out.^2)]),n)));
else
    Flux = 0.5*((dot(([0.5*(u_in.^2);0.5*(u_in.^2)])+([0.5*(u_out.^2);0.5*(u_out.^2)]),n))+sigma*(u_in-u_out));
end
