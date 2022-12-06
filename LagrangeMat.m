function Lmat = LagrangeMat(Xe,X,xc,rho,op,ScaleOrd,RBFinfo)
% This function computes the Lagrange function for polyharmonic splines (PHS) 
%  kernels by applying the scaling rule 
% Inputs:
%   Xe: evaluation points 
%   X: trial points (centers)
%   xc: the center of stencil
%   rho: the size of stencil
%   op: operator
%   ScaleOrd: the scaling order of op
%   RBFinfo: RBF informations, (type, par, poly order, ...)
% Outputs:
%   Lmat: Lagrange matrix 
%%
numop = length(op);
Lmat = cell(1,numop);
X = (X-xc)/rho;
Xe = (Xe-xc)/rho;
P = PolyMat(X,'1',RBFinfo.PolyOrder);
np = size(P,2);
A = KerMat(X, X, '1',RBFinfo.type,RBFinfo.par);
AP = [A P;P' zeros(np,np)];
[LL,UU,PP] = lu(AP);
for k=1:numop
    R = [KerMat(Xe, X, op{k},RBFinfo.type,RBFinfo.par)'; PolyMat(Xe,op{k},RBFinfo.PolyOrder)'];
    opts.LT = true;
    yy = linsolve(LL,PP*R,opts);
    opt.UT = true;
    Lm = linsolve(UU,yy,opt);
    Lmat{k} = Lm(1:end-np,:)'/rho^ScaleOrd(k);
end
