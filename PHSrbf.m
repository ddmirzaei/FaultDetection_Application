function [InvAP,A,Q,rbfpar] = PHSrbf(X,ScalePar)
% This function computes RBF matrix interpolation and the inverse of its
% corresponding saddle point matrix. The RBF is polyharmonic spline (PHS) kernel r^2 log r
% 
% Input:
% X: interpolation points
% ScalePar: scaling parameter 

% Output:
% InvAP: inverse of the saddlepoint matrix [A P;P' 0]
% A: interpolation matrix 
% Q: dimension of the polynomial space 
% rbfpar: r^rbfpar*lor(r) which is used for reconstruction
%
xc = X(1,:);           
X = (X-xc)/ScalePar;
rbfpar = 2; polyorder = 2;
A = KerMat(X,X,'1','tp',rbfpar);
P = PolyMat(X,'1',polyorder);
Q = size(P,2);
AP = [A,P;[P' zeros(Q,Q)]];
InvAP = AP\eye(size(AP));

