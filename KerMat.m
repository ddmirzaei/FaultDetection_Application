function K = KerMat(X,Y,op,RBFtype,RBFpar)
% This function computes the kernel matrix for different operator 
%  on two point sets X and Y
% Inputs:
%   X: test points of size (m x dim)
%   Y: trial points of size (n x dim)
%   op: operator: '0' for identity, 'x' for first derivative respect to x,
%       etc ...
%   RBFtype: type of RBF
%   RBFpar: the RBF parameter
% Outputs  
%   K: kernel matrix of size (m x n)

dim = size(X,2); % dimension
s = DistMatSqH(X,Y);
switch (op)
    case('1')
        K = Frbf(s,0,RBFtype,RBFpar);
    case('x')
        x = DiffMat(X(:,1),Y(:,1));
        K = x.*Frbf(s,1,RBFtype,RBFpar);
    case('y')
        y = DiffMat(X(:,2),Y(:,2));
        K = y.*Frbf(s,1,RBFtype,RBFpar);
    case('z')
        z = DiffMat(X(:,3),Y(:,3));
        K = z.*Frbf(s,1,RBFtype,RBFpar);
    case('xx')
        x = DiffMat(X(:,1),Y(:,1));
        K = Frbf(s,1,RBFtype,RBFpar)+x.^2.*Frbf(s,2,RBFtype,RBFpar);
    case('yy')
        y = DiffMat(X(:,2),Y(:,2));
        K = Frbf(s,1,RBFtype,RBFpar)+y.^2.*Frbf(s,2,RBFtype,RBFpar);
    case('zz')
        z = DiffMat(X(:,3),Y(:,3));
        K = Frbf(s,1,RBFtype,RBFpar)+z.^2.*Frbf(s,2,RBFtype,RBFpar);
    case ('L')   %\Delta
        K = dim*Frbf(s,1,RBFtype,RBFpar)+2*s.*Frbf(s,2,RBFtype,RBFpar);
    otherwise
        error('this type of KerMat operator (char argument) is not implemented')
end