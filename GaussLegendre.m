function [x,w]=GaussLegendre(n)
% Integration points and weights for the Gauss-Legendre formula

if n <= 1, error('n must be >1'); end
[a,b] = coeflegendre(n);
JacM = diag(a)+diag(sqrt(b(2:n)),1)+diag(sqrt(b(2:n)),-1);
[w,x] = eig(JacM); x = diag(x); w = w(1,:)'.^2*b(1);
[x,ind] = sort(x); w = w(ind);
end

function [a,b]=coeflegendre(n)
%Coefficients of Legendre polynomials.
if n <= 1, error('n must be >1'); end
a = zeros(n,1); b = a; b(1) = 2;
k = 2:n; b(k) = 1./(4-1./(k-1).^2);
end