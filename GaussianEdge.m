function GaussPW = GaussianEdge(p1,p2,NumGaussP)
% This function computes Gaussian points and weights on edge e = [p1-p2]

% Inputs:
%  p1 = [x1;y1] and p2 = [x2;y2] are two points in 2D

% Outputs:
%  GaussPW: Gaussian points and weights 
%
x1 = p1(1); x2 = p2(1); y1 = p1(2); y2 = p2(2);
[p,w] = GaussLegendre(NumGaussP);
m = (y2-y1)/(x2-x1);
m_inf = (y2-y1)/(round(x2,15)-round(x1,15));
GaussPW = [];
if  ((m_inf == Inf)|| (m_inf == -Inf))
    GaussPW(:,1) = 0.5*abs((y2-y1))*w;
    GaussPW(:,2) = x1;
    GaussPW(:,3) = 0.5*((y2-y1)*p+(y2+y1));
else
    GaussPW(:,1) = 0.5*sqrt(1+(m^2))*abs((x2-x1))*w;
    GaussPW(:,2) = (0.5)*((x2-x1)*p+(x2+x1));
    GaussPW(:,3) = y1+m*(GaussPW(:,2)-x1);
end
