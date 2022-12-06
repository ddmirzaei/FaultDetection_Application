function Y = ScatPoints2D(bmin,bmax,hcov)
% Domain [bmin,bmax]^2 for the problem.
% Y: covering centers of size (Nc x d)
% hcov: spacing distance between points in Y

[xc,yc] = meshgrid(bmin(1):hcov:bmax(1),bmin(2):hcov:bmax(2));
xc1 = xc(:); xc2 = yc(:); Y = [xc1 xc2];
end % end of function
