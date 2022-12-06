function U = InitialCondOnGauss(TR,IndxGhost,x,y)
% This function computes initial values at the Gaussian points

% Inputs: 
% TR: triangulation structure for the mesh
% IndxGhost: Indices of ghost triangles that share an edge on the mesh boundary
% x, y: x and y coordinates of Gaussian points.

% Outputs:
% U: vector of initial values at the Gaussian points
%
sz = size(TR);
U = [];
for tri = [1:sz(1),IndxGhost]
    for edge = 1:3
        for i = 1:size(x,1)
            U(i,edge,tri) = InitialCond([x(i,edge,tri) y(i,edge,tri)]);
        end
    end   
end
