function [w,x,y,Edge] = GaussianEdges(TR,Info,NumGaussP)
% This function Produces Gaussian points on each edge of triangles 

% Inputs:
%  TR,Info: triangulation structure for the mesh

% Ouputs: 
%  x:  x_coordinate of Gaussian points
%  y:  y_coordinate of Gaussian points
%  w:  gaussian weights.
%  Edge: indices of points on edges
%

sz = size(TR);
NumTri = sz(1);
NumEdges = sz(2);
Edge = EdgeTriangle(TR,Info);
w = []; x = []; y = [];
for i = 1:NumTri
    for j=1:NumEdges
        P1 = Info.Nodes(1:2,Edge(1,j,i));
        P2 = Info.Nodes(1:2,Edge(2,j,i));
        Gauss = GaussianEdge(P1,P2,NumGaussP);
        w(:,j,i)=Gauss(:,1);
        x(:,j,i)=Gauss(:,2);
        y(:,j,i)=Gauss(:,3);
    end
end