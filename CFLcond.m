function [dtCFL,AreaMesh] = CFLcond(TR,Info,Edge,Normals,NumCells,U0)
% This function computes time step that satisfies the CFL condition and areas of triangles

% Inputs:
%  TR,Info: triangulation structure for the mesh.
%  Edge: indices of the ordered nodes on edges 
%  Normals: unit outward normals at edges 
%  NumCells: number of triangles
%  U0: initial condition values at Gaussian points

% Outputs:
%  dtCFL: time step that satisfy cfl condition
%  AreaMesh: areas of all triangles in the mesh
%
sz = size(TR);
NumTri = sz(1);
RadIncircle = [];
for k = 1:NumTri
    Dist(:,:,k) = Info.Nodes(:,Edge(2,:,k))-Info.Nodes(:,Edge(1,:,k));
    for j = 1:3
        Distance(j,k) = norm(Dist(:,j,k),2);
    end
end
for j = 1:NumTri
    LenEdge1 = Distance(1,j); LenEdge2 = Distance(2,j); LenEdge3 = Distance(3,j);
    s = 0.5*(LenEdge1+LenEdge2+LenEdge3);
    RadIncircle(1,j) = (LenEdge1*LenEdge2*LenEdge3)/(4*sqrt(s*(s-LenEdge1)*(s-LenEdge2)*(s-LenEdge3)));
    AreaMesh(1,j) = sqrt(s*(s-LenEdge1)*(s-LenEdge2)*(s-LenEdge3));
end

% computing cfl condition for F(u)=1/2*[u^2;u^2]
DerFlux = @(u) [u;u];
MaxDerFlux = [];
for k = 1:NumCells
    for j = 1:3
       MaxDerFluxEdge(1,j) = max(abs(Normals(:,j,k)'*DerFlux(U0(:,j,k)')));
    end
    MaxDerFlux(1,k) = max(MaxDerFluxEdge);
end
dtCFL = min(RadIncircle(1,1:NumCells)./MaxDerFlux);
