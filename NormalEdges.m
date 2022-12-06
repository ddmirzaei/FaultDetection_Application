function NorE = NormalEdges(TR,Info,Edge)
% This function computes unit outward normal vector at edges of triangles

% Inputs:
%  TR, Info: triangulation structure for the mesh.
%  Edge: indices of edges of triangulation.

% Outputs:
%  NorE: normal vector at edges 
%
sz = size(TR);
NumTri = sz(1);
for k = 1:NumTri
    for j = 1:3
        NorVec(:,j) = Info.Nodes(:,Edge(2,j,k))-Info.Nodes(:,Edge(1,j,k));
    end
    NorVec([1,2],:) = NorVec([2,1],:);
    NorVec(1,:) = -NorVec(1,:);
    NorVec = -NorVec;
    for j = 1:3
        NorVec(:,j) = NorVec(:,j)/norm(NorVec(:,j));
    end
    NorE(:,:,k) = NorVec;
end