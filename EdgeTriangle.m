function Edge = EdgeTriangle(TR,Info)
% This function determines indices of the ordered nodes on edges of mesh triangles

% Inputs:
% T R,Info: triangulation structure for the mesh

% Outputs: 
%  Edge: indices of the ordered Nodes corresponding to the edges of mesh triangles.
%
sz = size(TR);
Edge = [];
for i = 1:sz(1)
    t = 1;
    for j = 1:sz(2)-1
        for k = j+1:sz(2)
            Edge(1,t,i) = Info.Elements(j,i);
            Edge(2,t,i) = Info.Elements(k,i);
            t = t+1;
        end
    end   
end
Edge(:,[2,3],:) = Edge(:,[3,2],:);
Edge([1,2],3,:) = Edge([2,1],3,:);
