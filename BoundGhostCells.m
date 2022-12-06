function GhostTri = BoundGhostCells(TR,Info,ExTR,Rect)
% This function determines ghost triangles that share an edge on the boundary

% Inputs:
%  TR,Info: triangulation structure for the initial mesh
%  ExTR: triangulation structure for the extended mesh (with ghost cells)
%  Rect: domain of the problem

% Outputs:
%  GhostTri:  Indices of ghost triangles that share an edge on the boundary
%
sz = size(TR);
NumTri = sz(1);
[~,LeftNodes] = find(Info.Nodes(1,:) == Rect.xl);
[~,RightNodes] = find(Info.Nodes(1,:) == Rect.xr);
[~,DownNodes] = find(Info.Nodes(2,:) == Rect.yl);
[~,UpNodes] = find(Info.Nodes(2,:) == Rect.yr);
[~,Ind] = sort(Info.Nodes(2,LeftNodes));
LeftNodes = LeftNodes(Ind);
[~,Ind] = sort(Info.Nodes(2,RightNodes));
RightNodes = RightNodes(Ind);
[~,Ind] = sort(Info.Nodes(1,DownNodes));
DownNodes = DownNodes(Ind);
[~,Ind] = sort(Info.Nodes(1,UpNodes));
UpNodes = UpNodes(Ind);
BoundaryNodes = {LeftNodes;RightNodes;DownNodes;UpNodes};
BoundGhostTri = {};
for i = 1:4
    LenNods = length(BoundaryNodes{i,1});
    BDnode = BoundaryNodes{i,1};
    count = 0;
    for j = 1:LenNods-1
        NodeStart = BDnode(j);
        NodeEnd = BDnode(j+1);
        IndTri = edgeAttachments(ExTR,NodeStart,NodeEnd);
        if IndTri{1,1}(1) > NumTri
            count = count+1;
            BoundGhostTri{i,1}(1,count) = IndTri{1,1}(1);
        else
            count = count+1;
            BoundGhostTri{i,1}(1,count) = IndTri{1,1}(2);            
        end
    end
end
LeftGhostTri  = BoundGhostTri{1,1};
RightGhostTri = BoundGhostTri{2,1};
DownGhostTri  = BoundGhostTri{3,1};
UpGhostTri    = BoundGhostTri{4,1};
GhostTri = [LeftGhostTri,RightGhostTri,DownGhostTri,UpGhostTri];
