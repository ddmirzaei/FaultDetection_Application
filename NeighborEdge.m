function [NeighbEdge] = NeighborEdge(TR,NumCells,Edge)
% Looping on triangulation TR for edges 'e1' of the triangle 'T1' this function determines
%  neighboring tringle 'T2' with edge 'e2' so that 'e1' and 'e2' are the same

% Inputs:
%  TR: triangulation structure for the mesh
%  NumCells: number of triangles in the mesh
%  Edge: indices of edges in triangulation TR

% Output
%  NeighbEdge: indices of neighboring edges 
%
NeighbEdge = [];
for tri = 1:NumCells
    for j = 1:3
        StartNode = Edge(1,j,tri);
        EndNode = Edge(2,j,tri);
        IndTri = edgeAttachments(TR,StartNode,EndNode);
        if IndTri{1}(1) == tri
            NeighbEdge(1,j,tri) = IndTri{1}(2);
        else
            NeighbEdge(1,j,tri) = IndTri{1}(1);
        end
        count = 0;
        for i = 1:3
            if (StartNode == Edge(1,i, NeighbEdge(1,j,tri))) && (EndNode == Edge(2,i, NeighbEdge(1,j,tri)))
                count = 1;
                NeighbEdge(2,j,tri) = i;
                break;
            end
        end
        if count == 0
            for i = 1:3
                if (StartNode == Edge(2,i, NeighbEdge(1,j,tri))) && (EndNode == Edge(1,i, NeighbEdge(1,j,tri)))
                    NeighbEdge(2,j,tri) = i;
                    break
                end
            end
        end  
    end
end