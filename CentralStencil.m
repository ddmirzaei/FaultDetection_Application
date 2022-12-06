function CntStencil = CentralStencil(TR,tri,StSize)
% This function generates a central stencils of size 'StSize' for a reference triangle

% Inputs:
%  TR: triangulation structure for the mesh.
%  tri: The index of the triangle that we want to get its corresponding stencils.
%  size: size of the each stenils (number of triangles in each stencil).

% Outputs:
%  CntStencil: a central stencil of size 'StSize'.
%
Neighbor = [];
CntStencil = tri;
Size = 1;
Rem = StSize-Size;
NextTrNeighbor = tri;
kk = 0;
while Size < StSize
    kk = kk+1;
    NeighbLevel = neighbors(TR,NextTrNeighbor');
    NeighbLevel = reshape(NeighbLevel',1,[]);
    NeighbLevel(isnan(NeighbLevel)) = [];
    Ab=zeros(1,length(NeighbLevel));
    for i = length(NeighbLevel):-1:2
        Ab(i) = sum(NeighbLevel(i) == NeighbLevel(i-1:-1:1));
    end
    NeighbLevel = NeighbLevel(Ab == 0);
    NeighbLevel = NeighbLevel';
    NeighbLevel = NeighbLevel(~ismember(NeighbLevel,CntStencil));
    AcceptedNeighbor = NeighbLevel';
    LenAccepted = length(AcceptedNeighbor);
    if Rem == LenAccepted
        Neighbor = [Neighbor,AcceptedNeighbor];
        CntStencil = [tri,Neighbor];
        break
    elseif Rem < LenAccepted
        Neighbor = [Neighbor,AcceptedNeighbor(1:Rem)];
        CntStencil = [tri,Neighbor];
        break
    elseif Rem > LenAccepted
        Neighbor = [Neighbor,AcceptedNeighbor];
        CntStencil = [tri,Neighbor];
        Rem = Rem-LenAccepted;
        Size = Size+LenAccepted;
        NextTrNeighbor = NeighbLevel';
    end        
end
