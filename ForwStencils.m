function ForwSt = ForwStencils(TR,Info,tri,StSize,ExTriCnt)
% This function generates three forward sector stencils of size 'StSize' for a reference triangle

% Inputs:
%  TR, Info: triangulation structure for the mesh
%  tri: index of the reference triangle 
%  StSize: size of the each stenils (number of triangles in each stencil)
%  ExTriCnt: barycenters of triangles in the extended mesh

% Outputs:
%  ForwSt: three forward stencils
%
NodeTri = Info.Elements(:,tri);
ThreeKindNode = [];
for i=1:3
    b = NodeTri;
    b(i) = [];
    ThreeKindNode(i,:) = [NodeTri(i),b'];
end
for num = 1:3
    Neighbor = [];
    Fstencil = tri;
    sizeS = 1;
    Rem = StSize-sizeS;
    NextTriNeigh = tri;
    v1 = Info.Nodes(:,ThreeKindNode(num,1));
    d_v2_v1 = Info.Nodes(:,ThreeKindNode(num,2))-v1;
    d_v3_v1 = Info.Nodes(:,ThreeKindNode(num,3))-v1;
    detA = (d_v2_v1(1)*d_v3_v1(2))-(d_v2_v1(2)*d_v3_v1(1));
    invA = (1/detA)*[d_v3_v1(2) -d_v3_v1(1);-d_v2_v1(2) d_v2_v1(1)];
    kk = 0;
    while sizeS < StSize
        kk = kk+1;
        NeighborLevel = neighbors(TR,NextTriNeigh');
        NeighborLevel = reshape(NeighborLevel',1,[]);
        NeighborLevel(isnan(NeighborLevel)) = [];
        Ab=zeros(1,length(NeighborLevel)); 
        for i = length(NeighborLevel):-1:2
            Ab(i) = sum(NeighborLevel(i) == NeighborLevel(i-1:-1:1));
        end
        NeighborLevel = NeighborLevel(Ab == 0);
        NeighborLevel = NeighborLevel';
        NeighborLevel = NeighborLevel(~ismember(NeighborLevel,Fstencil));
        LenNeighbor = length(NeighborLevel);
        AcceptedNeighbor = [];
        for i = 1:LenNeighbor
            alpha = ExTriCnt(NeighborLevel(i),1);
            beta = ExTriCnt(NeighborLevel(i),2);
            d_p_v1 = [alpha;beta]-v1;
            c = invA*d_p_v1;
            if ((c(1) >= 0) && (c(2) >= 0))
                AcceptedNeighbor = [AcceptedNeighbor,NeighborLevel(i)];
            end
        end
        LenAccepted = length(AcceptedNeighbor);
        if Rem == LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor];
            Fstencil = [tri,Neighbor];
            break
        elseif Rem < LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor(1:Rem)];
            Fstencil = [tri,Neighbor];
            break
        elseif Rem > LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor];
            Fstencil = [tri,Neighbor];
            Rem = Rem-LenAccepted;
            sizeS = sizeS+LenAccepted;
            if isempty(AcceptedNeighbor)
                NextTriNeigh = NeighborLevel';
            else
                NextTriNeigh = AcceptedNeighbor;
            end
             if isempty(NextTriNeigh)
                 break
             end
             if kk > 30
                 break
             end
        end
    end
    ForwSt(num,:) = Fstencil;
end
 