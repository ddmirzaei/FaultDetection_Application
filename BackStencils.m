function BackwardStencils = BackStencils(TR,Info,tri,StSize,ExTriCnt)
% This function generates three backward sector stencils of size 'StSize' for triangle tri

% Inputs:
%  TR, Info: triangulation structure for the mesh.
%  tri: the refrence triangle index
%  StSize: size of the each stenils 
%  ExTriCnt: center of triangles in the extended mesh.

% Outputs:
%  BackwardStencils: three backward stencils
%
NodeTri = Info.Elements(:,tri);
MidNode = [];
LabelMidNode = [];
for i = 1:3         
    if i == 3
        j = 1;
    else
        j = i+1;
    end
    point1 = Info.Nodes(:,NodeTri(i));
    point2 = Info.Nodes(:,NodeTri(j));
    MidNode(:,i) = (point1+point2)/2;
    LabelMidNode(i,1) = i;
end
ThreeKindNode = []; 
for i = 1:3
    b = LabelMidNode;
    b(i) = [];
    ThreeKindNode(i,:)=[LabelMidNode(i),b'];
end
for num = 1:3
    Neighbor = [];
    Bstencil = tri;
    sizeS = 1;
    Rem = StSize - sizeS;
    NextTriNeighbor = tri;
    v1 = MidNode(:,ThreeKindNode(num,1));
    d_v2_v1 = MidNode(:,ThreeKindNode(num,2))-v1;
    d_v3_v1 = MidNode(:,ThreeKindNode(num,3))-v1;
    detA = (d_v2_v1(1)*d_v3_v1(2))-(d_v2_v1(2)*d_v3_v1(1));
    invA = (1/detA)*[d_v3_v1(2) -d_v3_v1(1);-d_v2_v1(2) d_v2_v1(1)];
    kk = 0;
    while sizeS < StSize
        kk = kk+1;
        NeighLevel = neighbors(TR,NextTriNeighbor');
        NeighLevel = reshape(NeighLevel',1,[]);
        NeighLevel(isnan(NeighLevel)) = [];
        z = length(NeighLevel);
        Ab = zeros(1,z);
        for i = z:-1:2
            Ab(i) = sum(NeighLevel(i) == NeighLevel(i-1:-1:1));
        end
        NeighLevel = NeighLevel(Ab == 0);
        NeighLevel = NeighLevel';
        NeighLevel = NeighLevel(~ismember(NeighLevel,Bstencil));       
        AcceptedNeighbor = [];  
        for i = 1:length(NeighLevel)     
            alpha = ExTriCnt(NeighLevel(i),1);
            beta = ExTriCnt(NeighLevel(i),2);
            d_p_v1 = [alpha;beta]-v1;
            c = invA*d_p_v1;
            if ((c(1) >= 0) && (c(2) >= 0))
                AcceptedNeighbor = [AcceptedNeighbor,NeighLevel(i)];
            end
        end       
        LenAccepted = length(AcceptedNeighbor);
        if Rem == LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor];
            Bstencil= [tri,Neighbor];
            break
        elseif Rem < LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor(1:Rem)];
            Bstencil = [tri,Neighbor];
            break
        elseif Rem > LenAccepted
            Neighbor = [Neighbor,AcceptedNeighbor];
            Bstencil = [tri,Neighbor];
            Rem = Rem-LenAccepted;
            sizeS = sizeS+LenAccepted;                                 
            if isempty(AcceptedNeighbor)    
                NextTriNeighbor = NeighLevel';
            else
                NextTriNeighbor = AcceptedNeighbor;
            end
            if isempty(NextTriNeighbor)
                break                
            end
            if kk > 30
                break
            end
        end
    end
       BackwardStencils(num,:) = Bstencil;
end
