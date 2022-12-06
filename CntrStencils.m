function CntStencils = CntrStencils(TR,tri,StencilType)
% This function generates three centered stencils of size 7 for a reference triangle

% Inputs:
%  TR: triangulation structure for the mesh
%  tri: index of the reference triangle
%  StencilType:  type of stencils (1 or 2)

% Outputs
% CntStencils: three centered stencil of size 7.
%
if StencilType == 1
    Neighbor = [];
    Neighbor0 = neighbors(TR,tri);
    for i = 1:3
        Neighbor1 = neighbors(TR,Neighbor0(i));
        Neighbor1(Neighbor1 == tri) = [];
        if i == 3
            Neighbor2 = neighbors(TR,Neighbor0(1));
            Neighbor2(Neighbor2 == tri) = [];
            if any(Neighbor2(1) == Neighbor1)
                Neighbor2 = Neighbor2(2);
            else
                Neighbor2 = Neighbor2(1);
            end
        else
            Neighbor2 = neighbors(TR,Neighbor0(i+1));
            Neighbor2(Neighbor2 == tri) = [];
            if any(Neighbor2(1) == Neighbor1)
                Neighbor2 = Neighbor2(2);
            else
                Neighbor2 = Neighbor2(1);
            end
        end
        Neighbor(i,:) = [Neighbor0,Neighbor1,Neighbor2];
        CntStencils(i,:) = [tri,Neighbor(i,:)];
    end 
end 
