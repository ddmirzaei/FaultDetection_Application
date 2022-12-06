function StencilsForMesh = StencilSelect(TR,Info,ExTriCnt,NumCells,GhostTri,StencilType,SizeCentSt)
% This function generates  centered, forward sector and backward sector stencils out of mesh triangles

% Inputs: 
% TR, Info: triangulation structure for the mesh
% ExTriCnt: center of triangles in the extended mesh
% NumCells: number of triangles
% GhostTri:  Indices of ghost triangles that share an edge on the mesh boundary
% StencilType: Type of the stencils for triangles
%    StencilType = 1:   7 stencils(1 centerd,3 forward,3 backward) of size 7
%    StencilType = 2:   a central stencil of size SizeCentSt
% SizeCentSt: size of central stencil in the case of StencilType = 2

% Output:
% StencilsForMesh: stencils for all mesh triangles
%
for tri=[1:NumCells,GhostTri]
    if StencilType == 1     % 7 stencils(1 centerd,3 forward,3 backward) of size 7
        [C_Stencil] = CntrStencils(TR,tri,StencilType);
        [F_Stencil] = ForwStencils(TR,Info,tri,7,ExTriCnt);
        [B_Stencil] = BackStencils(TR,Info,tri,7,ExTriCnt);
        Stencils{1,1} = C_Stencil(1,:);
        Stencils{1,2} = F_Stencil(1,:);
        Stencils{1,3} = F_Stencil(2,:);
        Stencils{1,4} = F_Stencil(3,:);
        Stencils{1,5} = B_Stencil(1,:);
        Stencils{1,6} = B_Stencil(2,:);
        Stencils{1,7} = B_Stencil(3,:);            
    elseif StencilType == 2 % Central stencil.
        [Central_Stencil] = CentralStencil(TR,tri,SizeCentSt);
        Stencils{1,1} = Central_Stencil;
    else
        error('StencilType is not implemented ...')
    end   
    StencilsForMesh{1,tri} = Stencils;
end
