function [InvAP,A,Q,rbfpar] = PHSrbfMesh(TR,Stencils,GhostTri,ExTriCnt,ScalePar)
% This function gethers all RBF matrices on the mesh

% Inputs:
% TR: triangulation structure for the mesh
% Stencils: stencils for all mesh triangles
% GhostTri:  indices of ghost triangles that share an edge on the mesh boundary
% ExTriCnt: Barycenter of tringles in the extended mesh
% ScalePar: scaling parameter 

% Outputs:
% InvAP: inverse of the saddlepoint matrix [A P;P' 0]
% A: interpolation matrix 
% Q: dimension of the polynomial space 
% rbdpar: r^rbfpar*lor(r) which is used for reconstruction.
%
sz = size(TR);                                          
for tri = [1:sz(1),GhostTri]                                             
    St = Stencils{1,tri};
    SizeSt = size(St); NumSt = SizeSt(2);
    TriInvAP = [];
    TriA = [];
    for i = 1:NumSt                         
        Stencil = St{1,i};
        dsites = ExTriCnt(Stencil,:);                  
        [InvAPTri0,A0,Q,rbfpar] = PHSrbf(dsites,ScalePar); 
        TriInvAP{1,i} = InvAPTri0;
        TriA{1,i} = A0;  
    end
    InvAP{1,tri} = TriInvAP;
    A{1,tri} = TriA;
end