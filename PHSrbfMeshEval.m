function AP = PHSrbfMeshEval(TR,x,y,AllStencils,ExTriCnt,GhostTri,ScalePar)
% This function computes the evaluation RBF matrix at all Gaussian points on the mesh 


% Inputs:
% TR: triangulation structure for the mesh.
% x, y: x and y coordinates of Gaussian points
% AllStencils: All stencils on the mesh 
% GhostTri:  Indices of ghost triangles 
% ExTriCnt: Barycenter of tringles in the extended mesh
% ScalePar: scaling parameter 

% Outputs:
% AP: Evaluation matrices for all stencils
%
sz = size(TR);
AP = {};
for tri = [1:sz(1),GhostTri]                                   
    xEpoints = x(:,:,tri);
    yEpoints = y(:,:,tri);
    Stencils = AllStencils{1,tri};
    TriAP = {};
    SizeSt = size(Stencils);
    xEpointsTri = [xEpoints(:,1);xEpoints(:,2);xEpoints(:,3)]; 
    yEpointsTri = [yEpoints(:,1);yEpoints(:,2);yEpoints(:,3)];  
    for i = 1:SizeSt(2)                                            
        Stencil = Stencils{1,i};
        Epoints = [xEpointsTri,yEpointsTri];
        ctrs = ExTriCnt(Stencil,:);
        EM = PHSrbfEval(Epoints,ctrs,ScalePar);           
        TriAP{1,i} = EM;
    end
    AP{1,tri} = TriAP;
end
