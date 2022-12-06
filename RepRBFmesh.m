function SolValMesh = RepRBFmesh(MeshAP,MeshInvAP,MeshA,Q,rbfpar,TR,AllStencils,GhostTri,U,eps1,rho,NumGaussP,ScalePar)
% This function computes solution values at Gaussain edge points based on
% cell averages using WENO method

% Inputs:
% MeshAP: evaluation matrices on the mesh
% MeshInvAP: inverse of the interpolation matrices on the mesh
% MeshA: interpolation matrices on the mesh 
% Q: dimension of the polynomial space
% rbfpar: r^rbfpar*lor(r) 
% TR: triangulation structure for the initial mesh.
% AllStencils: all stencils on the mesh
% GhostTri:  Indices of ghost triangles that share an edge on the mesh boundary.
% U: cell average values on the mesh
% eps1, rho: parameters for the weno reconstruction indicator
% ScalePar: scaling parameter for PHS kernel
% NumGaussP: number of Gaussian integration points on each edge of triangles

% Outputs: 
% SolValMesh: reconstructed solution based on cell average values
%
sz = size(TR);
for tri = [1:sz(1),GhostTri]                    
    TriInvAP = MeshInvAP{1,tri};
    TriA = MeshA{1,tri};
    Stencils = AllStencils{1,tri};
    TriAP = MeshAP{1,tri};
    SizeSt = size(Stencils);
    SolValSt = [];
    OscillationSt = [];
    for i = 1:SizeSt(2)                      
        InvAP = TriInvAP{1,i};
        A = TriA{1,i};
        AP = TriAP{1,i};
        uval = U([Stencils{1,i}],:);
        
        Coef = InvAP*[uval;zeros(Q,1)];           
        CoeRBF = Coef(1:(end-Q),:);
        Oscillation = CoeRBF'*A*CoeRBF;         
        Oscillation = Oscillation*(ScalePar^(-rbfpar));  
        Solval = AP*Coef;                                         
        SolValSt(i,:) = Solval;                  
        OscillationSt(i,:) = Oscillation;    
    end
    w0 = (eps1+OscillationSt).^(-rho);    
    w = w0/sum(w0);
    SolValEdge1 = w'*SolValSt(:,1:NumGaussP);      
    SolValEdge2 = w'*SolValSt(:,NumGaussP+1:2*NumGaussP);
    SolValEdge3 = w'*SolValSt(:,2*NumGaussP+1:3*NumGaussP);
    SolVal = [SolValEdge1',SolValEdge2',SolValEdge3'];    
    SolValMesh(:,:,tri) = SolVal;      
end

