function  [U,CPUtime]= ...
          SolveEqWENO(h,Rect,t0,tfinal,dt,GhostLen,eps1,rho,StencilTypeWeno,SizeCentSt,NumGaussP,ScalePar)

% This function computes cell average values at mesh triangles at different time levels
%
% Inputs:
%  h: Mesh size.
%  Rect: triangular domain of the problem
%  t0: Initial time.
%  tfinal: Final time.
%  dt: time difference length.
%  GhostLen: length of ghost cells added on each side of the domain
%  eps1, rho: parameters for the WENO reconstruction
%  StencilTypeWeno: kind of stencils is used for fault triangles
%  StencilTypeCnt: kind of stencils is used for non fault triangles

% Outputs:
%  U:  Cell average values 
%  MeanFaultTri: Mean number of fault triangles in each time step.
%  CPUtime: Cpu time used for the hybrid method.
%
tic  
[TR,Info] = MeshGen(h,Rect);        % triangulation
NumCells = size(TR,1);              % Number of triangles                        
TriCnt = incenter(TR);              % Barycenter of triangles
U = InitialCond(TriCnt);   % Initial condition

% Introducing ghost cells added to the original mesh TR for imposing periodic BC
[GhostCells,ExtendedInfo,ExtendedTR] = BoundMesh(Info,TriCnt,h,Rect,GhostLen); 
U = BoundCond(NumCells,GhostCells,U);        
ExTriCnt = incenter(ExtendedTR);             
[wGauss,xGauss,yGauss,Edge] = GaussianEdges(ExtendedTR,ExtendedInfo,NumGaussP); 
NorEdge = NormalEdges(ExtendedTR,ExtendedInfo,Edge);        
GhostTri = BoundGhostCells(TR,Info,ExtendedTR,Rect);        

U0Gauss = InitialCondOnGauss(TR,GhostTri,xGauss,yGauss);    % Initial condition at the gaussian points
NeighbEdge = NeighborEdge(ExtendedTR,NumCells,Edge);
[TimeCFL,AreaMesh] = CFLcond(ExtendedTR,ExtendedInfo,Edge,NorEdge,NumCells,U0Gauss); % CFL condition and area of triangles

fprintf('\n The selected dt for time discretization should be less than: %d .\n',TimeCFL);
fprintf('\n Construction of stencils for all triangles in the mesh ... \n');

Stencils = ...
   StencilSelect(ExtendedTR,ExtendedInfo,ExTriCnt,NumCells,GhostTri,StencilTypeWeno,SizeCentSt); % 7 stencils of size 7 for mesh triangles

fprintf('\n Computation of interpolation and evaluation matrices based on WENO stencils ... \n');
[MeshInvAP,MeshA,Q,rbfpar] = PHSrbfMesh(TR,Stencils,GhostTri,ExTriCnt,ScalePar); 
MeshAP = PHSrbfMeshEval(TR,xGauss,yGauss,Stencils,ExTriCnt,GhostTri,ScalePar);

% initial setting for the Time discretization
t = t0;                          
FirstIter = 0;                 
NumIter = 0;                   
PlotFig(U,ExTriCnt,'hybrid',0);

% Time and space discretization
while round(t,7) < tfinal       

    % SolValMesh: reconstructing solution values at the all gaussian points based on cell average values U at each time step.
    if FirstIter == 0            
        SolValMesh = U0Gauss; 
        FirstIter = 1;
    else                     
     SolValMesh = RepRBFmesh(MeshAP,MeshInvAP,MeshA,Q,rbfpar,TR,Stencils,GhostTri,U,eps1,rho,NumGaussP,ScalePar);
    end
       
    % advancing in time using Strong Stability Preserving RK method  SSPRK(3) 
    t = t+dt;     % updating time step
    % First step of ssprk(3) 
    for tri = 1:NumCells                       
        ReconstTriOut(:,1) = SolValMesh(NumGaussP:-1:1,NeighbEdge(2,1,tri),NeighbEdge(1,1,tri));
        ReconstTriOut(:,2) = SolValMesh(NumGaussP:-1:1,NeighbEdge(2,2,tri),NeighbEdge(1,2,tri));
        ReconstTriOut(:,3) = SolValMesh(NumGaussP:-1:1,NeighbEdge(2,3,tri),NeighbEdge(1,3,tri));       
        ReconstTriIn = SolValMesh(:,:,tri);
        NormalTri = NorEdge(:,:,tri);     % Normal vectors at the edges of a triangle.
        wTri = wGauss(:,:,tri);           % Gaussian weights on the edges of a triangle.       
        AreaTri = AreaMesh(:,tri);        % Area of a triangle.
        SumFluxTri = FluxTri(wTri,NormalTri,ReconstTriIn,ReconstTriOut,AreaTri,NumGaussP); % Totall flux over edges of a triangle
        Ustep1(tri,1) = U(tri,1)+dt*SumFluxTri;     
    end
    Ustep1 = BoundCond(NumCells,GhostCells,Ustep1); 
    SolValMesh1 = RepRBFmesh(MeshAP,MeshInvAP,MeshA,Q,rbfpar,TR,Stencils,GhostTri,Ustep1,eps1,rho,NumGaussP,ScalePar);   
    
    % Second step of ssprk(3) 
    for tri = 1:NumCells
        ReconstTriOut(:,1) = SolValMesh1(NumGaussP:-1:1,NeighbEdge(2,1,tri),NeighbEdge(1,1,tri));
        ReconstTriOut(:,2) = SolValMesh1(NumGaussP:-1:1,NeighbEdge(2,2,tri),NeighbEdge(1,2,tri));
        ReconstTriOut(:,3) = SolValMesh1(NumGaussP:-1:1,NeighbEdge(2,3,tri),NeighbEdge(1,3,tri));
        ReconstTriIn = SolValMesh1(:,:,tri);
        NormalTri = NorEdge(:,:,tri);
        wTri = wGauss(:,:,tri);
        AreaTri = AreaMesh(:,tri);
        SumFluxTri = FluxTri(wTri,NormalTri,ReconstTriIn,ReconstTriOut,AreaTri,NumGaussP);
        Ustep2(tri,1) = (3/4)*U(tri,1)+(1/4)*Ustep1(tri,1)+(1/4)*dt*SumFluxTri;  
    end    
    Ustep2 = BoundCond(NumCells,GhostCells,Ustep2);
    SolValMesh2 = RepRBFmesh(MeshAP,MeshInvAP,MeshA,Q,rbfpar,TR,Stencils,GhostTri,Ustep2,eps1,rho,NumGaussP,ScalePar);
    
    %  Third step of ssprk(3) 
    for tri=1:NumCells
        ReconstTriOut(:,1) = SolValMesh2(NumGaussP:-1:1,NeighbEdge(2,1,tri),NeighbEdge(1,1,tri));
        ReconstTriOut(:,2) = SolValMesh2(NumGaussP:-1:1,NeighbEdge(2,2,tri),NeighbEdge(1,2,tri));
        ReconstTriOut(:,3) = SolValMesh2(NumGaussP:-1:1,NeighbEdge(2,3,tri),NeighbEdge(1,3,tri));
        ReconstTriIn = SolValMesh2(:,:,tri);
        NormalTri = NorEdge(:,:,tri);
        wTri = wGauss(:,:,tri);
        AreaTri = AreaMesh(:,tri);
        SumFluxTri = FluxTri(wTri,NormalTri,ReconstTriIn,ReconstTriOut,AreaTri,NumGaussP);
        Ustep3(tri,1) = (1/3)*U(tri,1)+(2/3)*Ustep2(tri,1)+(2/3)*dt*SumFluxTri; 
    end
    Ustep3 = BoundCond(NumCells,GhostCells,Ustep3);

    U = Ustep3;     
    disp(['Time level = ',num2str(t)])
    
    % Plotting at some time step
    tt = round(t,7);
    if tt == 0.3 || tt == 0.5 || tt == 1 
        PlotFig(U,ExTriCnt,'hybrid',tt);
    end
end
fprintf('\n');
CPUtime = toc;
