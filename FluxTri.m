function SumFluxTri = FluxTri(WeightTri,NormalTri,ReconstTriIn,ReconstTriOut,AreaTri,NumGaussP)
% This function computes the total flux on the 3 edges of a triangle

% Inputs: 
%  WeightTri: Gaussian weights on edges of a triangle
%  NormalTri: Normal vectors at edges of a triangle
%  ReconstTriIn: Inside reconstruction at Gaussian points on edges of a triangle
%  ReconstTriOut: Outside reconstruction at Gaussian points on edges of a triangle
%  AreaTri: Area of a triangle
%  NumGaussP: Number of Gauss integration points on each edge

% Outputs:
%  SumFluxTri: Totall flux over edges of a triangle

%%
SumTri = 0;
for edge = 1:3
    WeightEdge = WeightTri(:,edge);
    NormalEdge = NormalTri(:,edge);
    ReconstEdgeIn = ReconstTriIn(:,edge);
    ReconstEdgeOut = ReconstTriOut(:,edge);
    [SumFluxEdge] = FluxEdges(WeightEdge,NormalEdge,ReconstEdgeIn,ReconstEdgeOut,NumGaussP);
    SumTri = SumTri+SumFluxEdge;
end
SumFluxTri = -(1/AreaTri)*SumTri;
