function SumFluxEdge = FluxEdges(WeightEdge,NormalEdge,ReconstEdgeIn,ReconstEdgeOut,NumGaussP)
% This function computes sum of the flux on gaussian points of a edge of a triangle

% Inputs:
%  WeightEdge: gaussian weights on the edge of a triangle
%  NormalEdge: normal vector
%  ReconstEdgeIn: left reconstruction at the gaussian points on a edge
%  ReconstEdgeOut: right reconstruction at the gaussian points on a edge

% Outputs:
%  SumFluxEdge: sum of the flux on gaussian points of a edge of a triangle
%
SumEdge=0;
for gauss = 1:NumGaussP
    w = WeightEdge(gauss,1);
    u_in = ReconstEdgeIn(gauss,1);
    u_out = ReconstEdgeOut(gauss,1);
    flux = FluxFunct(u_in,u_out,NormalEdge);
    SumEdge = SumEdge+w*flux;
end
SumFluxEdge = SumEdge;
